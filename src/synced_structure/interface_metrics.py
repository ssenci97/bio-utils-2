"""Pairwise interface metrics for all chain pairs: buried surface area, contacts, distances, energies, and partners, keyed by standard and author chain IDs."""

###################################### IMPORTS AND CONSTANTS

import itertools

import numpy as np
import pandas as pd
from Bio.PDB.Model import Model
from Bio.PDB.SASA import ShrakeRupley
from tqdm import tqdm

from .constants import (
    CONTACTS_CUTOFF_NM, HYDRO, INTERFACE_COMPLEX_COLUMNS, INTERFACE_RESIDUE_COLUMNS, POLYMER_MOIETIES, TIEN, TO1,
)
from .residues import is_heavy, residue_moiety
from .utils import log_info


def split_partners(text):
    """Split a semicolon-separated partner string into a list of non-empty labels."""
    return [p for p in str(text).split(";") if p]

###################################### SASA AND ATOM GEOMETRY


class InterfaceMixin:
    @staticmethod
    def _heavy_atoms(residues):
        """Return heavy-atom coordinates and the residue index owning each atom.

        Args:
            residues: List of Biopython Residues.
        """
        coords, owners = [], []
        for index, residue in enumerate(residues):
            for atom in residue:
                if is_heavy(atom):
                    coords.append(atom.coord)
                    owners.append(index)
        return np.array(coords, float).reshape(-1, 3), np.array(owners, int)

    @staticmethod
    def _residue_sasa(chains):
        """Return {(chain_id, residue_id): SASA in A^2} for polymer residues, computed on copies of the chains.

        Args:
            chains: List of Biopython Chains; ligands, ions, and water are excluded from the calculation.
        """
        model = Model(0)
        for chain in chains:
            copied = chain.copy()
            for residue in [r for r in copied if residue_moiety(r) not in POLYMER_MOIETIES]:
                copied.detach_child(residue.id)
            model.add(copied)
        ShrakeRupley().compute(model, level="R")
        return {(chain.id, residue.id): residue.sasa for chain in model for residue in chain}

    @staticmethod
    def _atom_contacts(xa, xb, cutoff, score_cutoff, block):
        """Return contact count, score, minimum distance, and per-atom contact masks for two atom sets.

        Args:
            xa: Coordinates of atoms in the first set.
            xb: Coordinates of atoms in the second set.
            cutoff: Contact distance cutoff in A.
            score_cutoff: Distance cutoff in A for the inverse-distance score.
            block: Number of rows processed per block.
        """
        n, score, dmin = 0, 0.0, np.inf
        ha, hb = np.zeros(len(xa), bool), np.zeros(len(xb), bool)
        for start in range(0, len(xa), block):
            d = np.linalg.norm(xa[start:start + block, None] - xb[None], axis=-1)
            m = d <= cutoff
            n += int(m.sum())
            w = (d <= score_cutoff) & (d > 0)
            score += float((1.0 / d[w]).sum())
            dmin = min(dmin, float(d.min()))
            ha[start:start + block] |= m.any(axis=1)
            hb |= m.any(axis=0)
        return n, score, dmin, ha, hb

    @staticmethod
    def _residue_min_distances(xa, oa, xb, ob, na, nb, block):
        """Return an (na, nb) matrix of minimum heavy-atom distances between residue pairs.

        Args:
            xa: Coordinates of atoms in the first set.
            oa: Residue index owning each atom in the first set.
            xb: Coordinates of atoms in the second set.
            ob: Residue index owning each atom in the second set.
            na: Number of residues in the first set.
            nb: Number of residues in the second set.
            block: Number of rows processed per block.
        """
        out = np.full((na, nb), np.inf, np.float32)
        for start in range(0, len(xa), block):
            d = np.linalg.norm(xa[start:start + block, None] - xb[None], axis=-1).astype(np.float32)
            np.minimum.at(out, (np.repeat(oa[start:start + block], len(ob)), np.tile(ob, d.shape[0])), d.ravel())
        return out

    @staticmethod
    def _hydrophobicity(residue):
        """Return the hydropathy value of a residue, or 0.0 if unknown.

        Args:
            residue: Biopython Residue.
        """
        return HYDRO.get(TO1.get(residue.get_resname().upper(), "X"), 0.0)

###################################### PER-CHAIN GEOMETRY

    def _chain_geometry(self, chain_key):
        """Return cached per-chain geometry: polymer residues, heavy atoms, isolated SASA, and hydropathy.

        Args:
            chain_key: Chain key as held by the structure.
        """
        key = ("geometry", chain_key)
        if key not in self._interface_cache:
            chain = self.chains[chain_key]["chain"]
            residues = self.chains[chain_key]["residues"]
            xyz, owner = self._heavy_atoms(residues)
            sasa = self._residue_sasa([chain])
            self._interface_cache[key] = {
                "id": chain.id,
                "chain_id": self.chains[chain_key]["chain_id"],
                "residues": residues,
                "xyz": xyz,
                "owner": owner,
                "hydro": np.array([self._hydrophobicity(r) for r in residues]),
                "sasa": sasa,
                "total_sasa": float(sum(sasa.values())),
            }
        return self._interface_cache[key]

###################################### INTERFACE COMPUTATION

    def _compute_interface(self, chain_a, chain_b):
        """Compute complex-level and residue-level interface metrics for two chains.

        Args:
            chain_a: First chain key.
            chain_b: Second chain key.
        """
        cutoff = self.config["interface_contact_cutoff_a"]
        score_cutoff = self.config["interface_score_cutoff_a"]
        block = self.config["interface_block_size"]
        ga, gb = self._chain_geometry(chain_a), self._chain_geometry(chain_b)
        if not len(ga["xyz"]) or not len(gb["xyz"]):
            raise ValueError(f"No heavy atoms in chain {chain_a} or {chain_b}.")
        cpx = self._residue_sasa([self.chains[chain_a]["chain"], self.chains[chain_b]["chain"]])
        initial = ga["total_sasa"] + gb["total_sasa"]
        dasa = initial - sum(cpx.values())
        n, score, dmin_c, ha, hb = self._atom_contacts(ga["xyz"], gb["xyz"], cutoff, score_cutoff, block)
        dmin = self._residue_min_distances(ga["xyz"], ga["owner"], gb["xyz"], gb["owner"],
                                           len(ga["residues"]), len(gb["residues"]), block)
        contact = dmin <= cutoff
        mask = dmin <= score_cutoff
        with np.errstate(invalid="ignore", divide="ignore"):
            energy = np.where(mask, -0.5 * (ga["hydro"][:, None] + gb["hydro"][None]) * (1 - dmin / score_cutoff), 0.0)
        pair = f"{ga['chain_id']}-{gb['chain_id']}"
        complex_row = {
            "structure": self.structure_id,
            "chain_a": ga["chain_id"], "auth_chain_a": self.chains[chain_a]["auth_chain_id"],
            "chain_b": gb["chain_id"], "auth_chain_b": self.chains[chain_b]["auth_chain_id"],
            "len_res_a": len(ga["residues"]), "len_res_b": len(gb["residues"]),
            "len_atoms_a": len(ga["xyz"]), "len_atoms_b": len(gb["xyz"]),
            "dASA_A2": round(dasa, 2),
            "initial_sasa_A2": round(initial, 2),
            "frac_dASA": round(dasa / initial, 4) if initial > 0 else 0.0,
            "n_contacts": n,
            "iface_res_a": int(np.unique(ga["owner"][ha]).size),
            "iface_res_b": int(np.unique(gb["owner"][hb]).size),
            "min_dist_A": round(dmin_c, 3),
            "contact_score": round(score, 4),
        }
        sides = (
            (chain_a, ga, gb, energy.sum(axis=1), mask.sum(axis=1), contact.sum(axis=1), dmin.min(axis=1), contact),
            (chain_b, gb, ga, energy.sum(axis=0), mask.sum(axis=0), contact.sum(axis=0), dmin.min(axis=0), contact.T),
        )
        rows = []
        for cid, g, other, en_sum, nscore, npart, dm, pairs in sides:
            for k, residue in enumerate(g["residues"]):
                if not npart[k]:
                    continue
                key = (g["id"], residue.id)
                dsasa = float(g["sasa"][key] - cpx[key])
                maxa = TIEN.get(residue.get_resname().upper())
                en = float(en_sum[k])
                partners = [other["residues"][j] for j in np.flatnonzero(pairs[k])]
                labels = ";".join(f"{other['chain_id']}{p.id[1]} {p.get_resname()}" for p in partners)
                rows.append({
                    "structure": self.structure_id,
                    "pair": pair,
                    "chain": g["chain_id"],
                    "auth_chain": self.chains[cid]["auth_chain_id"],
                    "resid": residue.id[1],
                    "resname": residue.get_resname(),
                    "dSASA_A2": round(dsasa, 2),
                    "isolated_sasa_A2": round(float(g["sasa"][key]), 2),
                    "rel_dSASA": round(dsasa / maxa, 4) if maxa else float("nan"),
                    "n_partners": int(npart[k]),
                    "min_dist_A": round(float(dm[k]), 3),
                    "contact_energy": round(en, 4),
                    "energy_per_partner": round(en / nscore[k], 4),
                    "partners": labels,
                })
        return {"complex": complex_row, "residues": pd.DataFrame(rows, columns=INTERFACE_RESIDUE_COLUMNS)}

###################################### INTERFACE ACCESSORS

    def interface(self, chain_a="A", chain_b="B"):
        """Return {'complex': dict, 'residues': DataFrame} for two chains, using the cache when possible.

        Args:
            chain_a: First chain identifier, standard or author.
            chain_b: Second chain identifier, standard or author.
        """
        chain_a, chain_b = self._require_chain(chain_a), self._require_chain(chain_b)
        if chain_a == chain_b:
            raise ValueError("chain_a and chain_b must be different chains.")
        key = (chain_a, chain_b, self.config["interface_contact_cutoff_a"],
               self.config["interface_score_cutoff_a"], self.config["interface_block_size"])
        if key not in self._interface_cache:
            self._interface_cache[key] = self._compute_interface(chain_a, chain_b)
        cached = self._interface_cache[key]
        return {"complex": dict(cached["complex"]), "residues": cached["residues"].copy()}

    def interface_all(self, chain_pairs=None):
        """Return unified complex and residue tables for every chain pair, skipping pairs that cannot be computed.

        Args:
            chain_pairs: Iterable of (chain_a, chain_b) tuples, or None for all pairs.
        """
        pairs = list(chain_pairs or itertools.combinations(self.chain_ids, 2))
        complexes, residues = [], []
        for chain_a, chain_b in tqdm(pairs, desc="pairs", unit="pair", leave=False):
            try:
                result = self.interface(chain_a, chain_b)
            except ValueError as ex:
                log_info(f"{self.structure_id}: pair {chain_a}-{chain_b} skipped, {ex}")
                continue
            complexes.append(result["complex"])
            if not result["residues"].empty:
                residues.append(result["residues"])
        complex_frame = pd.DataFrame(complexes, columns=INTERFACE_COMPLEX_COLUMNS)
        residue_frame = (pd.concat(residues, ignore_index=True) if residues
                         else pd.DataFrame(columns=INTERFACE_RESIDUE_COLUMNS))
        return {"complex": complex_frame, "residues": residue_frame}

    def chain_partner_interfaces(self, chain_id, cutoff=CONTACTS_CUTOFF_NM):
        """Return {'{structure_id}_{chain}_{partner}': interface dict} for every contact partner of one chain.

        Args:
            chain_id: Chain identifier, standard or author.
            cutoff: Contact cutoff in nm, passed to chain_contact_partners.
        """
        chain_id = self._require_chain(chain_id)
        out = {}
        for partner in self.chain_contact_partners(cutoff=cutoff).get(chain_id, []):
            try:
                out[f"{self.structure_id}_{chain_id}_{partner}"] = self.interface(chain_id, partner)
            except ValueError as ex:
                log_info(f"{self.structure_id}: pair {chain_id}-{partner} skipped, {ex}")
        return out

    def chain_partner_frame(self, chain_id, cutoff=CONTACTS_CUTOFF_NM, level="complex"):
        """Return a DataFrame with one row per interface pair of a chain, or per residue of those pairs.

        Args:
            chain_id: Chain identifier, standard or author.
            cutoff: Contact cutoff in nm, passed to chain_partner_interfaces.
            level: 'complex' for one row per pair, or 'residues' for the residue rows of every pair.
        """
        if level not in ("complex", "residues"):
            raise ValueError("level must be 'complex' or 'residues'.")
        interfaces = self.chain_partner_interfaces(chain_id, cutoff=cutoff)
        if level == "complex":
            return pd.DataFrame([v["complex"] for v in interfaces.values()], columns=INTERFACE_COMPLEX_COLUMNS)
        frames = [v["residues"] for v in interfaces.values() if not v["residues"].empty]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=INTERFACE_RESIDUE_COLUMNS)

    def interface_complex(self, chain_a="A", chain_b="B"):
        """Return complex-level interface metrics as a dict.

        Args:
            chain_a: First chain identifier.
            chain_b: Second chain identifier.
        """
        return self.interface(chain_a, chain_b)["complex"]

    def interface_residues(self, chain_a="A", chain_b="B"):
        """Return residue-level interface metrics as a DataFrame.

        Args:
            chain_a: First chain identifier.
            chain_b: Second chain identifier.
        """
        return self.interface(chain_a, chain_b)["residues"]

    def buried_area(self, chain_a="A", chain_b="B"):
        """Return the buried surface area in A^2 between two chains."""
        return self.interface_complex(chain_a, chain_b)["dASA_A2"]

    def initial_sasa(self, chain_a="A", chain_b="B"):
        """Return the sum of isolated SASA of both chains in A^2."""
        return self.interface_complex(chain_a, chain_b)["initial_sasa_A2"]

    def buried_fraction(self, chain_a="A", chain_b="B"):
        """Return buried area divided by initial SASA."""
        return self.interface_complex(chain_a, chain_b)["frac_dASA"]

    def contact_count(self, chain_a="A", chain_b="B"):
        """Return the number of heavy-atom pairs within the contact cutoff."""
        return self.interface_complex(chain_a, chain_b)["n_contacts"]

    def min_distance(self, chain_a="A", chain_b="B"):
        """Return the minimum heavy-atom distance between two chains in A."""
        return self.interface_complex(chain_a, chain_b)["min_dist_A"]

    def contact_score(self, chain_a="A", chain_b="B"):
        """Return the sum of inverse heavy-atom distances within the score cutoff."""
        return self.interface_complex(chain_a, chain_b)["contact_score"]

    def partners(self, chain_a="A", chain_b="B"):
        """Return a dict mapping interface residue labels to their partner labels for one chain pair.

        Args:
            chain_a: First chain identifier.
            chain_b: Second chain identifier.
        """
        frame = self.interface_residues(chain_a, chain_b)
        return {
            f"{row.chain}{row.resid} {row.resname}": split_partners(row.partners)
            for row in frame.itertuples(index=False)
        }