#!/usr/bin/env python
"""Pairwise interface metrics: buried surface area, heavy-atom contacts, distances, energies, and partners."""

###################################### IMPORTS AND CONSTANTS

import numpy as np
import pandas as pd
from Bio.PDB.Model import Model
from Bio.PDB.SASA import ShrakeRupley

from constants import HYDRO, INTERFACE_RESIDUE_COLUMNS, TIEN, TO1

###################################### SASA AND ATOM GEOMETRY


class InterfaceMixin:
    @staticmethod
    def _standard_residues(chain):
        """Return the standard (non-hetero) residues of a chain."""
        return [residue for residue in chain if residue.id[0] == " "]

    @staticmethod
    def _heavy_atoms(residues):
        """Return heavy-atom coordinates and the residue index owning each atom."""
        coords, owners = [], []
        for index, residue in enumerate(residues):
            for atom in residue:
                if (atom.element or atom.get_name()[:1]).strip().upper() != "H":
                    coords.append(atom.coord)
                    owners.append(index)
        return np.array(coords, float).reshape(-1, 3), np.array(owners, int)

    @staticmethod
    def _total_sasa(chains):
        """Return the total solvent-accessible surface area of the given chains in A^2."""
        model = Model(0)
        for chain in chains:
            model.add(chain.copy())
        ShrakeRupley().compute(model, level="A")
        return float(sum(atom.sasa for atom in model.get_atoms()))

    @staticmethod
    def _residue_sasa(chains):
        """Return a dict mapping (chain_id, resid) to per-residue SASA for the given chains."""
        model = Model(0)
        for chain in chains:
            model.add(chain.copy())
        ShrakeRupley().compute(model, level="R")
        return {(chain.id, residue.id[1]): residue.sasa
                for chain in model for residue in chain if residue.id[0] == " "}

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
        """Return the hydropathy value of a residue, or 0.0 if unknown."""
        return HYDRO.get(TO1.get(residue.get_resname().upper(), "X"), 0.0)

    @staticmethod
    def _split_partners(text):
        """Split a semicolon-separated partner string into a list of non-empty labels."""
        return [p for p in str(text).split(";") if p]

###################################### INTERFACE COMPUTATION

    def _compute_interface(self, chain_a, chain_b):
        """Compute complex-level and residue-level interface metrics for two chains.

        Args:
            chain_a: First chain identifier.
            chain_b: Second chain identifier.
        """
        cutoff = self.config["interface_contact_cutoff_a"]
        score_cutoff = self.config["interface_score_cutoff_a"]
        block = self.config["interface_block_size"]
        A, B = self.chains[chain_a]["chain"], self.chains[chain_b]["chain"]
        ra, rb = self._standard_residues(A), self._standard_residues(B)
        xa, oa = self._heavy_atoms(ra)
        xb, ob = self._heavy_atoms(rb)
        if not len(xa) or not len(xb):
            raise ValueError(f"No heavy atoms in chain {chain_a} or {chain_b}.")
        sa, sb, sab = self._total_sasa([A]), self._total_sasa([B]), self._total_sasa([A, B])
        dasa, initial = sa + sb - sab, sa + sb
        n, score, dmin_c, ha, hb = self._atom_contacts(xa, xb, cutoff, score_cutoff, block)
        dmin = self._residue_min_distances(xa, oa, xb, ob, len(ra), len(rb), block)
        contact = dmin <= cutoff
        mask = dmin <= score_cutoff
        ha_h = np.array([self._hydrophobicity(r) for r in ra])
        hb_h = np.array([self._hydrophobicity(r) for r in rb])
        with np.errstate(invalid="ignore", divide="ignore"):
            energy = np.where(mask, -0.5 * (ha_h[:, None] + hb_h[None]) * (1 - dmin / score_cutoff), 0.0)
        complex_row = {
            "structure": self.structure_id,
            "len_res_a": len(ra), "len_res_b": len(rb),
            "len_atoms_a": len(xa), "len_atoms_b": len(xb),
            "dASA_A2": round(dasa, 2),
            "initial_sasa_A2": round(initial, 2),
            "frac_dASA": round(dasa / initial, 4) if initial > 0 else 0.0,
            "n_contacts": n,
            "iface_res_a": int(np.unique(oa[ha]).size),
            "iface_res_b": int(np.unique(ob[hb]).size),
            "min_dist_A": round(dmin_c, 3),
            "contact_score": round(score, 4),
        }
        iso_a = self._residue_sasa([A])
        iso_b = self._residue_sasa([B])
        cpx = self._residue_sasa([A, B])
        sides = (
            (chain_a, chain_b, ra, rb, iso_a, energy.sum(axis=1), contact.sum(axis=1), dmin.min(axis=1), contact),
            (chain_b, chain_a, rb, ra, iso_b, energy.sum(axis=0), contact.sum(axis=0), dmin.min(axis=0), contact.T),
        )
        rows = []
        for cid, oid, residues, other, iso, en_sum, npart, dm, pairs in sides:
            for k, residue in enumerate(residues):
                if not npart[k]:
                    continue
                key = (cid, residue.id[1])
                dsasa = float(iso[key] - cpx[key])
                maxa = TIEN.get(residue.get_resname().upper())
                en = float(en_sum[k])
                labels = [f"{oid}{other[j].id[1]} {other[j].get_resname()}" for j in np.flatnonzero(pairs[k])]
                rows.append({
                    "structure": self.structure_id,
                    "chain": cid,
                    "resid": residue.id[1],
                    "resname": residue.get_resname(),
                    "dSASA_A2": round(dsasa, 2),
                    "isolated_sasa_A2": round(float(iso[key]), 2),
                    "rel_dSASA": round(dsasa / maxa, 4) if maxa else float("nan"),
                    "n_partners": int(npart[k]),
                    "min_dist_A": round(float(dm[k]), 3),
                    "contact_energy": round(en, 4),
                    "energy_per_partner": round(en / npart[k], 4),
                    "partners": ";".join(labels),
                })
        return {"complex": complex_row, "residues": pd.DataFrame(rows, columns=INTERFACE_RESIDUE_COLUMNS)}

###################################### INTERFACE ACCESSORS

    def interface(self, chain_a="A", chain_b="B"):
        """Return {'complex': dict, 'residues': DataFrame} for two chains, using the cache when possible.

        Args:
            chain_a: First chain identifier.
            chain_b: Second chain identifier.
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

    def interface_complex(self, chain_a="A", chain_b="B"):
        """Return complex-level interface metrics as a dict."""
        return self.interface(chain_a, chain_b)["complex"]

    def interface_residues(self, chain_a="A", chain_b="B"):
        """Return residue-level interface metrics as a DataFrame."""
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
        """Return a dict mapping interface residue labels to their partner labels."""
        frame = self.interface_residues(chain_a, chain_b)
        return {
            f"{row.chain}{row.resid} {row.resname}": self._split_partners(row.partners)
            for row in frame.itertuples(index=False)
        }
