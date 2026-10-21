#!/usr/bin/env python
"""C-alpha and sidechain contact detection, contact strings, and contact summaries with standard and author chain IDs."""

###################################### IMPORTS AND CONSTANTS

import warnings

import mdtraj as md
import numpy as np
import pandas as pd

from constants import (
    BACKBONE_ATOMS, CONTACT_COLUMNS, CONTACT_SEQUENCE_COLUMNS, DEFAULT_CONTACT_MOIETIES,
    DEFAULT_CONTACT_TYPE, CONTACTS_CUTOFF_nm, MOIETIES, MOIETY_ION, MOIETY_LIGAND, MOIETY_PROTEIN,
    REPRESENTATIVE_ATOMS,
)
from helpers import _check_contact_args, _keep_pair, residue_moiety
###################################### C-ALPHA TRAJECTORY


class ContactMixin:
    def _ca_records(self, chain_ids=None, moieties=DEFAULT_CONTACT_MOIETIES):
        """Return (chain_key, sequence_index, residue) tuples for residues of the given moieties that have a contact point.

        Args:
            chain_ids: Chains to include, or None for all.
            moieties: Biomolecule types to include, from MOIETIES.
        """
        selected = None if chain_ids is None else {str(c) for c in chain_ids}
        records = []
        for chain_id, data in self.chains.items():
            if selected is not None and chain_id not in selected:
                continue
            for index, residue in enumerate(data["residues"]):
                moiety = residue_moiety(residue)
                if moiety in moieties and self._representative_coord(residue, moiety) is not None:
                    records.append((chain_id, index, residue))
        return records

    @staticmethod
    def _representative_coord(residue, moiety):
        """Return a residue's contact point in A, or None if it has none.

        Protein uses C-alpha, nucleic acids use C4' (or C1', then P), and ligands use the heavy-atom centroid.

        Args:
            residue: Biopython Residue.
            moiety: Biomolecule type from residue_moiety.
        """
        if moiety in (MOIETY_LIGAND, MOIETY_ION):
            coords = [atom.coord for atom in residue.get_atoms()
                      if (atom.element or atom.get_name()[:1]).strip().upper() != "H"]
            return np.mean(coords, axis=0) if coords else None
        for name in REPRESENTATIVE_ATOMS.get(moiety, ()):
            if name in residue:
                return residue[name].coord
        return None

    def _make_ca_trajectory(self, chain_ids=None, moieties=DEFAULT_CONTACT_MOIETIES):
        """Build an MDTraj trajectory of one contact point per residue in nm, plus per-residue records.

        mdtraj's 'ca' scheme selects atoms named CA, so every contact point carries that topology label
        whatever its moiety. The coordinate is the real contact point (see _representative_coord).

        Args:
            chain_ids: Chains to include, or None for all.
            moieties: Biomolecule types to include.
        """
        topology = md.Topology()
        md_chains = {}
        coordinates, records = [], []
        for chain_id, index, residue in self._ca_records(chain_ids, moieties):
            data = self.chains[chain_id]
            moiety = residue_moiety(residue)
            if chain_id not in md_chains:
                md_chains[chain_id] = topology.add_chain(chain_id=chain_id)
            name = residue.resname.strip().upper()
            resnum = int(residue.id[1])
            md_residue = topology.add_residue(name, md_chains[chain_id], resSeq=resnum)
            topology.add_atom("CA", md.element.carbon, md_residue)
            coordinates.append(np.asarray(self._representative_coord(residue, moiety), dtype=np.float32) / 10.0)
            records.append({"chain_id": data["chain_id"], "auth_chain_id": data["auth_chain_id"],
                            "sequence_index": index, "resnum": resnum, "resname": name, "moiety": moiety})
        if not coordinates:
            raise ValueError("No contact points found in the selected chains for the requested moieties.")
        trajectory = md.Trajectory(xyz=np.asarray(coordinates, dtype=np.float32)[None, :, :], topology=topology)
        return trajectory, records

    def to_ca_trajectory(self, chain_ids=None, moieties=DEFAULT_CONTACT_MOIETIES):
        """Return the contact-point MDTraj trajectory in nm for the selected chains.

        Args:
            chain_ids: Chains to include, or None for all.
            moieties: Biomolecule types to include.
        """
        if chain_ids is not None:
            chain_ids = [self._require_chain(chain_id) for chain_id in chain_ids]
        return self._make_ca_trajectory(chain_ids, moieties)[0]

    def _make_sidechain_trajectory(self, chain_ids=None, moieties=DEFAULT_CONTACT_MOIETIES):
        """Build a trajectory of sidechain atoms for protein and one contact point for other moieties.

        Protein residues contribute their heavy sidechain atoms (CA for Gly). Nucleic, ligand, and ion
        residues contribute the same contact point used by the 'ca' scheme. Each residue is one mdtraj
        residue, so the records stay aligned with the topology.

        Args:
            chain_ids: Chains to include, or None for all.
            moieties: Biomolecule types to include.
        """
        topology = md.Topology()
        md_chains = {}
        coordinates, records = [], []
        for chain_id, index, residue in self._ca_records(chain_ids, moieties):
            data = self.chains[chain_id]
            moiety = residue_moiety(residue)
            if chain_id not in md_chains:
                md_chains[chain_id] = topology.add_chain(chain_id=chain_id)
            name = residue.resname.strip().upper()
            resnum = int(residue.id[1])
            md_residue = topology.add_residue(name, md_chains[chain_id], resSeq=resnum)
            if moiety == MOIETY_PROTEIN:
                atoms = [atom for atom in residue.get_atoms()
                         if atom.get_id().strip() not in BACKBONE_ATOMS
                         and (atom.element or atom.get_id()[:1]).strip().upper() != "H"]
                if not atoms:
                    atoms = [residue["CA"]]
                for atom in atoms:
                    letters = "".join(char for char in atom.get_id() if char.isalpha())
                    symbol = (atom.element or letters[:1] or "C").strip().capitalize()
                    topology.add_atom(atom.get_id().strip(), md.element.get_by_symbol(symbol), md_residue)
                    coordinates.append(np.asarray(atom.coord, dtype=np.float32) / 10.0)
            else:
                topology.add_atom("CP", md.element.carbon, md_residue)
                coordinates.append(np.asarray(self._representative_coord(residue, moiety), dtype=np.float32) / 10.0)
            records.append({"chain_id": data["chain_id"], "auth_chain_id": data["auth_chain_id"],
                            "sequence_index": index, "resnum": resnum, "resname": name, "moiety": moiety})
        if not coordinates:
            raise ValueError("No sidechain or contact points found in the selected chains for the requested moieties.")
        trajectory = md.Trajectory(xyz=np.asarray(coordinates, dtype=np.float32)[None, :, :], topology=topology)
        return trajectory, records

###################################### CONTACTS

    def compute_contacts(self, cutoff=CONTACTS_CUTOFF_nm, contact_type=DEFAULT_CONTACT_TYPE, chain_ids=None,
                         scheme="ca", moieties=DEFAULT_CONTACT_MOIETIES):
        """Return a DataFrame of residue pairs closer than the cutoff in nm, with standard and author chain IDs.

        Args:
            cutoff: Distance cutoff in nm.
            contact_type: 'inter', 'intra', or 'all'.
            chain_ids: Chains to include, or None for all.
            scheme: 'ca' (one contact point per residue) or 'sidechain' (protein sidechain atoms plus one
                contact point per other moiety).
            moieties: Biomolecule types to include, from MOIETIES.
        """
        _check_contact_args(cutoff, contact_type)
        if scheme not in ("ca", "sidechain"):
            raise ValueError(f"Unsupported scheme '{scheme}'. Use 'ca' or 'sidechain'.")
        if not set(moieties) <= set(MOIETIES):
            raise ValueError(f"moieties must be a subset of {MOIETIES}.")
        if chain_ids is not None:
            chain_ids = [self._require_chain(chain_id) for chain_id in chain_ids]
        moieties = tuple(moieties)
        if scheme == "ca":
            trajectory, records = self._make_ca_trajectory(chain_ids, moieties)
            md_scheme = "ca"
        else:
            trajectory, records = self._make_sidechain_trajectory(chain_ids, moieties)
            md_scheme = "closest-heavy"
        pairs = [
            (i, j)
            for i in range(len(records))
            for j in range(i + 1, len(records))
            if _keep_pair(records[i]["chain_id"], records[j]["chain_id"], contact_type)
        ]
        if not pairs:
            return pd.DataFrame(columns=CONTACT_COLUMNS)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, message=".*No alpha carbon.*")
            distances, residue_pairs = md.compute_contacts(
                trajectory,
                contacts=np.asarray(pairs, dtype=int),
                scheme=md_scheme,
                ignore_nonprotein=False,
                periodic=False,
            )
        rows = []
        for (first, second), distance in zip(residue_pairs, distances[0]):
            if distance > cutoff:
                continue
            a, b = records[int(first)], records[int(second)]
            rows.append({
                "chain_1": a["chain_id"],
                "auth_chain_1": a["auth_chain_id"],
                "resnum_1": a["resnum"],
                "resname_1": a["resname"],
                "sequence_index_1": a["sequence_index"],
                "chain_2": b["chain_id"],
                "auth_chain_2": b["auth_chain_id"],
                "resnum_2": b["resnum"],
                "resname_2": b["resname"],
                "sequence_index_2": b["sequence_index"],
                "min_distance_nm": round(float(distance), 4),
                "contact_type": "intra" if a["chain_id"] == b["chain_id"] else "inter",
                "moiety_1": a["moiety"],
                "moiety_2": b["moiety"],
            })
        return pd.DataFrame(rows, columns=CONTACT_COLUMNS)

    def _string_indices(self, chain_id, moieties=DEFAULT_CONTACT_MOIETIES):
        """Return sequence indices of a chain's residues whose moiety is included in contact strings.

        Args:
            chain_id: Chain key.
            moieties: Biomolecule types included in the strings.
        """
        residues = self.chains[chain_id]["residues"]
        return [index for index, residue in enumerate(residues) if residue_moiety(residue) in moieties]

    def _observed_bits(self, chain_id, cutoff, contact_type, contacts, scheme, moieties):
        """Return one bit per observed residue of a chain: X, 0, or 1, before filtering by moiety.

        Args:
            chain_id: Chain key.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts DataFrame, or None.
            scheme: Contact scheme, used if contacts is None.
            moieties: Biomolecule types, used if contacts is None.
        """
        label = self.chains[chain_id]["chain_id"]
        if contacts is None:
            contacts = self.compute_contacts(cutoff=cutoff, contact_type=contact_type,
                                             scheme=scheme, moieties=moieties)
        bits = ["X"] * self.chains[chain_id]["sequence_length"]
        for _, index, _ in self._ca_records([chain_id], moieties):
            bits[index] = "0"
        touched = set()
        for row in contacts.itertuples(index=False):
            if row.chain_1 == label:
                touched.add(int(row.sequence_index_1))
            if row.chain_2 == label:
                touched.add(int(row.sequence_index_2))
        for index in touched:
            bits[index] = "1"
        return bits

    def get_contact_sequence(self, chain_id, cutoff=CONTACTS_CUTOFF_nm, contact_type=DEFAULT_CONTACT_TYPE,
                             contacts=None, scheme="ca", moieties=DEFAULT_CONTACT_MOIETIES):
        """Return a per-residue string of X (no contact point), 0 (no contact), or 1 (contact).

        Only residues whose moiety is in moieties appear in the string. Ions and water are left out by default.

        Args:
            chain_id: Chain identifier.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts DataFrame, or None.
            scheme: Contact scheme, used if contacts is None.
            moieties: Biomolecule types included in the string.
        """
        chain_id = self._require_chain(chain_id)
        bits = self._observed_bits(chain_id, cutoff, contact_type, contacts, scheme, moieties)
        return "".join(bits[index] for index in self._string_indices(chain_id, moieties))

    def get_reference_contact_sequence(self, chain_id, cutoff=CONTACTS_CUTOFF_nm, contact_type=DEFAULT_CONTACT_TYPE,
                                       contacts=None, scheme="ca", moieties=DEFAULT_CONTACT_MOIETIES):
        """Return the contact string laid out along the full reference sequence.

        Reference positions not mapped to an included observed residue are 'X'. Mapped positions are
        '1' if in contact and '0' otherwise, including residues without a contact point.

        Args:
            chain_id: Chain identifier.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts DataFrame, or None.
            scheme: Contact scheme, used if contacts is None.
            moieties: Biomolecule types included in the string.
        """
        chain_id = self._require_chain(chain_id)
        data = self.chains[chain_id]
        reference = data["reference_sequence"]
        if not reference:
            raise ValueError(f"Chain {data['chain_id']} has no reference sequence. Set one with set_reference_sequence() first.")
        observed = self._observed_bits(chain_id, cutoff, contact_type, contacts, scheme, moieties)
        kept = set(self._string_indices(chain_id, moieties))
        bits = ["X"] * len(reference)
        for seq_index, ref_index in data["sequence_index_to_ref_index"].items():
            if seq_index in kept and seq_index < len(observed) and 0 <= ref_index < len(bits):
                bits[ref_index] = "1" if observed[seq_index] == "1" else "0"
        return "".join(bits)

    def export_contact_sequences(self, output_path, cutoff=CONTACTS_CUTOFF_nm, contact_type=DEFAULT_CONTACT_TYPE,
                                 reference=False, scheme="ca", moieties=DEFAULT_CONTACT_MOIETIES):
        """Write each chain's sequence and contact string to a CSV file and return its path.

        Sequences and strings both include only residues whose moiety is in moieties.

        Args:
            output_path: Destination CSV path.
            cutoff: Distance cutoff in nm.
            contact_type: Contact type.
            reference: Whether to use reference sequences, skipping chains without one.
            scheme: Contact scheme ('ca' or 'sidechain').
            moieties: Biomolecule types included in the output.
        """
        contacts = self.compute_contacts(cutoff=cutoff, contact_type=contact_type, scheme=scheme, moieties=moieties)
        rows = []
        for chain_key, data in self.chains.items():
            if reference:
                if not data["reference_sequence"]:
                    continue
                sequence = data["reference_sequence"]
                contact_string = self.get_reference_contact_sequence(chain_key, cutoff=cutoff, contact_type=contact_type,
                                                                     contacts=contacts, scheme=scheme, moieties=moieties)
            else:
                indices = self._string_indices(chain_key, moieties)
                sequence = "".join(data["sequence"][index] for index in indices if index < len(data["sequence"]))
                contact_string = self.get_contact_sequence(chain_key, cutoff=cutoff, contact_type=contact_type,
                                                           contacts=contacts, scheme=scheme, moieties=moieties)
            rows.append({
                "structure_id": self.structure_id,
                "chain_id": data["chain_id"],
                "auth_chain_id": data["auth_chain_id"],
                "sequence": sequence,
                "contacts": contact_string,
                "contact_type": contact_type,
                "mdtraj_scheme": scheme,
                "cutoff_nm": cutoff,
            })
        output_path = self._prepare_path(output_path)
        pd.DataFrame(rows, columns=CONTACT_SEQUENCE_COLUMNS).to_csv(output_path, index=False)
        return output_path