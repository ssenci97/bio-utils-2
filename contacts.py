#!/usr/bin/env python
"""C-alpha and sidechain contact detection, contact strings, and contact summaries with standard and author chain IDs."""

###################################### IMPORTS AND CONSTANTS

import warnings

import mdtraj as md
import numpy as np
import pandas as pd

from constants import CONTACT_COLUMNS, CONTACT_SEQUENCE_COLUMNS, DEFAULT_CONTACT_TYPE, DEFAULT_CUTOFF
from helpers import _check_contact_args, _keep_pair

###################################### C-ALPHA TRAJECTORY


class ContactMixin:
    def _ca_records(self, chain_ids=None):
        """Return (chain_key, sequence_index, residue) tuples for protein residues that have a C-alpha atom.

        Args:
            chain_ids: Chains to include, or None for all.
        """
        selected = None if chain_ids is None else {str(c) for c in chain_ids}
        return [
            (chain_id, index, residue)
            for chain_id, data in self.chains.items()
            if selected is None or chain_id in selected
            for index, residue in enumerate(data["residues"])
            if "CA" in residue
        ]

    def _make_ca_trajectory(self, chain_ids=None):
        """Build an MDTraj C-alpha trajectory in nm and its per-residue records.

        Args:
            chain_ids: Chains to include, or None for all.
        """
        topology = md.Topology()
        md_chains = {}
        coordinates, records = [], []
        for chain_id, index, residue in self._ca_records(chain_ids):
            data = self.chains[chain_id]
            if chain_id not in md_chains:
                md_chains[chain_id] = topology.add_chain(chain_id=chain_id)
            name = residue.resname.strip().upper()
            resnum = int(residue.id[1])
            md_residue = topology.add_residue(name, md_chains[chain_id], resSeq=resnum)
            topology.add_atom("CA", md.element.carbon, md_residue)
            coordinates.append(np.asarray(residue["CA"].coord, dtype=np.float32) / 10.0)
            records.append({"chain_id": data["chain_id"], "auth_chain_id": data["auth_chain_id"],
                            "sequence_index": index, "resnum": resnum, "resname": name})
        if not coordinates:
            raise ValueError("No protein C-alpha atoms found in the selected chains.")
        trajectory = md.Trajectory(xyz=np.asarray(coordinates, dtype=np.float32)[None, :, :], topology=topology)
        return trajectory, records

    def to_ca_trajectory(self, chain_ids=None):
        """Return the C-alpha MDTraj trajectory in nm for the selected chains.

        Args:
            chain_ids: Chains to include, or None for all.
        """
        if chain_ids is not None:
            chain_ids = [self._require_chain(chain_id) for chain_id in chain_ids]
        return self._make_ca_trajectory(chain_ids)[0]

    def _make_sidechain_trajectory(self, chain_ids=None):
        """Build an MDTraj trajectory with all atoms of residues that have a C-alpha, plus per-residue records.

        Args:
            chain_ids: Chains to include, or None for all.
        """
        topology = md.Topology()
        md_chains = {}
        coordinates, records = [], []
        for chain_id, index, residue in self._ca_records(chain_ids):
            data = self.chains[chain_id]
            if chain_id not in md_chains:
                md_chains[chain_id] = topology.add_chain(chain_id=chain_id)
            name = residue.resname.strip().upper()
            resnum = int(residue.id[1])
            md_residue = topology.add_residue(name, md_chains[chain_id], resSeq=resnum)
            for atom in residue.get_atoms():
                letters = "".join(char for char in atom.get_id() if char.isalpha())
                symbol = (atom.element or letters[:1] or "C").strip().capitalize()
                topology.add_atom(atom.get_id().strip(), md.element.get_by_symbol(symbol), md_residue)
                coordinates.append(np.asarray(atom.coord, dtype=np.float32) / 10.0)
            records.append({"chain_id": data["chain_id"], "auth_chain_id": data["auth_chain_id"],
                            "sequence_index": index, "resnum": resnum, "resname": name})
        if not coordinates:
            raise ValueError("No protein atoms found in the selected chains for sidechain contacts.")
        trajectory = md.Trajectory(xyz=np.asarray(coordinates, dtype=np.float32)[None, :, :], topology=topology)
        return trajectory, records

###################################### CONTACTS

    def compute_contacts(self, cutoff=DEFAULT_CUTOFF, contact_type=DEFAULT_CONTACT_TYPE, chain_ids=None, scheme="ca"):
        """Return a DataFrame of residue pairs closer than the cutoff in nm, with standard and author chain IDs.

        Args:
            cutoff: Distance cutoff in nm.
            contact_type: 'inter', 'intra', or 'all'.
            chain_ids: Chains to include, or None for all.
            scheme: mdtraj contact scheme, 'ca' (C-alpha only) or 'sidechain'.
        """
        _check_contact_args(cutoff, contact_type)
        if scheme not in ("ca", "sidechain"):
            raise ValueError(f"Unsupported scheme '{scheme}'. Use 'ca' or 'sidechain'.")
        if chain_ids is not None:
            chain_ids = [self._require_chain(chain_id) for chain_id in chain_ids]
        if scheme == "ca":
            trajectory, records = self._make_ca_trajectory(chain_ids)
        else:
            trajectory, records = self._make_sidechain_trajectory(chain_ids)
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
                scheme=scheme,
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
            })
        return pd.DataFrame(rows, columns=CONTACT_COLUMNS)

###################################### CONTACT STRINGS

    def get_contact_sequence(self, chain_id, cutoff=DEFAULT_CUTOFF, contact_type=DEFAULT_CONTACT_TYPE, contacts=None, scheme="ca"):
        """Return a per-observed-residue string: X for no C-alpha, 0 for no contact, 1 for contact.

        Args:
            chain_id: Chain identifier.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts DataFrame, or None.
            scheme: Passed to compute_contacts when contacts is None.
        """
        chain_id = self._require_chain(chain_id)
        label = self.chains[chain_id]["chain_id"]
        if contacts is None:
            contacts = self.compute_contacts(cutoff=cutoff, contact_type=contact_type, scheme=scheme)
        bits = ["X"] * self.chains[chain_id]["sequence_length"]
        for _, index, _ in self._ca_records([chain_id]):
            bits[index] = "0"
        touched = set()
        for row in contacts.itertuples(index=False):
            if row.chain_1 == label:
                touched.add(int(row.sequence_index_1))
            if row.chain_2 == label:
                touched.add(int(row.sequence_index_2))
        for index in touched:
            if bits[index] == "0":
                bits[index] = "1"
        return "".join(bits)

    def get_reference_contact_sequence(self, chain_id, cutoff=DEFAULT_CUTOFF, contact_type=DEFAULT_CONTACT_TYPE, contacts=None, scheme="ca"):
        """Return the contact string laid out along the full reference sequence.

        Args:
            chain_id: Chain identifier.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts DataFrame, or None.
            scheme: Passed to compute_contacts when contacts is None.
        """
        chain_id = self._require_chain(chain_id)
        data = self.chains[chain_id]
        reference = data["reference_sequence"]
        if not reference:
            raise ValueError(f"Chain {data['chain_id']} has no reference sequence. Set one with set_reference_sequence() first.")
        observed_bits = self.get_contact_sequence(chain_id, cutoff=cutoff, contact_type=contact_type, contacts=contacts, scheme=scheme)
        bits = ["X"] * len(reference)
        for seq_index, ref_index in data["sequence_index_to_ref_index"].items():
            if seq_index < len(observed_bits) and 0 <= ref_index < len(bits):
                bits[ref_index] = observed_bits[seq_index]
        return "".join(bits)

###################################### CONTACT SUMMARY AND EXPORT

    def get_contact_summary(self, cutoff=DEFAULT_CUTOFF, contact_type=DEFAULT_CONTACT_TYPE, contacts=None, scheme="ca"):
        """Return a DataFrame with one row per chain of contact counts and the fraction of residues in contact.

        Args:
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts DataFrame, or None.
            scheme: Passed to compute_contacts when contacts is None.
        """
        if contacts is None:
            contacts = self.compute_contacts(cutoff=cutoff, contact_type=contact_type, scheme=scheme)
        rows = []
        for data in self.chains.values():
            label = data["chain_id"]
            involved = contacts[(contacts["chain_1"] == label) | (contacts["chain_2"] == label)]
            contacted = set(involved.loc[involved["chain_1"] == label, "sequence_index_1"])
            contacted |= set(involved.loc[involved["chain_2"] == label, "sequence_index_2"])
            length = data["sequence_length"]
            rows.append({
                "chain_id": label,
                "auth_chain_id": data["auth_chain_id"],
                "sequence_length": length,
                "reference_length": data["reference_length"],
                "contact_pairs": len(involved),
                "contacted_residues": len(contacted),
                "contact_fraction": len(contacted) / length if length else 0.0,
            })
        return pd.DataFrame(rows)

    def export_contact_sequences(self, output_path, cutoff=DEFAULT_CUTOFF, contact_type=DEFAULT_CONTACT_TYPE, reference=False, scheme="ca"):
        """Write each chain's sequence and contact string to a CSV file and return its path.

        Args:
            output_path: Destination CSV path.
            cutoff: Distance cutoff in nm.
            contact_type: Contact type.
            reference: Whether to use reference sequences, skipping chains without one.
            scheme: Contact scheme ('ca' or 'sidechain').
        """
        contacts = self.compute_contacts(cutoff=cutoff, contact_type=contact_type, scheme=scheme)
        rows = []
        for chain_key, data in self.chains.items():
            if reference:
                if not data["reference_sequence"]:
                    continue
                sequence = data["reference_sequence"]
                contact_string = self.get_reference_contact_sequence(chain_key, cutoff=cutoff, contact_type=contact_type, contacts=contacts, scheme=scheme)
            else:
                sequence = data["sequence"]
                contact_string = self.get_contact_sequence(chain_key, cutoff=cutoff, contact_type=contact_type, contacts=contacts, scheme=scheme)
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