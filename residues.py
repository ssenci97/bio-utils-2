#!/usr/bin/env python
"""Residue numbering, reference-index mapping, and per-residue tables."""

###################################### IMPORTS AND CONSTANTS

import copy

import pandas as pd

from constants import MAPPING_COLUMNS, RESIDUE_COLUMNS
from helpers import is_protein_residue

###################################### RESIDUE MAPPING


class ResidueMixin:
    def _residue_frame(self, chain_id):
        """Return a DataFrame with one row per observed residue of a chain and its reference index."""
        _, _, seq_to_ref = self._local_alignment(chain_id)
        rows = [
            {
                "structure_id": self.structure_id,
                "chain_id": chain_id,
                "pdb_resnum": residue.id[1],
                "resname": residue.resname.strip().upper(),
                "sequence_index": position,
                "reference_index": seq_to_ref.get(position),
            }
            for position, residue in enumerate(self.chains[chain_id]["residues"], start=1)
        ]
        frame = pd.DataFrame(rows, columns=RESIDUE_COLUMNS)
        frame["reference_index"] = frame["reference_index"].astype("Int64")
        return frame

    def get_residue_mapping(self, chain_id):
        """Return a DataFrame with pdb_resnum, sequence_index, and reference_index for one chain."""
        chain_id = self._require_chain(chain_id)
        return self._residue_frame(chain_id)[MAPPING_COLUMNS].reset_index(drop=True)

    @property
    def residue_map(self):
        """Return a DataFrame of all observed residues across chains with their indices."""
        if not self.chains:
            return pd.DataFrame(columns=RESIDUE_COLUMNS)
        return pd.concat([self._residue_frame(chain_id) for chain_id in self.chains], ignore_index=True)

    def to_dataframe(self):
        """Alias for residue_map."""
        return self.residue_map

###################################### RENUMBERING AND TRIMMING

    @staticmethod
    def _reassign(chain, pairs):
        """Set residue IDs from (residue, new_id) pairs in two phases to avoid transient collisions.

        Args:
            chain: Biopython Chain whose residues are renumbered.
            pairs: List of (residue, new_id) tuples.
        """
        for offset, (residue, _) in enumerate(pairs):
            residue.id = ("#", offset, " ")
        chain.child_dict = {residue.id: residue for residue in chain.child_list}
        for residue, new_id in pairs:
            residue.id = new_id
        chain.child_dict = {residue.id: residue for residue in chain.child_list}

    @staticmethod
    def _renumber(chain):
        """Renumber a chain's residues 1..N in place."""
        ResidueMixin._reassign(chain, [(residue, (" ", number, " "))
                                       for number, residue in enumerate(chain.get_residues(), start=1)])

    @classmethod
    def _keep_mapped(cls, chain, index_mapping):
        """Remove unmapped residues from a chain and renumber the kept ones by reference position.

        Args:
            chain: Biopython Chain to trim in place.
            index_mapping: Dict from observed to reference indices.
        """
        protein = [residue for residue in chain.get_residues() if is_protein_residue(residue)]
        pairs = [(protein[index], (" ", ref_index + 1, " "))
                 for index, ref_index in index_mapping.items() if index < len(protein)]
        kept = {id(residue) for residue, _ in pairs}
        for residue in list(chain.get_residues()):
            if id(residue) not in kept:
                chain.detach_child(residue.id)
        cls._reassign(chain, pairs)

    def apply_reference_numbering(self, chain_ids=None):
        """Trim unmapped residues and renumber mapped ones by reference position in place, and return self.

        Args:
            chain_ids: Chains to process, or None for all chains.
        """
        selected = self.chain_ids if chain_ids is None else [self._require_chain(c) for c in chain_ids]
        trimmed = []
        for chain_id in selected:
            mapping = self.chains[chain_id]["sequence_index_to_ref_index"]
            if mapping:
                self._keep_mapped(self.chains[chain_id]["chain"], mapping)
                trimmed.append(chain_id)
        self._build_chains()
        for chain_id in trimmed:
            data = self.chains[chain_id]
            self._store_mapping(data, {index: residue.id[1] - 1 for index, residue in enumerate(data["residues"])})
        return self

    def reindex_chain(self, chain_id, in_place=False):
        """Return a chain with residues renumbered 1..N, modifying the structure only if in_place is True.

        Args:
            chain_id: Chain identifier.
            in_place: Whether to renumber the stored chain itself.
        """
        chain_id = self._require_chain(chain_id)
        chain = self.chains[chain_id]["chain"] if in_place else copy.deepcopy(self.chains[chain_id]["chain"])
        self._renumber(chain)
        if in_place:
            self._build_chains()
        return chain
