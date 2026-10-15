#!/usr/bin/env python
"""Chain lifecycle, sequence mapping, chain renaming, and per-chain summaries."""

###################################### IMPORTS AND CONSTANTS

import copy

import pandas as pd

from helpers import clean_sequence, is_protein_residue, residue_to_letter

###################################### CHAIN MANAGEMENT


class ChainMixin:
    def _require_chain(self, chain_id):
        """Return the chain ID as a string, raising KeyError if it is not loaded."""
        chain_id = str(chain_id)
        if chain_id not in self.chains:
            raise KeyError(f"Unknown chain: {chain_id}")
        return chain_id

    @staticmethod
    def _protein_residues(chain):
        """Return the amino acid residues of a chain."""
        return [residue for residue in chain.get_residues() if is_protein_residue(residue)]

    def _map_to_reference(self, sequence, reference):
        """Return a dict mapping observed 0-based indices to reference 0-based indices via global alignment."""
        if not sequence or not reference:
            return {}
        alignment = self._global_aligner.align(reference, sequence)[0]
        return {
            int(seq_index): int(ref_index)
            for ref_index, seq_index in alignment.indices.T
            if ref_index >= 0 and seq_index >= 0
        }

    def _store_mapping(self, data, index_mapping):
        """Save an index mapping and its derived residue-ID mapping and count into chain data.

        Args:
            data: Chain data dict to update.
            index_mapping: Dict from observed to reference indices.
        """
        residues = data["residues"]
        data["sequence_index_to_ref_index"] = index_mapping
        data["resid_to_ref_index"] = {residues[seq_index].id: ref_index for seq_index, ref_index in index_mapping.items()}
        data["mapped_residues"] = len(index_mapping)

    def _make_chain_data(self, chain, reference_sequence=None, entity_id=None, metadata=None):
        """Build the data dict for a chain, including its sequence, reference, and mapping.

        Args:
            chain: Biopython Chain object.
            reference_sequence: Reference sequence, or None to use mmCIF metadata.
            entity_id: Entity ID, or None to resolve it automatically.
            metadata: Optional dict of user metadata.
        """
        chain_id = str(chain.id)
        residues = self._protein_residues(chain)
        if reference_sequence is None:
            reference_sequence = self._reference_for_chain(chain_id)
        data = {
            "chain": chain,
            "chain_id": chain_id,
            "entity_id": self._resolve_entity(chain_id, entity_id),
            "residues": residues,
            "sequence": "".join(residue_to_letter(residue) for residue in residues),
            "reference_sequence": clean_sequence(reference_sequence),
            "metadata": dict(metadata or {}),
        }
        data["sequence_length"] = len(data["sequence"])
        data["reference_length"] = len(data["reference_sequence"])
        self._store_mapping(data, self._map_to_reference(data["sequence"], data["reference_sequence"]))
        return data

    def _build_chains(self):
        """Rebuild chain data from the model, preserving reference, entity, and metadata of existing chains."""
        previous = self.chains
        rebuilt = {}
        for chain in self.model.get_chains():
            chain_id = str(chain.id)
            old = previous.get(chain_id, {})
            rebuilt[chain_id] = self._make_chain_data(
                chain,
                reference_sequence=old.get("reference_sequence"),
                entity_id=old.get("entity_id"),
                metadata=old.get("metadata"),
            )
        self.chains = rebuilt
        self._interface_cache.clear()

    @property
    def chain_ids(self):
        """List of loaded chain IDs."""
        return list(self.chains)

    def __getitem__(self, chain_id):
        """Return the data dict of a chain."""
        return self.chains[str(chain_id)]

    def __iter__(self):
        """Iterate over chain data dicts."""
        return iter(self.chains.values())

    def __len__(self):
        """Return the number of loaded chains."""
        return len(self.chains)

    def __contains__(self, chain_id):
        """Return True if the chain is loaded."""
        return str(chain_id) in self.chains

    def get_sequence(self, chain_id):
        """Return the observed one-letter sequence of a chain."""
        return self.chains[self._require_chain(chain_id)]["sequence"]

    def get_reference_sequence(self, chain_id):
        """Return the reference sequence of a chain."""
        return self.chains[self._require_chain(chain_id)]["reference_sequence"]

    def set_reference_sequence(self, chain_id, sequence, remap=True):
        """Set a chain's reference sequence, optionally rebuilding its mapping, and return self.

        Args:
            chain_id: Chain identifier.
            sequence: New reference sequence.
            remap: Whether to recompute the observed-to-reference mapping.
        """
        data = self.chains[self._require_chain(chain_id)]
        data["reference_sequence"] = clean_sequence(sequence)
        data["reference_length"] = len(data["reference_sequence"])
        index_mapping = self._map_to_reference(data["sequence"], data["reference_sequence"]) if remap else {}
        self._store_mapping(data, index_mapping)
        return self

    def set_reference_sequences(self, sequences, remap=True):
        """Set reference sequences from a {chain_id: sequence} dict and return self.

        Args:
            sequences: Dict mapping chain IDs to reference sequences.
            remap: Whether to recompute each mapping.
        """
        for chain_id, sequence in sequences.items():
            self.set_reference_sequence(chain_id, sequence, remap=remap)
        return self

    def remap_chain(self, chain_id):
        """Rebuild the mapping of one chain against its current reference and return self."""
        data = self.chains[self._require_chain(chain_id)]
        return self.set_reference_sequence(chain_id, data["reference_sequence"], remap=True)

    def remap_all_chains(self):
        """Rebuild the mappings of all chains and return self."""
        for chain_id in list(self.chains):
            self.remap_chain(chain_id)
        return self

    def _apply_chain_options(self, chain_id, reference_sequence, entity_id, metadata):
        """Apply optional reference, entity, and metadata values to a chain.

        Args:
            chain_id: Chain identifier.
            reference_sequence: Reference sequence or None.
            entity_id: Entity ID or None.
            metadata: Metadata dict or None.
        """
        if reference_sequence is not None:
            self.set_reference_sequence(chain_id, reference_sequence)
        if entity_id is not None:
            self.chains[chain_id]["entity_id"] = str(entity_id)
        if metadata:
            self.chains[chain_id]["metadata"].update(metadata)

    def add_chain(self, chain, reference_sequence=None, entity_id=None, metadata=None, replace=False):
        """Add a copy of a Biopython chain to the structure and return self.

        Args:
            chain: Biopython Chain to add.
            reference_sequence: Optional reference sequence for the chain.
            entity_id: Optional entity ID.
            metadata: Optional metadata dict.
            replace: Whether to overwrite an existing chain with the same ID.
        """
        chain_id = str(chain.id)
        if chain_id in self.chains:
            if not replace:
                raise ValueError(f"Chain already exists: {chain_id}")
            self.model.detach_child(chain_id)
        self.model.add(copy.deepcopy(chain))
        self._build_chains()
        self._apply_chain_options(chain_id, reference_sequence, entity_id, metadata)
        return self

    def remove_chain(self, chain_id):
        """Remove a chain from the structure and return self."""
        chain_id = self._require_chain(chain_id)
        self.model.detach_child(chain_id)
        self.chains.pop(chain_id)
        self._alignment_cache.pop(chain_id, None)
        self._interface_cache.clear()
        return self

    def replace_chain(self, chain_id, chain, reference_sequence=None, entity_id=None, metadata=None):
        """Replace an existing chain with a copy of a new Biopython chain and return self.

        Args:
            chain_id: ID of the chain to replace.
            chain: Biopython Chain providing the new coordinates.
            reference_sequence: Optional reference sequence.
            entity_id: Optional entity ID.
            metadata: Optional metadata dict.
        """
        chain_id = self._require_chain(chain_id)
        new_chain = copy.deepcopy(chain)
        new_chain.id = chain_id
        self.model.detach_child(chain_id)
        self.model.add(new_chain)
        self._build_chains()
        self._apply_chain_options(chain_id, reference_sequence, entity_id, metadata)
        return self

###################################### SEQUENCE MATCHING AND RENAMING

    def _pair_identity(self, seq_a, seq_b):
        """Return (identity_percent, aligned_residues, coverage_percent) from a local alignment of two sequences."""
        if not seq_a or not seq_b:
            return 0.0, 0, 0.0
        alignments = self._local_aligner.align(seq_a, seq_b)
        if not alignments:
            return 0.0, 0, 0.0
        top, bottom = alignments[0]
        columns = [(x, y) for x, y in zip(top, bottom) if x != "-" and y != "-"]
        matches = sum(x == y for x, y in columns)
        aligned = len(columns)
        identity = round(matches / aligned * 100, 2) if aligned else 0.0
        coverage = round(aligned / len(seq_b) * 100, 2)
        return identity, aligned, coverage

    def match_chains(self, sequence, min_identity=90.0):
        """Return a DataFrame of chains aligning to a query sequence, best match first.

        Args:
            sequence: Query sequence.
            min_identity: Minimum identity percent to include a chain.
        """
        query = clean_sequence(sequence)
        rows = []
        for chain_id, data in self.chains.items():
            identity, aligned, coverage = self._pair_identity(data["sequence"], query)
            if aligned and identity >= min_identity:
                rows.append({
                    "chain_id": chain_id,
                    "identity_percent": identity,
                    "aligned_residues": aligned,
                    "query_coverage_percent": coverage,
                })
        columns = ["chain_id", "identity_percent", "aligned_residues", "query_coverage_percent"]
        return pd.DataFrame(rows, columns=columns).sort_values("identity_percent", ascending=False, ignore_index=True)

    def rename_chains(self, mapping):
        """Rename chains from an {old_id: new_id} dict, keeping their data, and return self."""
        mapping = {self._require_chain(old): str(new) for old, new in mapping.items()}
        final_ids = [mapping.get(chain_id, chain_id) for chain_id in self.chains]
        if len(set(final_ids)) != len(final_ids):
            raise ValueError("Renaming would create duplicate chain IDs.")
        moved = {old: self.model[old] for old in mapping}
        for old in mapping:
            self.model.detach_child(old)
        for old, new in mapping.items():
            moved[old].id = new
            self.model.add(moved[old])
        self.chains = {mapping.get(c, c): data for c, data in self.chains.items()}
        self._alignment_cache.clear()
        self._build_chains()
        return self

    def rename_by_sequence(self, targets, min_identity=90.0):
        """Rename chains one-to-one to match {new_id: sequence} targets, highest identity first.

        Args:
            targets: Dict mapping new chain IDs to target sequences.
            min_identity: Minimum identity percent required for a rename.

        Returns:
            DataFrame reporting each target and chain with its match status.
        """
        targets = {str(new): clean_sequence(seq) for new, seq in targets.items()}
        candidates = []
        for new_id, sequence in targets.items():
            for old_id, data in self.chains.items():
                identity, aligned, _ = self._pair_identity(data["sequence"], sequence)
                candidates.append((identity, aligned, old_id, new_id))
        candidates.sort(key=lambda row: -row[0])
        renames, used_old, used_new, report = {}, set(), set(), {}
        for identity, aligned, old_id, new_id in candidates:
            if new_id in used_new or old_id in used_old or identity < min_identity:
                continue
            renames[old_id] = new_id
            used_old.add(old_id)
            used_new.add(new_id)
            report[new_id] = (old_id, identity, "renamed")
        rows = []
        for new_id in targets:
            if new_id in report:
                old_id, identity, status = report[new_id]
            else:
                old_id, identity, status = None, None, "unmatched_target"
            rows.append({"new_chain_id": new_id, "old_chain_id": old_id, "identity_percent": identity, "status": status})
        for chain_id in self.chains:
            if chain_id not in used_old:
                rows.append({"new_chain_id": None, "old_chain_id": chain_id, "identity_percent": None, "status": "unmatched_chain"})
        if renames:
            self.rename_chains(renames)
            for old_id, new_id in renames.items():
                self.set_reference_sequence(new_id, targets[new_id])
        return pd.DataFrame(rows, columns=["new_chain_id", "old_chain_id", "identity_percent", "status"])

###################################### SUMMARY

    def summary(self):
        """Return a DataFrame with one row per chain of lengths and mapping coverage."""
        rows = []
        for chain_id, data in self.chains.items():
            length = data["sequence_length"]
            rows.append({
                "chain_id": chain_id,
                "entity_id": data["entity_id"],
                "sequence_length": length,
                "reference_length": data["reference_length"],
                "mapped_residues": data["mapped_residues"],
                "mapping_fraction": round(data["mapped_residues"] / length, 3) if length else 0.0,
            })
        return pd.DataFrame(rows)
