#!/usr/bin/env python
"""Chain lifecycle, sequence mapping, chain renaming, external references, and per-chain summaries keyed by author IDs with standard IDs alongside."""

###################################### IMPORTS AND CONSTANTS

import copy

import pandas as pd

from helpers import clean_sequence, is_protein_residue, residue_to_letter

###################################### CHAIN MANAGEMENT


class ChainMixin:
    def _require_chain(self, chain_id):
        """Return the chain key matching a key, standard ID, or author ID, raising KeyError if not loaded.

        Args:
            chain_id: Chain key, standard chain ID, or author chain ID.
        """
        chain_id = str(chain_id)
        if chain_id in self.chains:
            return chain_id
        for field in ("chain_id", "auth_chain_id"):
            for key, data in self.chains.items():
                if data[field] == chain_id:
                    return key
        raise KeyError(f"Unknown chain: {chain_id}")

    @staticmethod
    def _protein_residues(chain):
        """Return the amino acid residues of a chain."""
        return [residue for residue in chain.get_residues() if is_protein_residue(residue)]

    def _map_to_reference(self, sequence, reference):
        """Return {observed 0-based index: reference 0-based index} from a global alignment.

        Args:
            sequence: Observed sequence.
            reference: Reference sequence.
        """
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

    def _make_chain_data(self, chain, reference_sequence=None, entity_id=None, metadata=None, external_ref=None):
        """Build the data dict for a chain, with its key, standard ID, author ID, and optional external reference.

        Args:
            chain: Biopython Chain object, keyed by its chain key.
            reference_sequence: Reference sequence, or None to use mmCIF metadata.
            entity_id: Entity ID, or None to resolve it automatically.
            metadata: Optional dict of user metadata.
            external_ref: User-set external reference, or None to use the mmCIF reference.
        """
        key = str(chain.id)
        residues = self._protein_residues(chain)
        if reference_sequence is None:
            reference_sequence = self._reference_for_chain(key)
        data = {
            "chain": chain,
            "chain_id": self._label_of(key),
            "auth_chain_id": self._auth_of(key),
            "entity_id": self._resolve_entity(key, entity_id),
            "external_ref": external_ref,
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
        """Rebuild chain data from the model, preserving reference, entity, metadata, and external reference of existing chains."""
        previous = self.chains
        rebuilt = {}
        for chain in self.model.get_chains():
            key = str(chain.id)
            old = previous.get(key, {})
            rebuilt[key] = self._make_chain_data(
                chain,
                reference_sequence=old.get("reference_sequence"),
                entity_id=old.get("entity_id"),
                metadata=old.get("metadata"),
                external_ref=old.get("external_ref"),
            )
        self.chains = rebuilt
        self._interface_cache.clear()

    @property
    def chain_ids(self):
        """List of loaded chain keys."""
        return list(self.chains)

    def __getitem__(self, chain_id):
        """Return the data dict of a chain."""
        return self.chains[self._require_chain(chain_id)]

    def __iter__(self):
        """Iterate over chain data dicts."""
        return iter(self.chains.values())

    def __len__(self):
        """Return the number of loaded chains."""
        return len(self.chains)

    def __contains__(self, chain_id):
        """Return True if the chain is loaded."""
        try:
            self._require_chain(chain_id)
            return True
        except KeyError:
            return False

    def get_sequence(self, chain_id):
        """Return the observed one-letter sequence of a chain."""
        return self.chains[self._require_chain(chain_id)]["sequence"]

    def get_reference_sequence(self, chain_id):
        """Return the reference sequence of a chain."""
        return self.chains[self._require_chain(chain_id)]["reference_sequence"]

    def get_external_ref(self, chain_id):
        """Return the external reference of a chain: the user-set value if present, else the mmCIF accession.

        Args:
            chain_id: Chain identifier, standard or author.
        """
        data = self.chains[self._require_chain(chain_id)]
        if data["external_ref"]:
            return data["external_ref"]
        exts = self._external_refs_for_entity(data["entity_id"])
        return ", ".join(exts) if exts else None

    def _external_refs_for_entity(self, entity_id):
        """Return the mmCIF external accessions linked to an entity, or an empty list.

        Args:
            entity_id: Entity identifier.
        """
        if not hasattr(self, "_external_ref_map"):
            self._chain_entity_map()
        return [ext for ext, ent in self._external_ref_map.items() if ent == entity_id]

    def set_external_ref(self, chain_id, external_ref):
        """Set a chain's external reference, overriding the mmCIF value in summaries, and return self.

        Args:
            chain_id: Chain identifier, standard or author.
            external_ref: External accession string, or None to fall back to the mmCIF reference.
        """
        chain_id = self._require_chain(chain_id)
        self.chains[chain_id]["external_ref"] = None if external_ref is None else str(external_ref)
        return self

    def set_reference_sequence(self, chain_id, sequence, remap=True):
        """Set a chain's reference sequence, optionally rebuilding its mapping, and return self.

        Args:
            chain_id: Chain identifier, standard or author.
            sequence: New reference sequence.
            remap: Whether to rebuild the mapping.
        """
        chain_id = self._require_chain(chain_id)
        data = self.chains[chain_id]
        data["reference_sequence"] = clean_sequence(sequence)
        data["reference_length"] = len(data["reference_sequence"])
        index_mapping = self._map_to_reference(data["sequence"], data["reference_sequence"]) if remap else {}
        self._store_mapping(data, index_mapping)
        return self

    def set_reference_sequences(self, sequences, remap=True):
        """Set reference sequences from a {chain_id: sequence} dict and return self."""
        for chain_id, sequence in sequences.items():
            self.set_reference_sequence(chain_id, sequence, remap=remap)
        return self

    def remap_chain(self, chain_id):
        """Rebuild the mapping of one chain against its current reference and return self."""
        chain_id = self._require_chain(chain_id)
        data = self.chains[chain_id]
        print(f"remapped {chain_id}")
        return self.set_reference_sequence(chain_id, data["reference_sequence"], remap=True)

    def remap_all_chains(self):
        """Rebuild the mappings of all chains and return self."""
        for chain_id in list(self.chains):
            self.remap_chain(chain_id)
        return self

    def _apply_chain_options(self, chain_id, reference_sequence, entity_id, metadata):
        """Apply optional reference, entity, and metadata values to a chain.

        Args:
            chain_id: Chain identifier, standard or author.
            reference_sequence: Reference sequence, or None.
            entity_id: Entity ID, or None.
            metadata: Metadata dict, or None.
        """
        chain_id = self._require_chain(chain_id)
        if reference_sequence is not None:
            self.set_reference_sequence(chain_id, reference_sequence)
        if entity_id is not None:
            self.chains[chain_id]["entity_id"] = str(entity_id)
        if metadata:
            self.chains[chain_id]["metadata"].update(metadata)

    def add_chain(self, chain, reference_sequence=None, entity_id=None, metadata=None, replace=False):
        """Add a Biopython chain, optionally replacing an existing one, and return self.

        Args:
            chain: Biopython Chain object.
            reference_sequence: Optional reference sequence.
            entity_id: Optional entity ID.
            metadata: Optional metadata dict.
            replace: Whether to replace an existing chain with the same ID.
        """
        key = str(chain.id)
        try:
            existing = self._require_chain(key)
        except KeyError:
            existing = None
        if existing is not None:
            if not replace:
                raise ValueError(f"Chain already exists: {key}")
            self.remove_chain(existing)
        self.model.add(copy.deepcopy(chain))
        self._build_chains()
        self._apply_chain_options(key, reference_sequence, entity_id, metadata)
        return self

    def remove_chain(self, chain_id):
        """Remove a chain and return self."""
        chain_id = self._require_chain(chain_id)
        self.model.detach_child(chain_id)
        self.chains.pop(chain_id, None)
        self._alignment_cache.pop(chain_id, None)
        self._interface_cache.clear()
        return self

    def replace_chain(self, chain_id, chain, reference_sequence=None, entity_id=None, metadata=None):
        """Replace a chain with a new Biopython chain and return self.

        Args:
            chain_id: Chain identifier, standard or author.
            chain: Replacement Biopython Chain.
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
        """Return (identity_percent, aligned_residues, coverage_percent_of_b) from a local alignment.

        Args:
            seq_a: First sequence.
            seq_b: Second sequence.
        """
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
        """Return a DataFrame of chains aligning to a query sequence with standard and author IDs.

        Args:
            sequence: Query sequence.
            min_identity: Minimum identity percent to report.
        """
        query = clean_sequence(sequence)
        rows = []
        for data in self.chains.values():
            identity, aligned, coverage = self._pair_identity(data["sequence"], query)
            if aligned and identity >= min_identity:
                rows.append({
                    "chain_id": data["chain_id"],
                    "auth_chain_id": data["auth_chain_id"],
                    "identity_percent": identity,
                    "aligned_residues": aligned,
                    "query_coverage_percent": coverage,
                })
        columns = ["chain_id", "auth_chain_id", "identity_percent", "aligned_residues", "query_coverage_percent"]
        return pd.DataFrame(rows, columns=columns).sort_values("identity_percent", ascending=False, ignore_index=True)

    def rename_chains(self, mapping):
        """Rename chain keys with a {old: new} dict, keeping standard and author ID maps, entities, and sequences in sync.

        Args:
            mapping: Dict of {old_chain_identifier: new_key}.
        """
        resolved = {self._require_chain(old): str(new) for old, new in mapping.items()}
        final_ids = [resolved.get(key, key) for key in self.chains]
        if len(set(final_ids)) != len(final_ids):
            raise ValueError("Renaming would create duplicate chain IDs.")
        moved = {old: self.model[old] for old in resolved}
        for old in resolved:
            self.model.detach_child(old)
        key_map = self._auth_to_label if self.key_by == "auth" else self._label_to_auth
        other_map = self._label_to_auth if self.key_by == "auth" else self._auth_to_label
        for old, new in resolved.items():
            moved[old].id = new
            self.model.add(moved[old])
            other = key_map.pop(old, None)
            if other is not None:
                key_map[new] = other
                other_map[other] = new
            entity = self._entity_by_chain.pop(old, None)
            if entity is not None:
                self._entity_by_chain[new] = entity
        self.chains = {resolved.get(key, key): data for key, data in self.chains.items()}
        self._alignment_cache.clear()
        self._build_chains()
        return self

    def rename_by_sequence(self, targets, min_identity=90.0):
        """Rename chain keys to match a {new_id: sequence} dict and report standard and author IDs of the old chains.

        Args:
            targets: Dict of {new_key: sequence}.
            min_identity: Minimum identity percent required for a rename.
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
                old = self.chains[old_id]
                rows.append({"new_chain_id": new_id, "old_chain_id": old["chain_id"],
                             "old_auth_chain_id": old["auth_chain_id"], "identity_percent": identity, "status": status})
            else:
                rows.append({"new_chain_id": new_id, "old_chain_id": None, "old_auth_chain_id": None,
                             "identity_percent": None, "status": "unmatched_target"})
        for key, data in self.chains.items():
            if key not in used_old:
                rows.append({"new_chain_id": None, "old_chain_id": data["chain_id"],
                             "old_auth_chain_id": data["auth_chain_id"], "identity_percent": None,
                             "status": "unmatched_chain"})
        if renames:
            self.rename_chains(renames)
            for old_id, new_id in renames.items():
                self.set_reference_sequence(new_id, targets[new_id])
        return pd.DataFrame(rows, columns=["new_chain_id", "old_chain_id", "old_auth_chain_id", "identity_percent", "status"])

###################################### SUMMARY

    def summary(self):
        """Return a summary DataFrame with standard and author chain IDs, entity, description, external reference, and mapping coverage."""
        desc_map = self._chain_description_map()
        rows = []
        for key, data in self.chains.items():
            length = data["sequence_length"]
            rows.append({
                "chain_id": data["chain_id"],
                "auth_chain_id": data["auth_chain_id"],
                "entity_id": data["entity_id"],
                "description": desc_map.get(key),
                "external_ref": self.get_external_ref(key),
                "sequence_length": length,
                "reference_length": data["reference_length"],
                "mapped_residues": data["mapped_residues"],
                "mapping_fraction": round(data["mapped_residues"] / length, 3) if length else 0.0,
            })
        return pd.DataFrame(rows)