#!/usr/bin/env python
"""Reference sequence and entity lookups derived from mmCIF metadata."""

###################################### IMPORTS AND CONSTANTS

from helpers import cif_sequence, to_list

###################################### METADATA


class MetadataMixin:
    def _entity_sequences(self):
        """Return a dict mapping entity IDs to cleaned reference sequences from mmCIF metadata."""
        entity_ids = to_list(self.mmcif_dict.get("_entity_poly.entity_id"))
        canonical = to_list(self.mmcif_dict.get("_entity_poly.pdbx_seq_one_letter_code_can"))
        original = to_list(self.mmcif_dict.get("_entity_poly.pdbx_seq_one_letter_code"))
        sequences = {}
        for index, entity_id in enumerate(entity_ids):
            raw = canonical[index] if index < len(canonical) else ""
            if not raw and index < len(original):
                raw = original[index]
            sequence = cif_sequence(raw) if raw else ""
            if sequence:
                sequences[entity_id] = sequence
        return sequences

    def _chain_entity_map(self):
        """Return a dict mapping chain IDs to entity IDs from mmCIF metadata."""
        chain_ids = to_list(self.mmcif_dict.get("_struct_asym.id"))
        entity_ids = to_list(self.mmcif_dict.get("_struct_asym.entity_id"))
        return dict(zip(chain_ids, entity_ids))

    def _reference_for_chain(self, chain_id):
        """Return the reference sequence of a chain's entity, or an empty string."""
        entity_id = self._entity_by_chain.get(str(chain_id))
        return self._reference_by_entity.get(entity_id, "")

    def _resolve_entity(self, chain_id, explicit=None):
        """Return the explicit entity ID, else the mmCIF entity, else the default entity ID.

        Args:
            chain_id: Chain identifier.
            explicit: Entity ID given by the caller, or None.
        """
        if explicit is not None:
            return str(explicit)
        return self._entity_by_chain.get(chain_id) or self.default_entity_id

    def set_chain_metadata(self, chain_id, key, value):
        """Store a user-defined value on a chain and return self.

        Args:
            chain_id: Chain identifier.
            key: Metadata key.
            value: Value to store.
        """
        self.chains[self._require_chain(chain_id)]["metadata"][str(key)] = value
        return self

    def set_entity_id(self, chain_id, entity_id):
        """Set a chain's entity ID and return self.

        Args:
            chain_id: Chain identifier.
            entity_id: New entity identifier.
        """
        self.chains[self._require_chain(chain_id)]["entity_id"] = str(entity_id)
        return self
