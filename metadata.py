#!/usr/bin/env python
"""Reference sequence, entity, and description lookups from mmCIF metadata, resolving label IDs for each chain key."""

###################################### IMPORTS AND CONSTANTS

from collections import Counter

from helpers import cif_sequence, to_list

###################################### METADATA


class MetadataMixin:
    def _entity_sequences(self):
        """Return {entity_id: cleaned reference sequence} from mmCIF metadata.

        Returns:
            Dict of {entity_id: sequence} for entities with a non-empty sequence.
        """
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

    def _clean_val(self, val):
        """Return None for mmCIF missing-value markers ('?' or '.'), otherwise the value unchanged.

        Args:
            val: Raw mmCIF value.
        """
        if val in ("?", ".", None):
            return None
        return val

    def _atom_pair_maps(self):
        """Build label-to-author and author-to-label maps from atom_site, preferring polymer atoms.

        Each author ID maps to the label ID that occurs most often among its polymer atoms, so a
        single author ID resolves to one label ID. Hetero atoms only fill in IDs with no polymer atoms.
        """
        labels = to_list(self.mmcif_dict.get("_atom_site.label_asym_id"))
        auths = to_list(self.mmcif_dict.get("_atom_site.auth_asym_id"))
        groups = to_list(self.mmcif_dict.get("_atom_site.group_PDB")) or ["ATOM"] * len(labels)
        polymer, hetero = Counter(), Counter()
        for group, label, auth in zip(groups, labels, auths):
            label, auth = self._clean_val(label), self._clean_val(auth)
            if label and auth:
                (polymer if group == "ATOM" else hetero)[(label, auth)] += 1
        self._auth_to_label, self._label_to_auth = {}, {}
        for counts in (polymer, hetero):
            for (label, auth), _ in counts.most_common():
                self._auth_to_label.setdefault(auth, label)
                self._label_to_auth.setdefault(label, auth)

    def _chain_entity_map(self):
        """Return {chain_key: entity_id} for every chain key, resolving each key through its label ID.

        Returns:
            Dict of {chain_key: entity_id}. Keys are author IDs unless key_by is 'label'.
        """
        if not hasattr(self, "_auth_to_label"):
            self._atom_pair_maps()
        labels = to_list(self.mmcif_dict.get("_struct_asym.id"))
        entities = to_list(self.mmcif_dict.get("_struct_asym.entity_id"))
        label_entity = {
            self._clean_val(label): self._clean_val(entity)
            for label, entity in zip(labels, entities)
            if self._clean_val(label) and self._clean_val(entity)
        }
        if self.key_by == "auth":
            pairs = [(key, label) for key, label in self._auth_to_label.items()]
        else:
            pairs = [(key, key) for key in self._label_to_auth]
        chain_map = {}
        for key, label in pairs:
            entity = label_entity.get(label)
            if entity:
                chain_map[key] = entity

        ref_ents = to_list(self.mmcif_dict.get("_struct_ref.entity_id"))
        db_accs = to_list(self.mmcif_dict.get("_struct_ref.pdbx_db_accession"))
        external_map = {}
        for entity, accession in zip(ref_ents, db_accs):
            entity_clean = self._clean_val(entity)
            accession_clean = self._clean_val(accession)
            if entity_clean and accession_clean:
                external_map[accession_clean] = entity_clean
        self._external_ref_map = external_map
        return chain_map

    def _label_of(self, chain_key):
        """Return the label chain ID for a chain key, or the key itself if unknown.

        Args:
            chain_key: Chain key as held by the structure.
        """
        if not hasattr(self, "_auth_to_label"):
            self._atom_pair_maps()
        key = str(chain_key)
        return key if self.key_by == "label" else self._auth_to_label.get(key, key)

    def _auth_of(self, chain_key):
        """Return the author chain ID for a chain key, or the key itself if unknown.

        Args:
            chain_key: Chain key as held by the structure.
        """
        if not hasattr(self, "_label_to_auth"):
            self._atom_pair_maps()
        key = str(chain_key)
        return key if self.key_by == "auth" else self._label_to_auth.get(key, key)

    def _chain_description_map(self):
        """Return {chain_key: entity description} for every chain key.

        Returns:
            Dict of {chain_key: description_or_None}.
        """
        ent_ids = to_list(self.mmcif_dict.get("_entity.id"))
        ent_desc = to_list(self.mmcif_dict.get("_entity.pdbx_description"))
        ent_map = {
            self._clean_val(entity): self._clean_val(desc)
            for entity, desc in zip(ent_ids, ent_desc)
            if self._clean_val(entity)
        }
        return {key: ent_map.get(entity) for key, entity in self._chain_entity_map().items()}

    def _reference_for_chain(self, chain_id):
        """Return the reference sequence of a chain's entity, or an empty string.

        Args:
            chain_id: Chain key as held by the structure.
        """
        entity_id = self._entity_by_chain.get(str(chain_id))
        return self._reference_by_entity.get(entity_id, "")

    def _resolve_entity(self, chain_id, explicit=None):
        """Return the explicit entity ID, else the mmCIF entity, else the default entity ID.

        Args:
            chain_id: Chain key as held by the structure.
            explicit: Entity ID given by the caller, or None.
        """
        if explicit is not None:
            return str(explicit)
        return self._entity_by_chain.get(str(chain_id)) or self.default_entity_id

    def set_chain_metadata(self, chain_id, key, value):
        """Store a user-defined value on a chain and return self.

        Args:
            chain_id: Chain identifier, standard or author.
            key: Metadata key.
            value: Value to store.
        """
        self.chains[self._require_chain(chain_id)]["metadata"][str(key)] = value
        return self

    def set_entity_id(self, chain_id, entity_id):
        """Set a chain's entity ID and return self.

        Args:
            chain_id: Chain identifier, standard or author.
            entity_id: New entity identifier.
        """
        self.chains[self._require_chain(chain_id)]["entity_id"] = str(entity_id)
        return self