"""Reference sequence, entity, external reference, and description lookups from mmCIF metadata, resolving label IDs for each chain key."""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

from collections import Counter

from .base import StructureState
from .sequence import clean_cif_sequence, clean_cif_value, to_list

###################################### METADATA


class MetadataMixin(StructureState):
    """Lookups from mmCIF metadata: entities, label and author ID maps, external references, and descriptions."""

    def _field(self, name: str) -> list[str]:
        """Return an mmCIF field as a list of strings.

        Args:
            name: mmCIF item name, such as '_entity.id'.

        Returns:
            The field values; empty if the field is missing.
        """
        return to_list(self.mmcif_dict.get(name))

    def _entity_sequences(self) -> dict[str, str]:
        """Return the cleaned reference sequence of each entity that has one.

        Returns:
            {entity_id: sequence}.
        """
        canonical = self._field("_entity_poly.pdbx_seq_one_letter_code_can")
        original = self._field("_entity_poly.pdbx_seq_one_letter_code")
        sequences: dict[str, str] = {}
        for index, entity_id in enumerate(self._field("_entity_poly.entity_id")):
            raw = canonical[index] if index < len(canonical) else ""
            if not raw and index < len(original):
                raw = original[index]
            sequence = clean_cif_sequence(raw) if raw else ""
            if sequence:
                sequences[entity_id] = sequence
        return sequences

    def _build_chain_id_maps(self) -> None:
        """Build label-to-author and author-to-label maps from atom_site, preferring polymer atoms.

        Each author ID maps to the label ID that occurs most often among its polymer atoms, so a
        single author ID resolves to one label ID. Hetero atoms only fill in IDs with no polymer atoms.
        """
        labels = self._field("_atom_site.label_asym_id")
        auths = self._field("_atom_site.auth_asym_id")
        groups = self._field("_atom_site.group_PDB") or ["ATOM"] * len(labels)
        polymer: Counter[tuple[str, str]] = Counter()
        hetero: Counter[tuple[str, str]] = Counter()
        for group, raw_label, raw_auth in zip(groups, labels, auths):
            label, auth = clean_cif_value(raw_label), clean_cif_value(raw_auth)
            if label and auth:
                (polymer if group == "ATOM" else hetero)[(label, auth)] += 1
        self._auth_to_label, self._label_to_auth = {}, {}
        for counts in (polymer, hetero):
            for (label, auth), _ in counts.most_common():
                self._auth_to_label.setdefault(auth, label)
                self._label_to_auth.setdefault(label, auth)

    def _chain_entity_map(self) -> dict[str, str]:
        """Return the entity of every chain key, resolving each key through its label ID.

        Returns:
            {chain_key: entity_id}.
        """
        label_entity: dict[str, str] = {}
        for raw_label, raw_entity in zip(self._field("_struct_asym.id"), self._field("_struct_asym.entity_id")):
            label, entity = clean_cif_value(raw_label), clean_cif_value(raw_entity)
            if label and entity:
                label_entity[label] = entity
        if self.key_by == "auth":
            pairs = list(self._auth_to_label.items())
        else:
            pairs = [(key, key) for key in self._label_to_auth]
        return {key: label_entity[label] for key, label in pairs if label in label_entity}

    def _entity_accessions(self) -> dict[str, list[str]]:
        """Return the external accessions of each entity from struct_ref.

        Returns:
            {entity_id: [accessions]}.
        """
        refs: dict[str, list[str]] = {}
        for raw_entity, raw_accession in zip(self._field("_struct_ref.entity_id"),
                                             self._field("_struct_ref.pdbx_db_accession")):
            entity, accession = clean_cif_value(raw_entity), clean_cif_value(raw_accession)
            if entity and accession and accession not in refs.setdefault(entity, []):
                refs[entity].append(accession)
        return refs

    def _entity_descriptions(self) -> dict[str, str]:
        """Return the description of each entity from the entity category.

        Returns:
            {entity_id: description}; entities without a description are left out.
        """
        descriptions: dict[str, str] = {}
        for raw_entity, raw_description in zip(self._field("_entity.id"), self._field("_entity.pdbx_description")):
            entity, description = clean_cif_value(raw_entity), clean_cif_value(raw_description)
            if entity and description:
                descriptions[entity] = description
        return descriptions

    def _label_id(self, chain_key: str) -> str:
        """Return the label chain ID of a chain key.

        Args:
            chain_key: Chain key as held by the structure.

        Returns:
            The label ID, or the key itself if unknown.
        """
        key = str(chain_key)
        return key if self.key_by == "label" else self._auth_to_label.get(key, key)

    def _auth_id(self, chain_key: str) -> str:
        """Return the author chain ID of a chain key.

        Args:
            chain_key: Chain key as held by the structure.

        Returns:
            The author ID, or the key itself if unknown.
        """
        key = str(chain_key)
        return key if self.key_by == "auth" else self._label_to_auth.get(key, key)

    def _entity_reference(self, chain_id: str) -> str:
        """Return the reference sequence of a chain's entity.

        Args:
            chain_id: Chain key as held by the structure.

        Returns:
            The entity sequence, or an empty string.
        """
        entity = self._entity_by_chain.get(str(chain_id))
        return self._reference_by_entity.get(entity, "") if entity else ""

    def _resolve_entity_id(self, chain_id: str, explicit: str | None = None) -> str | None:
        """Return the explicit entity ID, else the mmCIF entity, else the default entity ID.

        Args:
            chain_id: Chain key as held by the structure.
            explicit: Entity ID given by the caller, or None.

        Returns:
            The entity ID, or None if none is known.
        """
        if explicit is not None:
            return str(explicit)
        return self._entity_by_chain.get(str(chain_id)) or self.default_entity_id
