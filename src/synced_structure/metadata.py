"""Reference sequence, entity, external reference, and description lookups from mmCIF metadata, resolving label IDs for each chain key."""

###################################### IMPORTS AND CONSTANTS

from collections import Counter

from .sequence import cif_sequence, clean_val, to_list

###################################### METADATA


class MetadataMixin:
    """Lookups from mmCIF metadata: entities, label and author ID maps, external references, and descriptions."""

    def _field(self, name):
        """Return an mmCIF field as a list of strings."""
        return to_list(self.mmcif_dict.get(name))

    def _entity_sequences(self):
        """Return {entity_id: cleaned reference sequence} for entities with a non-empty sequence."""
        canonical = self._field("_entity_poly.pdbx_seq_one_letter_code_can")
        original = self._field("_entity_poly.pdbx_seq_one_letter_code")
        sequences = {}
        for index, entity_id in enumerate(self._field("_entity_poly.entity_id")):
            raw = canonical[index] if index < len(canonical) else ""
            if not raw and index < len(original):
                raw = original[index]
            sequence = cif_sequence(raw) if raw else ""
            if sequence:
                sequences[entity_id] = sequence
        return sequences

    def _atom_pair_maps(self):
        """Build label-to-author and author-to-label maps from atom_site, preferring polymer atoms.

        Each author ID maps to the label ID that occurs most often among its polymer atoms, so a
        single author ID resolves to one label ID. Hetero atoms only fill in IDs with no polymer atoms.
        """
        labels = self._field("_atom_site.label_asym_id")
        auths = self._field("_atom_site.auth_asym_id")
        groups = self._field("_atom_site.group_PDB") or ["ATOM"] * len(labels)
        polymer, hetero = Counter(), Counter()
        for group, label, auth in zip(groups, labels, auths):
            label, auth = clean_val(label), clean_val(auth)
            if label and auth:
                (polymer if group == "ATOM" else hetero)[(label, auth)] += 1
        self._auth_to_label, self._label_to_auth = {}, {}
        for counts in (polymer, hetero):
            for (label, auth), _ in counts.most_common():
                self._auth_to_label.setdefault(auth, label)
                self._label_to_auth.setdefault(label, auth)

    def _chain_entity_map(self):
        """Return {chain_key: entity_id} for every chain key, resolving each key through its label ID."""
        label_entity = {
            clean_val(label): clean_val(entity)
            for label, entity in zip(self._field("_struct_asym.id"), self._field("_struct_asym.entity_id"))
            if clean_val(label) and clean_val(entity)
        }
        if self.key_by == "auth":
            pairs = self._auth_to_label.items()
        else:
            pairs = ((key, key) for key in self._label_to_auth)
        return {key: label_entity[label] for key, label in pairs if label in label_entity}

    def _entity_external_refs(self):
        """Return {entity_id: [external accessions]} from struct_ref."""
        refs = {}
        for entity, accession in zip(self._field("_struct_ref.entity_id"), self._field("_struct_ref.pdbx_db_accession")):
            entity, accession = clean_val(entity), clean_val(accession)
            if entity and accession and accession not in refs.setdefault(entity, []):
                refs[entity].append(accession)
        return refs

    def _entity_descriptions(self):
        """Return {entity_id: description} from the entity category."""
        return {
            clean_val(entity): clean_val(desc)
            for entity, desc in zip(self._field("_entity.id"), self._field("_entity.pdbx_description"))
            if clean_val(entity)
        }

    def _label_of(self, chain_key):
        """Return the label chain ID for a chain key, or the key itself if unknown."""
        key = str(chain_key)
        return key if self.key_by == "label" else self._auth_to_label.get(key, key)

    def _auth_of(self, chain_key):
        """Return the author chain ID for a chain key, or the key itself if unknown."""
        key = str(chain_key)
        return key if self.key_by == "auth" else self._label_to_auth.get(key, key)

    def _reference_for_chain(self, chain_id):
        """Return the reference sequence of a chain's entity, or an empty string."""
        return self._reference_by_entity.get(self._entity_by_chain.get(str(chain_id)), "")

    def _resolve_entity(self, chain_id, explicit=None):
        """Return the explicit entity ID, else the mmCIF entity, else the default entity ID.

        Args:
            chain_id: Chain key as held by the structure.
            explicit: Entity ID given by the caller, or None.
        """
        if explicit is not None:
            return str(explicit)
        return self._entity_by_chain.get(str(chain_id)) or self.default_entity_id
