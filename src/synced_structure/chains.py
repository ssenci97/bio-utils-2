"""Chain lifecycle, sequence mapping, chain renaming, external references, and per-chain summaries keyed by author IDs with standard IDs alongside."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

import copy
import itertools
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from ._types import ChainRefs, Moieties, Self
from .alignment import query_coverage
from .base import StructureState
from .constants import (
    ANGSTROM_PER_NM, CONTACTS_CUTOFF_NM, COORD_BLOCK_ROWS, DEFAULT_CONTACT_MOIETIES, MATCH_COLUMNS, MOIETY_OTHER,
    POLYMER_MOIETIES, RENAME_COLUMNS,
)
from .residues import split_residues
from .sequence import clean_sequence, one_letter_code
from .utils import log_info

if TYPE_CHECKING:
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue

###################################### CHAIN RECORD


@dataclass
class ChainRecord:
    """State of one loaded chain.

    Sequence indices are 0-based. Values of sequence_index_to_ref_index are 0-based reference indices.

    Attributes:
        chain: Biopython chain.
        chain_id: Label (standard) chain ID.
        auth_chain_id: Author chain ID.
        entity_id: Entity ID, or None if unknown.
        residues: Polymer residues in chain order.
        moieties: Moiety of each polymer residue.
        hetero: (residue, moiety) pairs of ligands, ions, and water.
        sequence: Observed one-letter sequence.
        reference_sequence: Reference one-letter sequence.
        metadata: User-defined values.
        external_ref: User-set external reference, or None.
        refseq_id: User-set reference sequence ID, or None.
        sequence_index_to_ref_index: Observed-to-reference index map.
    """

    chain: Chain
    chain_id: str
    auth_chain_id: str
    entity_id: str | None
    residues: list[Residue]
    moieties: list[str]
    hetero: list[tuple[Residue, str]]
    sequence: str
    reference_sequence: str
    metadata: dict[str, Any] = field(default_factory=dict)
    external_ref: str | None = None
    refseq_id: str | None = None
    sequence_index_to_ref_index: dict[int, int] = field(default_factory=dict)

    @property
    def sequence_length(self) -> int:
        """Length of the observed polymer sequence."""
        return len(self.sequence)

    @property
    def reference_length(self) -> int:
        """Length of the reference sequence."""
        return len(self.reference_sequence)

    @property
    def mapped_residues(self) -> int:
        """Number of observed residues with a reference position."""
        return len(self.sequence_index_to_ref_index)

###################################### CHAIN MANAGEMENT


class ChainMixin(StructureState):
    """Chain lifecycle, reference sequence mapping, chain renaming, and per-chain summaries."""

    def _chain_key(self, chain_id: str) -> str:
        """Return the chain key matching a key, standard ID, or author ID.

        Args:
            chain_id: Chain key, standard chain ID, or author chain ID.

        Returns:
            The chain key.

        Raises:
            KeyError: If no loaded chain matches.
        """
        chain_id = str(chain_id)
        if chain_id in self.chains:
            return chain_id
        for attribute in ("chain_id", "auth_chain_id"):
            for key, data in self.chains.items():
                if getattr(data, attribute) == chain_id:
                    return key
        raise KeyError(f"Unknown chain: {chain_id}")

    def _reference_index_map(self, sequence: str, reference: str) -> dict[int, int]:
        """Map observed to reference positions with a global alignment; the only such mapping in the package.

        Args:
            sequence: Observed sequence.
            reference: Reference sequence.

        Returns:
            {observed 0-based index: reference 0-based index}; public tables report them 1-based.
        """
        if not sequence or not reference:
            return {}
        alignment = self._global_aligner.align(reference, sequence)[0]
        return {
            int(seq_index): int(ref_index)
            for ref_index, seq_index in alignment.indices.T
            if ref_index >= 0 and seq_index >= 0
        }

    def _make_chain_record(self, chain: Chain, reference_sequence: str | None = None, entity_id: str | None = None,
                         metadata: Mapping[str, Any] | None = None, external_ref: str | None = None,
                         refseq_id: str | None = None) -> ChainRecord:
        """Build the record of a chain from its residues and metadata.

        Args:
            chain: Biopython chain, keyed by its chain key.
            reference_sequence: Reference sequence, or None to use mmCIF metadata.
            entity_id: Entity ID, or None to resolve it automatically.
            metadata: User metadata, or None.
            external_ref: User-set external reference, or None to use the mmCIF reference.
            refseq_id: User-set reference sequence ID, or None to fall back to the external reference.

        Returns:
            The chain record, with its observed-to-reference mapping.
        """
        key = str(chain.id)
        polymer, hetero = split_residues(chain)
        if reference_sequence is None:
            reference_sequence = self._entity_reference(key)
        sequence = "".join(one_letter_code(residue, moiety) for residue, moiety in polymer)
        reference_sequence = clean_sequence(reference_sequence)
        record = ChainRecord(
            chain=chain,
            chain_id=self._label_id(key),
            auth_chain_id=self._auth_id(key),
            entity_id=self._resolve_entity_id(key, entity_id),
            residues=[residue for residue, _ in polymer],
            moieties=[moiety for _, moiety in polymer],
            hetero=hetero,
            sequence=sequence,
            reference_sequence=reference_sequence,
            metadata=dict(metadata or {}),
            external_ref=external_ref,
            refseq_id=refseq_id,
        )
        record.sequence_index_to_ref_index = self._reference_index_map(sequence, reference_sequence)
        return record

    def _build_chains(self) -> None:
        """Rebuild chain records from the model, keeping the user-set values of existing chains."""
        previous = self.chains
        rebuilt = {}
        for chain in self.model.get_chains():
            key = str(chain.id)
            old = previous.get(key)
            kept: dict[str, Any] = {} if old is None else {
                "reference_sequence": old.reference_sequence,
                "entity_id": old.entity_id,
                "metadata": old.metadata,
                "external_ref": old.external_ref,
                "refseq_id": old.refseq_id,
            }
            rebuilt[key] = self._make_chain_record(chain, **kept)
        self.chains = rebuilt
        self._interface_cache.clear()

    @property
    def chain_ids(self) -> list[str]:
        """List of loaded chain keys."""
        return list(self.chains)

    def __getitem__(self, chain_id: str) -> ChainRecord:
        """Return the record of a chain.

        Args:
            chain_id: Chain key, standard ID, or author ID.

        Returns:
            The chain record.
        """
        return self.chains[self._chain_key(chain_id)]

    def __iter__(self) -> Iterator[ChainRecord]:
        """Iterate over chain records.

        Returns:
            An iterator over the chain records.
        """
        return iter(self.chains.values())

    def __len__(self) -> int:
        """Return the number of loaded chains.

        Returns:
            The number of chains.
        """
        return len(self.chains)

    def __contains__(self, chain_id: object) -> bool:
        """Return True if the chain is loaded.

        Args:
            chain_id: Chain key, standard ID, or author ID.

        Returns:
            Whether a loaded chain matches.
        """
        try:
            self._chain_key(str(chain_id))
            return True
        except KeyError:
            return False

    def get_sequence(self, chain_id: str) -> str:
        """Return the observed one-letter sequence of a chain.

        Args:
            chain_id: Chain identifier, standard or author.

        Returns:
            The observed sequence.
        """
        return self.chains[self._chain_key(chain_id)].sequence

    def get_reference_sequence(self, chain_id: str) -> str:
        """Return the reference sequence of a chain.

        Args:
            chain_id: Chain identifier, standard or author.

        Returns:
            The reference sequence.
        """
        return self.chains[self._chain_key(chain_id)].reference_sequence

    def get_accession(self, chain_id: str) -> str | None:
        """Return the external reference of a chain: the user-set value, else the mmCIF accession.

        Args:
            chain_id: Chain identifier, standard or author.

        Returns:
            The reference, or None if there is none.
        """
        data = self.chains[self._chain_key(chain_id)]
        if data.external_ref:
            return data.external_ref
        return ", ".join(self._external_refs.get(data.entity_id or "", [])) or None

    def set_accession(self, chain_id: str, external_ref: str | None) -> Self:
        """Set a chain's external reference, overriding the mmCIF value in summaries.

        Args:
            chain_id: Chain identifier, standard or author.
            external_ref: External accession, or None to fall back to the mmCIF reference.

        Returns:
            This structure.
        """
        chain_id = self._chain_key(chain_id)
        self.chains[chain_id].external_ref = None if external_ref is None else str(external_ref)
        return self

    def get_reference_id(self, chain_id: str) -> str | None:
        """Return the ID of a chain's reference sequence: the user-set value, else the external reference.

        Args:
            chain_id: Chain identifier, standard or author.

        Returns:
            The reference sequence ID, or None.
        """
        data = self.chains[self._chain_key(chain_id)]
        return data.refseq_id or self.get_accession(chain_id)

    def set_reference_id(self, chain_id: str, refseq_id: str | None) -> Self:
        """Set the ID of a chain's reference sequence, which may differ from its external reference.

        Args:
            chain_id: Chain identifier, standard or author.
            refseq_id: Reference sequence ID, or None to fall back to the external reference.

        Returns:
            This structure.
        """
        chain_id = self._chain_key(chain_id)
        self.chains[chain_id].refseq_id = None if refseq_id is None else str(refseq_id)
        return self

    def set_chain_metadata(self, chain_id: str, key: str, value: Any) -> Self:
        """Store a user-defined value on a chain.

        Args:
            chain_id: Chain identifier, standard or author.
            key: Metadata key.
            value: Value to store.

        Returns:
            This structure.
        """
        self.chains[self._chain_key(chain_id)].metadata[str(key)] = value
        return self

    def set_entity_id(self, chain_id: str, entity_id: str) -> Self:
        """Set a chain's entity ID.

        Args:
            chain_id: Chain identifier, standard or author.
            entity_id: New entity identifier.

        Returns:
            This structure.
        """
        self.chains[self._chain_key(chain_id)].entity_id = str(entity_id)
        return self

    def set_reference_sequence(self, chain_id: str, sequence: str, remap: bool = True) -> Self:
        """Set a chain's reference sequence, optionally rebuilding its mapping.

        Args:
            chain_id: Chain identifier, standard or author.
            sequence: New reference sequence.
            remap: Whether to rebuild the mapping.

        Returns:
            This structure.
        """
        chain_id = self._chain_key(chain_id)
        data = self.chains[chain_id]
        data.reference_sequence = clean_sequence(sequence)
        if remap:
            data.sequence_index_to_ref_index = self._reference_index_map(data.sequence, data.reference_sequence)
        else:
            data.sequence_index_to_ref_index = {}
            log_info(f"{self.structure_id}: {chain_id} has no residue mapping until remapped.")
        return self

    def set_reference_sequences(self, sequences: Mapping[str, str], remap: bool = True) -> Self:
        """Set the reference sequences of several chains.

        Args:
            sequences: {chain identifier: reference sequence}.
            remap: Whether to rebuild the mappings.

        Returns:
            This structure.
        """
        for chain_id, sequence in sequences.items():
            self.set_reference_sequence(chain_id, sequence, remap=remap)
        return self

    def remap_chain(self, chain_id: str) -> Self:
        """Rebuild the mapping of one chain against its current reference.

        Args:
            chain_id: Chain identifier, standard or author.

        Returns:
            This structure.
        """
        chain_id = self._chain_key(chain_id)
        self.set_reference_sequence(chain_id, self.chains[chain_id].reference_sequence, remap=True)
        log_info(f"{self.structure_id}: remapped {chain_id}")
        return self

    def remap_all_chains(self) -> Self:
        """Rebuild the mappings of all chains.

        Returns:
            This structure.
        """
        for chain_id in list(self.chains):
            self.remap_chain(chain_id)
        return self

    def _apply_chain_options(self, chain_id: str, reference_sequence: str | None, entity_id: str | None,
                             metadata: Mapping[str, Any] | None) -> None:
        """Apply optional reference, entity, and metadata values to a chain.

        Args:
            chain_id: Chain identifier, standard or author.
            reference_sequence: Reference sequence, or None.
            entity_id: Entity ID, or None.
            metadata: Metadata values, or None.
        """
        chain_id = self._chain_key(chain_id)
        if reference_sequence is not None:
            self.set_reference_sequence(chain_id, reference_sequence)
        if entity_id is not None:
            self.set_entity_id(chain_id, entity_id)
        for key, value in (metadata or {}).items():
            self.set_chain_metadata(chain_id, key, value)

    def add_chain(self, chain: Chain, reference_sequence: str | None = None, entity_id: str | None = None,
                  metadata: Mapping[str, Any] | None = None, replace: bool = False) -> Self:
        """Add a copy of a Biopython chain, optionally replacing an existing one.

        Args:
            chain: Biopython chain.
            reference_sequence: Optional reference sequence.
            entity_id: Optional entity ID.
            metadata: Optional metadata values.
            replace: Whether to replace an existing chain with the same ID.

        Returns:
            This structure.

        Raises:
            ValueError: If the chain exists and replace is False.
        """
        key = str(chain.id)
        try:
            existing: str | None = self._chain_key(key)
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

    def remove_chain(self, chain_id: str) -> Self:
        """Remove a chain.

        Args:
            chain_id: Chain identifier, standard or author.

        Returns:
            This structure.
        """
        chain_id = self._chain_key(chain_id)
        self.model.detach_child(chain_id)
        self.chains.pop(chain_id, None)
        self._alignment_cache.pop(chain_id, None)
        self._interface_cache.clear()
        return self

    def subset_chains(self, chain_ids: ChainRefs, structure_id: str | None = None) -> Self:
        """Return a copy holding only the given chains, with their references, IDs, and metadata.

        The parent is unchanged; the mmCIF dict and aligners are shared with it, not copied.

        Args:
            chain_ids: Chain identifiers to keep, standard or author.
            structure_id: ID of the copy, or None for '{structure_id}_{key}_{key}...'.

        Returns:
            The new structure.

        Raises:
            ValueError: If no chain is given.
        """
        keep = [self._chain_key(chain_id) for chain_id in chain_ids]
        if not keep:
            raise ValueError("subset_chains needs at least one chain.")
        shared = (self.mmcif_dict, self._global_aligner, self._local_aligner)
        child = copy.deepcopy(self, memo={id(obj): obj for obj in shared})
        for key in self.chains:
            if key not in keep:
                child.remove_chain(key)
        child.structure_id = str(structure_id or "_".join([self.structure_id, *dict.fromkeys(keep)]))
        return child

    def replace_chain(self, chain_id: str, chain: Chain, reference_sequence: str | None = None,
                      entity_id: str | None = None, metadata: Mapping[str, Any] | None = None) -> Self:
        """Replace a chain with a copy of a Biopython chain, keeping its key.

        Args:
            chain_id: Chain identifier, standard or author.
            chain: Replacement Biopython chain.
            reference_sequence: Optional reference sequence.
            entity_id: Optional entity ID.
            metadata: Optional metadata values.

        Returns:
            This structure.
        """
        chain_id = self._chain_key(chain_id)
        new_chain = copy.deepcopy(chain)
        new_chain.id = chain_id
        self.model.detach_child(chain_id)
        self.model.add(new_chain)
        self._build_chains()
        self._apply_chain_options(chain_id, reference_sequence, entity_id, metadata)
        return self

###################################### SEQUENCE MATCHING AND RENAMING

    def _local_identity(self, seq_a: str, seq_b: str) -> tuple[float, int, float]:
        """Compare two sequences with a local alignment.

        Identity is matches over columns where both sequences have a residue, gaps excluded.

        Args:
            seq_a: First sequence.
            seq_b: Second sequence.

        Returns:
            (identity percent, aligned residues, percent of seq_b covered).
        """
        if not seq_a or not seq_b:
            return 0.0, 0, 0.0
        alignments = self._local_aligner.align(seq_a, seq_b)
        if not alignments:
            return 0.0, 0, 0.0
        top, bottom = alignments[0]
        columns = [(char_a, char_b) for char_a, char_b in zip(top, bottom) if char_a != "-" and char_b != "-"]
        matches = sum(char_a == char_b for char_a, char_b in columns)
        aligned = len(columns)
        identity = round(matches / aligned * 100, 2) if aligned else 0.0
        return identity, aligned, query_coverage(aligned, len(seq_b))

    def find_chains_by_sequence(self, sequence: str, min_identity: float = 90.0,
                                min_coverage: float = 50.0) -> pd.DataFrame:
        """Find chains whose observed sequence aligns to a query sequence.

        Args:
            sequence: Query sequence.
            min_identity: Minimum identity percent to report.
            min_coverage: Minimum percent of the query that must align, so a few matching letters are not a match.

        Returns:
            One row per matching chain, best identity first, then most aligned residues.
        """
        query = clean_sequence(sequence)
        rows = []
        for data in self.chains.values():
            identity, aligned, coverage = self._local_identity(data.sequence, query)
            if aligned and identity >= min_identity and coverage >= min_coverage:
                rows.append({
                    "chain_id": data.chain_id,
                    "auth_chain_id": data.auth_chain_id,
                    "identity_percent": identity,
                    "aligned_residues": aligned,
                    "query_coverage_percent": coverage,
                })
        frame = pd.DataFrame(rows, columns=MATCH_COLUMNS)
        return frame.sort_values(["identity_percent", "aligned_residues"], ascending=False, ignore_index=True)

    def rename_chains(self, mapping: Mapping[str, str]) -> Self:
        """Rename chain keys, keeping standard and author ID maps, entities, and sequences in sync.

        Args:
            mapping: {old chain identifier: new key}.

        Returns:
            This structure.

        Raises:
            ValueError: If renaming would create duplicate chain IDs.
        """
        resolved = {self._chain_key(old): str(new) for old, new in mapping.items()}
        final_ids = [resolved.get(key, key) for key in self.chains]
        if len(set(final_ids)) != len(final_ids):
            raise ValueError("Renaming would create duplicate chain IDs.")
        moved = {old: self.model[old] for old in resolved}
        for old in resolved:
            self.model.detach_child(old)
        key_to_other = self._auth_to_label if self.key_by == "auth" else self._label_to_auth
        other_to_key = self._label_to_auth if self.key_by == "auth" else self._auth_to_label
        for old, new in resolved.items():
            moved[old].id = new
            self.model.add(moved[old])
            other = key_to_other.pop(old, None)
            if other is not None:
                key_to_other[new] = other
                other_to_key[other] = new
            entity = self._entity_by_chain.pop(old, None)
            if entity is not None:
                self._entity_by_chain[new] = entity
        self.chains = {resolved.get(key, key): data for key, data in self.chains.items()}
        self._alignment_cache.clear()
        self._build_chains()
        return self

    def _match_targets(self, targets: Mapping[str, str], min_identity: float,
                       min_coverage: float) -> tuple[dict[str, str], dict[str, tuple[str, float]]]:
        """Pair each target with the most similar chain not yet used, greedily by identity, then aligned length.

        Args:
            targets: {new key: cleaned sequence}.
            min_identity: Minimum identity percent required for a pair.
            min_coverage: Minimum percent of the target sequence that must align.

        Returns:
            (renames, matched): old key to new key, and new key to (old key, identity).
        """
        candidates = []
        for new_id, sequence in targets.items():
            for old_id, record in self.chains.items():
                identity, aligned, coverage = self._local_identity(record.sequence, sequence)
                if coverage >= min_coverage:
                    candidates.append((identity, aligned, old_id, new_id))
        candidates.sort(key=lambda row: (-row[0], -row[1]))
        renames: dict[str, str] = {}
        matched: dict[str, tuple[str, float]] = {}
        used_new: set[str] = set()
        for identity, _, old_id, new_id in candidates:
            if new_id in used_new or old_id in renames or identity < min_identity:
                continue
            renames[old_id] = new_id
            matched[new_id] = (old_id, identity)
            used_new.add(new_id)
        return renames, matched

    def rename_chains_by_sequence(self, targets: Mapping[str, str], min_identity: float = 90.0,
                                  min_coverage: float = 50.0) -> pd.DataFrame:
        """Rename chain keys to the target whose sequence they match best.

        Args:
            targets: {new key: sequence}.
            min_identity: Minimum identity percent required for a rename.
            min_coverage: Minimum percent of the target sequence that must align, so a few matching letters are
                not a match.

        Returns:
            One row per target and unmatched chain, with the old standard and author IDs and a status.
        """
        targets = {str(new): clean_sequence(seq) for new, seq in targets.items()}
        renames, matched = self._match_targets(targets, min_identity, min_coverage)
        rows = []
        for new_id in targets:
            if new_id in matched:
                old_id, identity = matched[new_id]
                old = self.chains[old_id]
                rows.append({"new_chain_id": new_id, "old_chain_id": old.chain_id,
                             "old_auth_chain_id": old.auth_chain_id, "identity_percent": identity, "status": "renamed"})
            else:
                rows.append({"new_chain_id": new_id, "old_chain_id": None, "old_auth_chain_id": None,
                             "identity_percent": None, "status": "unmatched_target"})
        for key, record in self.chains.items():
            if key not in renames:
                rows.append({"new_chain_id": None, "old_chain_id": record.chain_id,
                             "old_auth_chain_id": record.auth_chain_id, "identity_percent": None,
                             "status": "unmatched_chain"})
        if renames:
            self.rename_chains(renames)
            for new_id in renames.values():
                self.set_reference_sequence(new_id, targets[new_id])
        return pd.DataFrame(rows, columns=RENAME_COLUMNS)

###################################### CONTACT PARTNERS

    @staticmethod
    def _coords_touch(coords_a: np.ndarray, coords_b: np.ndarray, cutoff: float,
                      block_rows: int = COORD_BLOCK_ROWS) -> bool:
        """Return True if any point of coords_a lies within cutoff of any point of coords_b.

        Args:
            coords_a: Coordinates of the first set.
            coords_b: Coordinates of the second set, in the same unit.
            cutoff: Distance cutoff, in the same unit.
            block_rows: Number of rows of coords_a compared at once.

        Returns:
            Whether the sets touch.
        """
        if (coords_a.min(axis=0) > coords_b.max(axis=0) + cutoff).any() or (
                coords_b.min(axis=0) > coords_a.max(axis=0) + cutoff).any():
            return False
        limit = cutoff * cutoff
        for start in range(0, len(coords_a), block_rows):
            if (((coords_a[start:start + block_rows, None] - coords_b[None]) ** 2).sum(axis=-1) <= limit).any():
                return True
        return False

    def chain_partner_map(self, cutoff: float = CONTACTS_CUTOFF_NM,
                               moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> dict[str, list[str]]:
        """Find which chains touch, from inter-chain contact-point distances.

        Args:
            cutoff: Contact cutoff in nm.
            moieties: Biomolecule types whose contact points are compared.

        Returns:
            {chain key: sorted keys of the chains in contact}.
        """
        points: dict[str, list[np.ndarray]] = {}
        for record in self._contact_records(moieties=moieties):
            points.setdefault(record["key"], []).append(record["point"])
        coords = {key: np.asarray(values, dtype=np.float64) / ANGSTROM_PER_NM for key, values in points.items()}
        partners: dict[str, list[str]] = {key: [] for key in self.chains}
        for chain_a, chain_b in itertools.combinations(coords, 2):
            if self._coords_touch(coords[chain_a], coords[chain_b], cutoff):
                partners[chain_a].append(chain_b)
                partners[chain_b].append(chain_a)
        return {key: sorted(values) for key, values in partners.items()}

###################################### SUMMARY

    def chain_table(self, partners: bool = False, cutoff: float = CONTACTS_CUTOFF_NM,
                      moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> pd.DataFrame:
        """Summarize every chain: structure and chain IDs, entity, polymer moiety, description, references, and mapping coverage.

        The moiety column names the polymer type, 'protein' or 'nucleic' ('+'-joined when mixed), or 'other' when
        the chain has no polymer; ligands, ions, and water are counted only in moiety_table.

        Args:
            partners: Whether to add a 'partners' column of ';'-joined contact partner keys; computes contacts.
            cutoff: Contact cutoff in nm, used only with partners=True.
            moieties: Biomolecule types whose contact points are compared, used only with partners=True.

        Returns:
            One row per chain.
        """
        contacts = self.chain_partner_map(cutoff=cutoff, moieties=moieties) if partners else None
        descriptions = self._entity_descriptions()
        rows = []
        for key, data in self.chains.items():
            length = data.sequence_length
            present = set(data.moieties)
            row: dict[str, Any] = {
                "structure_id": self.structure_id,
                "chain_id": data.chain_id,
                "auth_chain_id": data.auth_chain_id,
                "entity_id": data.entity_id,
                "moiety": "+".join(moiety for moiety in POLYMER_MOIETIES if moiety in present) or MOIETY_OTHER,
                "description": descriptions.get(data.entity_id or ""),
                "external_ref": self.get_accession(key),
                "refseq_id": self.get_reference_id(key),
                "sequence_length": length,
                "reference_length": data.reference_length,
                "mapped_residues": data.mapped_residues,
                "mapping_fraction": round(data.mapped_residues / length, 3) if length else 0.0,
            }
            if contacts is not None:
                row["partners"] = ";".join(contacts.get(key, []))
            rows.append(row)
        return pd.DataFrame(rows)
