"""Shared state of SyncedStructure, declared once so that every mixin is type-checked against it."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ._types import ChainPair, ChainRefs, Moieties, PathLike, Record, Self

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    import pandas as pd
    from Bio.Align import PairwiseAligner
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure
    from IPython.display import HTML

    from .chains import ChainRecord

###################################### SHARED STATE


class StructureState:
    """Attributes set by SyncedStructure.__init__ and the methods mixins call on each other.

    The method stubs exist only for type checkers; at runtime each name resolves to the mixin that defines it.

    Attributes:
        structure: Underlying Biopython structure.
        model: First model of the structure, the one analysed.
        structure_id: Identifier used in tables and output names.
        file_format: Source format, 'cif' or 'pdb'.
        mmcif_dict: Parsed mmCIF dictionary, empty for PDB input.
        default_entity_id: Entity ID used when mmCIF metadata is missing.
        config: Resolved configuration values.
        key_by: 'auth' or 'label', the chain ID used as key.
        chains: Chain records by chain key.
    """

    structure: Structure
    model: Model
    structure_id: str
    file_format: str
    mmcif_dict: dict[str, Any]
    default_entity_id: str | None
    config: dict[str, Any]
    key_by: str
    chains: dict[str, ChainRecord]
    _global_aligner: PairwiseAligner
    _local_aligner: PairwiseAligner
    _alignment_cache: dict[str, Any]
    _interface_cache: dict[tuple[Any, ...], Any]
    _auth_to_label: dict[str, str]
    _label_to_auth: dict[str, str]
    _reference_by_entity: dict[str, str]
    _entity_by_chain: dict[str, str]
    _external_refs: dict[str, list[str]]

    if TYPE_CHECKING:
        # chains.py
        def _require_chain(self, chain_id: str) -> str: ...
        def _build_chains(self) -> None: ...
        @property
        def chain_ids(self) -> list[str]: ...
        def get_sequence(self, chain_id: str) -> str: ...
        def get_reference_sequence(self, chain_id: str) -> str: ...
        def subset_chains(self, chain_ids: ChainRefs, structure_id: str | None = ...) -> Self: ...
        def chain_contact_partners(self, cutoff: float = ..., moieties: Moieties = ...) -> dict[str, list[str]]: ...

        # metadata.py
        def _entity_descriptions(self) -> dict[str, str]: ...
        def _label_of(self, chain_key: str) -> str: ...
        def _auth_of(self, chain_key: str) -> str: ...
        def _reference_for_chain(self, chain_id: str) -> str: ...
        def _resolve_entity(self, chain_id: str, explicit: str | None = ...) -> str | None: ...

        # residues.py and alignment.py
        def get_residue_mapping(self, chain_id: str) -> pd.DataFrame: ...
        @property
        def residue_table(self) -> pd.DataFrame: ...
        @property
        def moiety_table(self) -> pd.DataFrame: ...
        @property
        def reference_alignment_stats(self) -> pd.DataFrame: ...
        @property
        def reference_alignment_positions(self) -> pd.DataFrame: ...

        # contacts.py and contact_matrix.py
        def _contact_records(self, chain_ids: ChainRefs | None = ..., moieties: Moieties = ...) -> list[Record]: ...
        def compute_contacts(self, cutoff: float = ..., contact_type: str = ..., chain_ids: ChainRefs | None = ...,
                             scheme: str = ..., moieties: Moieties = ...) -> pd.DataFrame: ...
        def _resolve_pair(self, chain_a: str, chain_b: str | None) -> tuple[str, str, str, str]: ...
        @staticmethod
        def _pair_rows(contacts: pd.DataFrame, label_a: str, label_b: str) -> pd.DataFrame: ...
        def _block_matrix(self, keys: Sequence[str], cutoff: float, scheme: str, moieties: Moieties) -> Record: ...
        def contact_matrix(self, chain_a: str = ..., chain_b: str | None = ..., cutoff: float = ...,
                           moieties: Moieties = ..., scheme: str = ...) -> Record: ...

        # interface_metrics.py, exporter.py, viz3d.py
        def interface_metrics_all(self, chain_pairs: Sequence[ChainPair] | None = ...) -> dict[str, pd.DataFrame]: ...
        @staticmethod
        def _prepare_path(path: PathLike) -> Path: ...
        def show_structure_3d(self, angles: Sequence[float] = ..., chain_id: str | None = ...,
                              colors: dict[str, str] | None = ..., title: str | None = ..., style: str = ...,
                              background: str = ..., size: int = ..., show_angles: bool = ...,
                              chain_palette: Sequence[str] = ..., viz_format: str = ..., js_url: str = ...,
                              label_js: str = ..., angles_js: str = ..., title_font_px: int = ...,
                              angles_font_px: int = ...) -> HTML: ...
