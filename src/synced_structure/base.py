"""Shared state of SyncedStructure, declared once so that every mixin is type-checked against it."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ._types import ChainPair, ChainRefs, Moieties, Record, Self

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    import pandas as pd
    from Bio.Align import PairwiseAligner
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure
    from IPython.display import HTML

    from .chains import ChainRecord
    from .structure import SyncedStructure

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
        def _chain_key(self, chain_id: str) -> str: ...
        def _build_chains(self) -> None: ...
        @property
        def chain_ids(self) -> list[str]: ...
        def get_sequence(self, chain_id: str) -> str: ...
        def get_reference_sequence(self, chain_id: str) -> str: ...
        def get_accession(self, chain_id: str) -> str | None: ...
        def subset_chains(self, chain_ids: ChainRefs, structure_id: str | None = ...) -> Self: ...
        def chain_partner_map(self, cutoff: float = ..., moieties: Moieties = ...) -> dict[str, list[str]]: ...
        def chain_table(self, partners: bool = ..., cutoff: float = ..., moieties: Moieties = ...) -> pd.DataFrame: ...

        # metadata.py
        def _entity_descriptions(self) -> dict[str, str]: ...
        def _label_id(self, chain_key: str) -> str: ...
        def _auth_id(self, chain_key: str) -> str: ...
        def _entity_reference(self, chain_id: str) -> str: ...
        def _resolve_entity_id(self, chain_id: str, explicit: str | None = ...) -> str | None: ...

        # exporter.py
        def _mmcif_text(self, chain_ids: ChainRefs | None = ..., key_as_author: bool = ...) -> str: ...

        # residues.py and alignment.py
        def chain_residue_table(self, chain_id: str) -> pd.DataFrame: ...
        def residue_table(self) -> pd.DataFrame: ...
        def moiety_table(self) -> pd.DataFrame: ...
        def alignment_table(self) -> pd.DataFrame: ...
        def alignment_column_table(self) -> pd.DataFrame: ...

        # contacts.py and contact_matrix.py
        def _contact_records(self, chain_ids: ChainRefs | None = ..., moieties: Moieties = ...) -> list[Record]: ...
        def contact_table(self, cutoff: float = ..., contact_type: str = ..., chain_ids: ChainRefs | None = ...,
                             scheme: str = ..., moieties: Moieties = ...) -> pd.DataFrame: ...
        def _resolve_pair(self, chain_a: str, chain_b: str | None) -> tuple[str, str, str, str]: ...
        @staticmethod
        def _pair_contacts(contacts: pd.DataFrame, key_a: str, key_b: str) -> pd.DataFrame: ...
        def _block_matrix(self, keys: Sequence[str], cutoff: float, scheme: str, moieties: Moieties) -> Record: ...
        def _moiety_sequence(self, chain_id: str, moieties: Moieties = ...) -> str: ...
        def contact_string(self, chain_id: str, cutoff: float = ..., contact_type: str = ...,
                           contacts: pd.DataFrame | None = ..., scheme: str = ..., moieties: Moieties = ...) -> str: ...
        def reference_contact_string(self, chain_id: str, cutoff: float = ..., contact_type: str = ...,
                                     contacts: pd.DataFrame | None = ..., scheme: str = ...,
                                     moieties: Moieties = ...) -> str: ...
        def contact_stat_tables(self, cutoff: float = ..., moieties: Moieties = ...,
                                scheme: str = ...) -> dict[str, pd.DataFrame]: ...
        def contact_matrix(self, chain_a: str, chain_b: str | None = ..., cutoff: float = ...,
                           moieties: Moieties = ..., scheme: str = ...) -> Record: ...

        # interface_metrics.py, exporter.py, viz3d.py
        def interface_tables(self, chain_pairs: Sequence[ChainPair] | None = ...) -> dict[str, pd.DataFrame]: ...
        def show_structure_3d(self, angles: Sequence[float] = ..., chain_id: str | None = ...,
                              colors: dict[str, str] | None = ..., title: str | None = ..., style: str = ...,
                              background: str = ..., size: int = ..., show_angles: bool = ...,
                              chain_palette: Sequence[str] = ..., viz_format: str = ..., js_url: str = ...,
                              label_js: str = ..., angles_js: str = ..., title_font_px: int = ...,
                              angles_font_px: int = ..., highlight: dict[str, Sequence[int] | str] | None = ...,
                              dim_alpha: float = ..., highlight_sticks: bool = ...,
                              extra_models: Sequence[Record] | None = ...) -> HTML: ...
        def _chain_colors(self, colors: dict[str, str] | None = ...,
                          chain_palette: Sequence[str] = ...) -> dict[str, str]: ...

        # comparison.py
        def rmsd_metrics(self, reference: SyncedStructure, chains: str | ChainRefs | None = ...,
                         chain_map: Mapping[str, str] | None = ..., atoms: str = ...,
                         fit_chains: str | ChainRefs | None = ..., superpose: bool = ...) -> Record: ...
        def apply_superposition(self, reference: SyncedStructure, chains: str | ChainRefs | None = ...,
                                chain_map: Mapping[str, str] | None = ..., atoms: str = ...) -> Self: ...
