"""SyncedStructure class assembled from the mixin modules, keyed by author chain IDs by default."""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

from typing import Any

from Bio.PDB.Structure import Structure

from ._types import PathLike
from .alignment import AlignmentMixin, make_global_aligner, make_local_aligner
from .chains import ChainMixin
from .constants import CONFIG_SECTION, DEFAULT_CONFIG, KEY_MODES
from .contact_matrix import ContactMatrixMixin
from .contact_plots import ContactPlotMixin
from .contact_stats import ContactStatsMixin
from .contacts import ContactMixin
from .docking import DockingMixin
from .exporter import ExportMixin
from .interface_metrics import InterfaceMixin
from .loaders import LoaderMixin
from .metadata import MetadataMixin
from .residues import ResidueMixin
from .utils import load_config
from .viz3d import Viz3DMixin

###################################### SYNCHRONIZED STRUCTURE


class SyncedStructure(LoaderMixin, MetadataMixin, ChainMixin, AlignmentMixin, ResidueMixin, ContactMixin,
                      ContactMatrixMixin, ContactStatsMixin, ContactPlotMixin, ExportMixin, InterfaceMixin, Viz3DMixin,
                      DockingMixin):
    """A structure whose chains are held by author or label IDs, with mapping, contacts, interfaces, and plots."""

    def __init__(self, structure: Structure, structure_id: str = "structure", mmcif_dict: dict[str, Any] | None = None,
                 file_format: str = "cif", config_path: PathLike | None = None, default_entity_id: str | None = None,
                 key_by: str = "auth") -> None:
        """Build the synchronized structure from a Biopython Structure and its optional YAML config.

        Args:
            structure: Biopython Structure.
            structure_id: Identifier used in output names.
            mmcif_dict: Parsed mmCIF dict providing reference metadata.
            file_format: Source format, 'cif' or 'pdb'.
            config_path: Path to the YAML config file, or None for defaults; never written.
            default_entity_id: Entity ID used when mmCIF metadata is missing.
            key_by: 'auth' (default) keys chains by author IDs, 'label' by label IDs.

        Raises:
            ValueError: If key_by is not 'auth' or 'label'.
        """
        if str(key_by).lower() not in KEY_MODES:
            raise ValueError(f"key_by must be one of {KEY_MODES}.")
        self.key_by = str(key_by).lower()
        self.structure = structure
        self.structure_id = str(structure_id)
        self.file_format = str(file_format).lower()
        self.mmcif_dict = mmcif_dict or {}
        self.default_entity_id = None if default_entity_id is None else str(default_entity_id)
        self.config = load_config(config_path, DEFAULT_CONFIG, CONFIG_SECTION)
        self.model = next(structure.get_models())
        self.chains = {}
        self._global_aligner = make_global_aligner()
        self._local_aligner = make_local_aligner(self.config)
        self._alignment_cache = {}
        self._interface_cache = {}
        self._build_chain_id_maps()
        self._reference_by_entity = self._entity_sequences()
        self._entity_by_chain = self._chain_entity_map()
        self._external_refs = self._entity_accessions()
        self._build_chains()
