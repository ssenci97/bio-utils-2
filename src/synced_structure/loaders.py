"""Constructors that load structures from RCSB, local PDB or mmCIF files, or Biopython objects, keyed by author or label chain IDs."""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

import io
import time
import urllib.request
from pathlib import Path
from typing import TYPE_CHECKING, Any

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Structure import Structure

from ._types import PathLike, Self
from .base import StructureState
from .constants import RCSB_URL, SUPPORTED_SUFFIXES
from .utils import log_completion, log_info

###################################### LOADERS


class LoaderMixin(StructureState):
    """Constructors that load structures from RCSB, local files, or Biopython objects."""

    if TYPE_CHECKING:
        def __init__(self, structure: Structure, structure_id: str = ..., mmcif_dict: dict[str, Any] | None = ...,
                     file_format: str = ..., config_path: PathLike | None = ..., default_entity_id: str | None = ...,
                     key_by: str = ...) -> None:
            """Declare SyncedStructure's constructor for type checkers, so that cls(...) is checked.

            Args:
                structure: Biopython Structure.
                structure_id: Identifier used in output names.
                mmcif_dict: Parsed mmCIF dict, or None.
                file_format: Source format, 'cif' or 'pdb'.
                config_path: Path to the YAML config file, or None.
                default_entity_id: Entity ID used when mmCIF metadata is missing.
                key_by: 'auth' or 'label'.
            """

    @classmethod
    def from_rcsb(cls, pdb_id: str, config_path: PathLike | None = None, key_by: str = "auth") -> Self:
        """Download an mmCIF entry from RCSB.

        Args:
            pdb_id: Four-character PDB identifier.
            config_path: Path to the YAML config file, or None for defaults.
            key_by: 'auth' (default) keys chains by author IDs, 'label' by label IDs.

        Returns:
            The loaded structure.
        """
        start = time.perf_counter()
        pdb_id = str(pdb_id).strip().upper()
        with urllib.request.urlopen(RCSB_URL.format(pdb_id), timeout=30) as response:
            text = response.read().decode("utf-8")
        structure = MMCIFParser(QUIET=True, auth_chains=key_by == "auth").get_structure(pdb_id, io.StringIO(text))
        obj = cls(
            structure,
            structure_id=pdb_id,
            mmcif_dict=MMCIF2Dict(io.StringIO(text)),
            file_format="cif",
            config_path=config_path,
            key_by=key_by,
        )
        log_completion(start, f"Loaded RCSB structure {pdb_id}")
        return obj

    @classmethod
    def from_path(cls, path: PathLike, config_path: PathLike | None = None, key_by: str = "auth") -> Self:
        """Load a local PDB or mmCIF file; the file stem becomes the structure and entity ID.

        Args:
            path: Path to a .pdb, .cif, or .mmcif file.
            config_path: Path to the YAML config file, or None for defaults.
            key_by: 'auth' (default) keys chains by author IDs, 'label' by label IDs.

        Returns:
            The loaded structure.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file format is not supported.
        """
        start = time.perf_counter()
        file_path = Path(path).expanduser()
        if not file_path.is_file():
            raise FileNotFoundError(f"Structure file not found: {file_path}")
        suffix = file_path.suffix.lower()
        if suffix not in SUPPORTED_SUFFIXES:
            raise ValueError(f"Unsupported structure format: {suffix}. Supported formats: {sorted(SUPPORTED_SUFFIXES)}")
        structure_id = file_path.stem
        mmcif_dict: dict[str, Any]
        if suffix == ".pdb":
            if key_by == "label":
                log_info(f"{structure_id}: key_by='label' is not available for PDB files; chains keep their PDB IDs.")
            structure = PDBParser(QUIET=True).get_structure(structure_id, str(file_path))
            mmcif_dict, file_format = {}, "pdb"
        else:
            structure = MMCIFParser(QUIET=True, auth_chains=key_by == "auth").get_structure(structure_id, str(file_path))
            mmcif_dict, file_format = MMCIF2Dict(str(file_path)), "cif"
        obj = cls(
            structure,
            structure_id=structure_id,
            mmcif_dict=mmcif_dict,
            file_format=file_format,
            config_path=config_path,
            default_entity_id=structure_id,
            key_by=key_by,
        )
        log_completion(start, f"Loaded local structure {file_path}")
        return obj

    @classmethod
    def from_biopython(cls, structure: Structure, structure_id: str | None = None,
                       mmcif_dict: dict[str, Any] | None = None, file_format: str = "pdb",
                       config_path: PathLike | None = None, key_by: str = "auth") -> Self:
        """Wrap an existing Biopython structure; the structure ID becomes the default entity ID.

        Args:
            structure: Biopython Structure.
            structure_id: ID to use; defaults to the structure's own ID.
            mmcif_dict: Parsed mmCIF dict providing reference metadata, or None.
            file_format: Source format, 'cif' or 'pdb'.
            config_path: Path to the YAML config file, or None for defaults.
            key_by: 'auth' (default) or 'label'; must match how the structure was parsed.

        Returns:
            The wrapped structure.
        """
        structure_id = str(structure_id or structure.id)
        return cls(
            structure,
            structure_id=structure_id,
            mmcif_dict=mmcif_dict,
            file_format=file_format,
            config_path=config_path,
            default_entity_id=structure_id,
            key_by=key_by,
        )
