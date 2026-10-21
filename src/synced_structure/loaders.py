"""Constructors that load structures from RCSB, local PDB or mmCIF files, or Biopython objects, keyed by author or label chain IDs."""

###################################### IMPORTS AND CONSTANTS

import io
import time
import urllib.request
import warnings
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from .constants import RCSB_URL, SUPPORTED_SUFFIXES
from .utils import log_completion, log_info

###################################### LOADERS


class LoaderMixin:
    """Constructors that load structures from RCSB, local files, or Biopython objects."""

    @classmethod
    def from_rcsb(cls, pdb_id, config_path=None, key_by="auth"):
        """Download an mmCIF structure from RCSB and return a SyncedStructure.

        Args:
            pdb_id: Four-character PDB identifier.
            config_path: Path to the YAML config file, or None for defaults.
            key_by: 'auth' (default) keys chains by author IDs, 'label' by label IDs.
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
    def from_rcsb_id(cls, *args, **kwargs):
        """Deprecated alias of from_rcsb()."""
        warnings.warn("from_rcsb_id() is deprecated; use from_rcsb().", DeprecationWarning, stacklevel=2)
        return cls.from_rcsb(*args, **kwargs)

    @classmethod
    def from_path(cls, path, config_path=None, key_by="auth"):
        """Load a local PDB or mmCIF file; the file stem becomes the structure and entity ID.

        Args:
            path: Path to a .pdb, .cif, or .mmcif file.
            config_path: Path to the YAML config file, or None for defaults.
            key_by: 'auth' (default) keys chains by author IDs, 'label' by label IDs.
        """
        start = time.perf_counter()
        path = Path(path).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Structure file not found: {path}")
        suffix = path.suffix.lower()
        if suffix not in SUPPORTED_SUFFIXES:
            raise ValueError(f"Unsupported structure format: {suffix}. Supported formats: {sorted(SUPPORTED_SUFFIXES)}")
        structure_id = path.stem
        if suffix == ".pdb":
            if key_by == "label":
                log_info(f"{structure_id}: key_by='label' is not available for PDB files; chains keep their PDB IDs.")
            structure = PDBParser(QUIET=True).get_structure(structure_id, str(path))
            mmcif_dict, file_format = {}, "pdb"
        else:
            structure = MMCIFParser(QUIET=True, auth_chains=key_by == "auth").get_structure(structure_id, str(path))
            mmcif_dict, file_format = MMCIF2Dict(str(path)), "cif"
        obj = cls(
            structure,
            structure_id=structure_id,
            mmcif_dict=mmcif_dict,
            file_format=file_format,
            config_path=config_path,
            default_entity_id=structure_id,
            key_by=key_by,
        )
        log_completion(start, f"Loaded local structure {path}")
        return obj

    @classmethod
    def from_biopython(cls, structure, structure_id=None, mmcif_dict=None, file_format="pdb",
                       config_path=None, key_by="auth"):
        """Wrap an existing Biopython Structure in a SyncedStructure; the structure ID becomes the default entity ID.

        Args:
            structure: Biopython Structure object.
            structure_id: Optional ID; defaults to the structure's own ID.
            mmcif_dict: Optional parsed mmCIF dict providing reference metadata.
            file_format: Source format, 'cif' or 'pdb'.
            config_path: Path to the YAML config file, or None for defaults.
            key_by: 'auth' (default) or 'label'; must match how the structure was parsed.
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