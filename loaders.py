#!/usr/bin/env python
"""Constructors that load structures from RCSB, local PDB or mmCIF files, or Biopython objects, keyed by author or label chain IDs."""

###################################### IMPORTS AND CONSTANTS

import io
import time
import urllib.request
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from constants import DEFAULT_CONFIG_PATH, RCSB_URL, SUPPORTED_SUFFIXES
from helpers import finish_log

KEY_MODES = ("auth", "label")

###################################### LOADERS


class LoaderMixin:
    @classmethod
    def from_rcsb_id(cls, pdb_id, config_path=DEFAULT_CONFIG_PATH, key_by="auth"):
        """Download an mmCIF structure from RCSB and return a SyncedStructure.

        Args:
            pdb_id: Four-character PDB identifier.
            config_path: Path to the YAML config file.
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
        finish_log(start, f"Loaded RCSB structure {pdb_id}")
        return obj

    @classmethod
    def from_Path(cls, path, config_path=DEFAULT_CONFIG_PATH, key_by="auth"):
        """Load a local PDB or mmCIF file; the file stem becomes the structure and entity ID.

        Args:
            path: Path to a .pdb, .cif, or .mmcif file.
            config_path: Path to the YAML config file.
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
        finish_log(start, f"Loaded local structure {path}")
        return obj

    @classmethod
    def from_path(cls, path, config_path=DEFAULT_CONFIG_PATH, key_by="auth"):
        """Alias for from_Path."""
        return cls.from_Path(path, config_path=config_path, key_by=key_by)

    @classmethod
    def from_biopython(cls, structure, structure_id=None, config_path=DEFAULT_CONFIG_PATH, key_by="auth"):
        """Wrap an existing Biopython Structure in a SyncedStructure.

        Args:
            structure: Biopython Structure object.
            structure_id: Optional ID; defaults to the structure's own ID.
            config_path: Path to the YAML config file.
            key_by: 'auth' (default) or 'label'.
        """
        return cls(
            structure,
            structure_id=structure_id or structure.id,
            file_format="pdb",
            config_path=config_path,
            key_by=key_by,
        )