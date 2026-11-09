"""HADDOCK3 docking set up from synchronized structures, run with the haddock3 command line in a subprocess.

Restraints are given in reference, sequence, or author numbering. Each docking partner is one rigid body and gets
one segment ID (A, B, ...): all of its chains are written under that ID, because haddock3 moves every segment ID as
a separate body. haddock_residue_map.tsv records how each written residue maps back to its chain, author number,
and sequence and reference positions.
"""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

import gzip
import importlib.util
import io
import itertools
import json
import math
import re
import shutil
import subprocess
import sys
import time
import warnings
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, TypedDict, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from ._types import ChainRef, PathLike, Record, Self
from .base import StructureState
from .chains import ChainMixin
from .constants import (
    CAPRI_LABELS, CONTACTS_CUTOFF_NM, EXPORT_SEP, HADDOCK_ALASCAN_COLUMNS, HADDOCK_ANALYSE_EXECUTABLE, HADDOCK_CHAIN_GAP, HADDOCK_CLEARANCE_A,
    HADDOCK_DEFAULT_SEED, HADDOCK_EXECUTABLE, HADDOCK_MAX_RESID, HADDOCK_METHODS, HADDOCK_NUMBERINGS,
    HADDOCK_DEFAULT_SAMPLING, HADDOCK_DEFAULT_SELECT, HADDOCK_DEFAULT_TOP_MODELS, HADDOCK_MIN_TOP_CLUSTERS,
    HADDOCK_PASSTHROUGH_MODULES, HADDOCK_PATCH_COLUMNS, HADDOCK_PATCH_GAP,
    HADDOCK_PLACEMENT_GROWTH, HADDOCK_PLACEMENT_TRIALS, HADDOCK_POLL_SECONDS, HADDOCK_PRESETS, HADDOCK_PROGRESS_SECONDS,
    HADDOCK_PROTECTED_GLOBALS, HADDOCK_REFINE_MODULES, HADDOCK_RESIDUE_MAP_COLUMNS, HADDOCK_RESTRAINT_COLUMNS, HADDOCK_RESTRAINT_DISTANCE,
    HADDOCK_RESTRAINT_KINDS, HADDOCK_RESTRAINT_MODULES, HADDOCK_SCORE_COLUMNS, HADDOCK_SEEDED_MODULES, HADDOCK_SEGIDS,
    HADDOCK_TBL_TERMS_PER_LINE, HADDOCK_TOP_CLUSTER_FRACTION, POLYMER_MOIETIES,
)
from .alignment import reference_index_map, sequence_index_map
from .plot_style import (
    CAPRI_ALPHA, CAPRI_AXIS_COLOR, CAPRI_CMAP, CAPRI_EDGE_COLOR, CAPRI_GRID_COLOR, CAPRI_LABEL_PT, CAPRI_MARKER_COLOR,
    CAPRI_MARKER_SIZE, CAPRI_SIZE_IN, CAPRI_TEXT_COLOR, CAPRI_TICK_PT, CAPRI_TITLE_PT,
)
from .residues import is_heavy_atom
from .utils import core_count, log_completion, log_info

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure, SubFigure

    from .structure import SyncedStructure

RUN_SENTINEL = ".synced_structure_run"
LOG_TAIL_LINES = 40
STEP_PATTERN = re.compile(r"^(\d+)_(\w+)$")
SCAN_CLUSTER_PATTERN = re.compile(r"^scan_clt_(.+)\.tsv$")
SCAN_MODEL_PATTERN = re.compile(r"^scan_(.+)\.tsv$")
# REMARK labels of HADDOCK model headers, beyond the score and energies lines, and their io.json names
REMARK_TERMS = {"Desolvation energy": "desolv", "buried surface area": "bsa"}

# a docking partner as given by the user, a chain as named in restraints, and a (partner index, chain key) pair
PartnerSpec = tuple["_structure.SyncedStructure", str | Sequence[str]]
ChainSpec = ChainRef | tuple[int, ChainRef]
ChainTarget = tuple[int, str]
# a reference given as a structure: the complex and {source chain key: (partner index, chain key)}
ReferenceSource = tuple["_structure.SyncedStructure", dict[str, ChainTarget]]
PatchPair = tuple[tuple[ChainSpec, Iterable[int]], tuple[ChainSpec, Iterable[int]]]
StepRef = str | int


class Partner(TypedDict):
    """One normalized docking partner.

    Attributes:
        structure: Structure holding the partner's chains.
        keys: Chain keys of the partner.
        segid: Segment ID written for all of its chains.
    """

    structure: _structure.SyncedStructure
    keys: list[str]
    segid: str


class ResidueEntry(TypedDict):
    """One written residue with its HADDOCK3 number and original positions.

    Attributes:
        resid: HADDOCK3 residue number.
        residue: Biopython residue.
        sequence_index: 0-based observed sequence index.
        reference_index: 0-based reference index, or None if unmapped.
    """

    resid: int
    residue: Residue
    sequence_index: int
    reference_index: int | None


class Restraint(TypedDict):
    """One restraint set between two partners, in HADDOCK3 residue numbers.

    Attributes:
        kind: 'ambiguous' or 'unambiguous'.
        origin: Where the restraint came from, such as 'custom' or 'contacts'.
        partner_a: Partner index of the first side.
        residues_a: Active (or patch) residues of the first side.
        passive_a: Passive residues of the first side.
        partner_b: Partner index of the second side.
        residues_b: Active (or patch) residues of the second side.
        passive_b: Passive residues of the second side.
        distance: (distance, lower correction, upper correction) in A.
    """

    kind: str
    origin: str
    partner_a: int
    residues_a: list[int]
    passive_a: list[int]
    partner_b: int
    residues_b: list[int]
    passive_b: list[int]
    distance: tuple[float, ...]


class Placement(TypedDict):
    """Rigid placement of a separated partner.

    Attributes:
        rotation: 3x3 rotation matrix.
        origin: Centroid the rotation is applied about.
        center: New centroid position.
    """

    rotation: np.ndarray
    origin: np.ndarray
    center: np.ndarray


class InputPaths(TypedDict):
    """Files written by HaddockDocking.export_inputs().

    Attributes:
        molecules: One PDB per partner; for 'alascan', one PDB of the complex.
        reference: Reference complex for caprieval.
        ambiguous: Ambiguous restraint file.
        unambiguous: Unambiguous restraint file.
        haddock_residue_map: Residue map table.
        config: haddock3 configuration.
        log: haddock3 log.
    """

    molecules: list[Path]
    reference: Path
    ambiguous: Path
    unambiguous: Path
    haddock_residue_map: Path
    config: Path
    log: Path

###################################### HELPERS


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Return a uniformly random 3x3 rotation matrix (Shoemake's method).

    Args:
        rng: Random number generator.

    Returns:
        Rotation matrix.
    """
    u1, u2, u3 = rng.random(3)
    x, y = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2), np.sqrt(1 - u1) * np.cos(2 * np.pi * u2)
    z, w = np.sqrt(u1) * np.sin(2 * np.pi * u3), np.sqrt(u1) * np.cos(2 * np.pi * u3)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def _random_direction(rng: np.random.Generator) -> np.ndarray:
    """Return a uniformly random unit vector.

    Args:
        rng: Random number generator.

    Returns:
        Unit vector of length 3.
    """
    vector = rng.normal(size=3)
    return vector / np.linalg.norm(vector)


def cluster_contact_patches(points: Sequence[tuple[float, float]], gap: float = HADDOCK_PATCH_GAP) -> list[list[int]]:
    """Group 2D points by single linkage, joining two points when both coordinates differ by at most gap.

    Args:
        points: (position_a, position_b) pairs, one per contact.
        gap: Largest difference along either axis that links two contacts.

    Returns:
        Clusters as sorted lists of point indices, ordered by their first index.
    """
    parent = list(range(len(points)))

    def root(index: int) -> int:
        """Return the cluster root of a point.

        Args:
            index: Point index.

        Returns:
            Index of the root point.
        """
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for i, j in itertools.combinations(range(len(points)), 2):
        if abs(points[i][0] - points[j][0]) <= gap and abs(points[i][1] - points[j][1]) <= gap:
            parent[root(i)] = root(j)
    clusters: dict[int, list[int]] = {}
    for index in range(len(points)):
        clusters.setdefault(root(index), []).append(index)
    return sorted(clusters.values(), key=lambda members: members[0])


def _toml_literal(value: Any) -> str:
    """Return a TOML literal for a bool, number, string, path, or list of these.

    Args:
        value: Value to write; numpy scalars are converted first.

    Returns:
        TOML literal.

    Raises:
        TypeError: If the value has no TOML form here.
    """
    if hasattr(value, "item") and not isinstance(value, (list, tuple)):
        value = value.item()
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, (str, Path)):
        return json.dumps(str(value))
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_toml_literal(item) for item in value) + "]"
    raise TypeError(f"Unsupported HADDOCK3 parameter value: {value!r}")


def _pdb_atom_line(serial: int, atom: Atom | DisorderedAtom, resname: str, segid: str, resid: int,
                   coord: np.ndarray) -> str:
    """Return one PDB ATOM record with the chain ID and segment ID both set to segid.

    Args:
        serial: Atom serial number.
        atom: Biopython atom; for alternate locations, the selected one.
        resname: Residue name.
        segid: One-character segment ID.
        resid: Residue number.
        coord: Coordinates in A.

    Returns:
        Fixed-width ATOM line without newline.
    """
    name = atom.get_id().strip()
    element = (atom.element or "".join(char for char in name if char.isalpha())[:1]).strip().upper()
    if len(name) < 4 and len(element) == 1:
        name = f" {name}"
    return (f"ATOM  {serial % 100000:>5} {name:<4} {resname:>3} {segid}{resid:>4}    "
            f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{1.0:6.2f}{atom.get_bfactor():6.2f}      "
            f"{segid:<4}{element:>2}")


def _resid_terms(resids: Iterable[int]) -> list[str]:
    """Return CNS residue terms, joining consecutive numbers as 'resi a:b'.

    Args:
        resids: Residue numbers, in any order.

    Returns:
        Terms such as 'resi 5' or 'resi 7:9', sorted.
    """
    terms: list[str] = []
    run: list[int] = []
    for resid in sorted(set(resids)):
        if run and resid != run[-1] + 1:
            terms.append(f"resi {run[0]}" if len(run) == 1 else f"resi {run[0]}:{run[-1]}")
            run = []
        run.append(resid)
    if run:
        terms.append(f"resi {run[0]}" if len(run) == 1 else f"resi {run[0]}:{run[-1]}")
    return terms


def _tbl_selection(segid: str, resids: Iterable[int], indent: str = "       ") -> str:
    """Return a CNS selection of residues in one segment, wrapped to stay under CNS's line limit.

    A long selection opens its residue list on the next line and ends each line with 'or': HADDOCK3's restraint
    validator miscounts a nested parenthesis left open on the selection's first line, and needs the closing line
    to start with a keyword.

    Args:
        segid: Segment ID.
        resids: Residue numbers.
        indent: Indentation of continuation lines.

    Returns:
        Parenthesized selection.
    """
    terms = _resid_terms(resids)
    if len(terms) == 1:
        return f"(segid {segid} and {terms[0]})"
    if len(terms) <= HADDOCK_TBL_TERMS_PER_LINE:
        return f"(segid {segid} and ({' or '.join(terms)}))"
    chunks = [" or ".join(terms[start:start + HADDOCK_TBL_TERMS_PER_LINE])
              for start in range(0, len(terms), HADDOCK_TBL_TERMS_PER_LINE)]
    body = f" or\n{indent}   ".join(chunks)
    return f"(segid {segid} and\n{indent}  ({body}))"


def _duration(seconds: float) -> str:
    """Format a duration such as '1h02m', '3m05s' or '42s'.

    Args:
        seconds: Duration in seconds.

    Returns:
        The formatted duration.
    """
    minutes, secs = divmod(round(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours}h{minutes:02d}m"
    return f"{minutes}m{secs:02d}s" if minutes else f"{secs}s"


def _log_tail(path: PathLike, lines: int = LOG_TAIL_LINES) -> str:
    """Return the last lines of a log file.

    Args:
        path: Log file.
        lines: Number of lines.

    Returns:
        The lines, newline-joined.
    """
    return "\n".join(Path(path).read_text(errors="replace").splitlines()[-lines:])


def _flatten(items: Any) -> Iterator[Any]:
    """Yield the leaves of nested lists, tuples, and dict values, as found in haddock3 io.json outputs.

    Args:
        items: Nested container.

    Yields:
        Each leaf item.
    """
    if isinstance(items, dict):
        items = items.values()
    for item in items:
        if isinstance(item, (list, tuple, dict)):
            yield from _flatten(item)
        else:
            yield item


def _open_text(path: Path) -> IO[str]:
    """Open a text file for reading, gzip-compressed or not.

    Args:
        path: File path; a '.gz' suffix selects gzip.

    Returns:
        Open text handle.
    """
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def _read_remark_scores(path: PathLike) -> dict[str, float]:
    """Return the HADDOCK score and energy terms from the REMARK header of a model PDB.

    Args:
        path: Model PDB, gzip-compressed or not.

    Returns:
        {term: value}, with the HADDOCK score as 'score'; empty if the header has none.
    """
    values: dict[str, float] = {}
    header: list[str] | None = None
    with _open_text(Path(path)) as handle:
        for line in handle:
            if not line.startswith("REMARK"):
                break
            text = line[len("REMARK"):].strip()
            if ":" not in text and "," in text:
                header = [name.strip() for name in text.split(",")]
                continue
            label, _, value = text.partition(":")
            label = label.strip()
            if label == "HADDOCK score":
                values["score"] = float(value)
            elif label == "energies" and header:
                values.update(zip(header, (float(item) for item in value.split(","))))
            elif label in REMARK_TERMS:
                values[REMARK_TERMS[label]] = float(value)
    return values


def _caprieval_reference(folder: PathLike) -> str:
    """Return what a caprieval folder compared models with.

    A workflow step records the reference in params.cfg; haddock3-analyse copies it into its folder as reference.pdb.

    Args:
        folder: Caprieval step or analysis folder.

    Returns:
        The reference file name, or 'best model'.
    """
    params = Path(folder) / "params.cfg"
    if params.is_file():
        match = re.search(r'^reference_fname\s*=\s*"([^"]*)"', params.read_text(), re.M)
        return Path(match.group(1)).name if match and match.group(1) else "best model"
    return "reference.pdb" if (Path(folder) / "reference.pdb").is_file() else "best model"


def _read_model(path: PathLike) -> Structure:
    """Parse a PDB file, gzip-compressed or not.

    Args:
        path: Model file.

    Returns:
        Biopython structure named after the file stem.
    """
    path = Path(path)
    with _open_text(path) as handle:
        text = handle.read()
    return PDBParser(QUIET=True).get_structure(path.name.split(".")[0], io.StringIO(text))


def _split_step_name(name: str) -> tuple[str, str]:
    """Split a step folder name such as '3_emref'.

    Args:
        name: Step folder name.

    Returns:
        (step number, module name), both as strings.

    Raises:
        ValueError: If the name is not a step folder name.
    """
    match = STEP_PATTERN.match(name)
    if match is None:
        raise ValueError(f"Not a haddock3 step folder: {name}")
    return match.group(1), match.group(2)

###################################### HADDOCK SETUP


class HaddockDocking:
    """Configure, write, run, and read back a HADDOCK3 docking of chains from one or more SyncedStructures.

    Partners are rigid bodies. Restraints are given on the original chains, in reference, sequence, or author
    numbering, and are written in HADDOCK3 numbering through the residue map. The 'alascan' method docks nothing: it
    refines the partners as they sit, written as one complex, and scans their interface.

    Attributes:
        method: Workflow preset, 'ab-initio', 'oracle', 'targeted', or 'alascan'.
        seed: Random seed of every CNS module and of the partner placement.
        ncores: Number of cores haddock3 may use.
        workdir: Directory holding inputs, restraints, config, logs, and the run.
        clearance_a: Smallest gap in A between separated partners.
        module_params: Per-module parameter overrides.
        global_params: Top-level haddock3 parameters.
        caprieval: Whether the workflow ends with caprieval.
        haddock3_exe: Path to the haddock3 executable, or None to search for it.
        partners: Normalized partners with their structure, chain keys, and segment ID.
        separated: Whether partners were moved apart before writing.
        placements: Rotation and translation of each moved partner.
    """

    method: str
    seed: int
    ncores: int
    workdir: Path
    clearance_a: float
    module_params: dict[str, dict[str, Any]]
    global_params: dict[str, Any]
    caprieval: bool
    haddock3_exe: PathLike | None
    partners: list[Partner]
    separated: bool
    placements: dict[int, Placement]
    _reference: PathLike | None
    _reference_source: ReferenceSource | None
    _aligned_pairs: set[tuple[str, str, ChainTarget]]
    _restraints: list[Restraint]
    _patches: list[Record]
    _chains: dict[ChainTarget, list[ResidueEntry]]
    _top_clusters: int | None

    def __init__(self, partners: Iterable[PartnerSpec], workdir: PathLike, method: str = "ab-initio",
                 seed: int = HADDOCK_DEFAULT_SEED, ncores: int | None = None, separate: str | bool = "auto",
                 clearance_a: float = HADDOCK_CLEARANCE_A, module_params: Mapping[str, Mapping[str, Any]] | None = None,
                 global_params: Mapping[str, Any] | None = None,
                 reference: PathLike | SyncedStructure | None = "auto", caprieval: bool = True,
                 haddock3_exe: PathLike | None = None) -> None:
        """Assign segment IDs and HADDOCK3 residue numbers, and place partners apart that could start docked.

        Partners may come from different structures. Their restraints, contacts, and reference complex are then
        matched through reference positions; see set_reference().

        Args:
            partners: (SyncedStructure, chain or list of chains) pairs, one per rigid body.
            workdir: Directory for inputs, restraints, config, log, and the run directory 'run'.
            method: 'ab-initio', 'oracle', 'targeted', or 'alascan' (no docking: refine, cluster, and alanine-scan
                the complex as given; see alascan_tables()).
            seed: Random seed, used as iniseed in every CNS module and for the placement of partners.
            ncores: Number of cores haddock3 may use, or None for every core this process may use (see
                core_count(): the affinity mask, cgroup quota, or batch-scheduler allocation, not the whole node).
            separate: 'auto' moves partners apart when two share a structure ID or touch (heavy atoms within
                clearance_a); True or False forces it. 'alascan' never moves partners.
            clearance_a: Smallest gap in A between the bounding spheres of separated partners.
            module_params: {module: {parameter: value}} overriding the preset of each module.
            global_params: Top-level haddock3 parameters, such as mode or max_cpus.
            reference: Reference complex for caprieval: 'auto' (the input complex when all partners share a
                structure ID), a SyncedStructure (chains matched by key; use set_reference() for a chain map), a PDB
                path already in partner segment IDs and HADDOCK3 numbering, or None.
            caprieval: Whether to end with caprieval, against the reference or, without one, the best-scoring model.
            haddock3_exe: Path to the haddock3 executable, or None to search the PATH and this Python's bin directory.

        Raises:
            ValueError: If method or separate is invalid, or separate is True for 'alascan'.
        """
        if method not in HADDOCK_METHODS:
            raise ValueError(f"method must be one of {HADDOCK_METHODS}.")
        if separate not in ("auto", True, False):
            raise ValueError("separate must be 'auto', True, or False.")
        if method == "alascan" and separate is True:
            raise ValueError("'alascan' scans the complex as given; partners cannot be separated.")
        self.method = method
        self.seed = int(seed)
        self.ncores = core_count() if ncores is None else int(ncores)
        self.workdir = Path(workdir).expanduser().resolve()
        self.clearance_a = float(clearance_a)
        self.module_params = {str(name): dict(params) for name, params in (module_params or {}).items()}
        self.global_params = dict(global_params or {})
        self.caprieval = bool(caprieval)
        self.haddock3_exe = haddock3_exe
        self._reference = "auto"
        self._reference_source = None
        self._aligned_pairs = set()
        self._restraints = []
        self._patches = []
        self._top_clusters = None
        self.partners = self._normalize_partners(partners)
        self._chains = self._assign_numbering()
        if method == "alascan":
            self.separated = False
            self._check_interface()
        else:
            self.separated = self._partners_overlap() if separate == "auto" else bool(separate)
        self.placements = self._place_partners() if self.separated else {}
        self._check_overrides()
        self.set_reference(reference)

###################################### PARTNERS AND NUMBERING

    def _normalize_partners(self, partners: Iterable[PartnerSpec]) -> list[Partner]:
        """Return one record per partner with its structure, chain keys, and segment ID.

        Args:
            partners: (SyncedStructure, chain or list of chains) pairs.

        Returns:
            Partner records in input order.

        Raises:
            ValueError: If there are too few or too many partners, or a chain is empty or used twice.
            TypeError: If a partner is not a (structure, chains) pair.
        """
        partners = list(partners)
        if len(partners) < 2:
            raise ValueError("Docking needs at least two partners.")
        if len(partners) > len(HADDOCK_SEGIDS):
            raise ValueError(f"At most {len(HADDOCK_SEGIDS)} partners are supported.")
        normalized: list[Partner] = []
        used: set[tuple[int, str]] = set()
        for index, item in enumerate(partners):
            if not (isinstance(item, (tuple, list)) and len(item) == 2 and hasattr(item[0], "chains")):
                raise TypeError("Each partner must be a (SyncedStructure, chain or list of chains) pair.")
            structure, chains = item
            names = [chains] if isinstance(chains, str) else list(chains)
            keys = list(dict.fromkeys(structure._chain_key(chain) for chain in names))
            if not keys:
                raise ValueError(f"Partner {index + 1} has no chains.")
            for key in keys:
                if (id(structure), key) in used:
                    raise ValueError(f"Chain {key} of {structure.structure_id} is used by more than one partner.")
                if not structure.chains[key].residues:
                    raise ValueError(f"Chain {key} of {structure.structure_id} has no polymer residues to dock.")
                used.add((id(structure), key))
            normalized.append({"structure": structure, "keys": keys, "segid": HADDOCK_SEGIDS[index]})
        return normalized

    def _assign_numbering(self) -> dict[ChainTarget, list[ResidueEntry]]:
        """Return the HADDOCK3 residue number of each polymer residue, per partner chain.

        Author numbers are kept when they are unique across the partner, have no insertion code, and fit in
        1..9999; otherwise the partner is renumbered from 1, leaving a gap between chains.

        Returns:
            {(partner index, chain key): residue entries in chain order}.

        Raises:
            ValueError: If a renumbered partner exceeds the PDB residue number range.
        """
        chains: dict[ChainTarget, list[ResidueEntry]] = {}
        for index, partner in enumerate(self.partners):
            records = [(key, partner["structure"].chains[key]) for key in partner["keys"]]
            ids = [residue.id for _, record in records for residue in record.residues]
            numbers = [number for _, number, _ in ids]
            keep = (all(icode == " " for _, _, icode in ids) and len(set(numbers)) == len(numbers)
                    and all(1 <= number <= HADDOCK_MAX_RESID for number in numbers))
            if not keep:
                log_info(f"HADDOCK3: partner {partner['segid']} renumbered from 1; see haddock_residue_map.tsv.")
            next_resid = 1
            for key, record in records:
                entries: list[ResidueEntry] = []
                for seq_index, residue in enumerate(record.residues):
                    resid = residue.id[1] if keep else next_resid + seq_index
                    entries.append({"resid": resid, "residue": residue, "sequence_index": seq_index,
                                    "reference_index": record.sequence_index_to_ref_index.get(seq_index)})
                next_resid += len(record.residues) + HADDOCK_CHAIN_GAP
                chains[(index, key)] = entries
            if not keep and next_resid - HADDOCK_CHAIN_GAP - 1 > HADDOCK_MAX_RESID:
                raise ValueError(f"Partner {partner['segid']} has too many residues for PDB numbering.")
        return chains

    def _partner_coords(self, index: int) -> np.ndarray:
        """Return the heavy-atom coordinates of a partner's polymer residues, before placement.

        Args:
            index: Partner index.

        Returns:
            (n_atoms, 3) array in A.
        """
        return np.array([atom.coord for key in self.partners[index]["keys"]
                         for entry in self._chains[(index, key)]
                         for atom in entry["residue"].get_atoms() if is_heavy_atom(atom)], dtype=np.float64)

    def _place_coords(self, index: int, coords: np.ndarray) -> np.ndarray:
        """Return coordinates moved by a partner's placement, or unchanged when it has none.

        Args:
            index: Partner index.
            coords: (n, 3) coordinates in A.

        Returns:
            Moved coordinates.
        """
        placement = self.placements.get(index)
        if placement is None:
            return coords
        return (coords - placement["origin"]) @ placement["rotation"].T + placement["center"]

    def _partners_overlap(self) -> bool:
        """Return whether two partners share a structure ID or have heavy atoms within clearance_a.

        Returns:
            Whether the partners should be moved apart.
        """
        structure_ids = [partner["structure"].structure_id for partner in self.partners]
        if len(set(structure_ids)) < len(structure_ids):
            return True
        coords = [self._partner_coords(index) for index in range(len(self.partners))]
        for i, j in itertools.combinations(range(len(coords)), 2):
            if ChainMixin._coords_touch(coords[i], coords[j], self.clearance_a):
                log_info(f"HADDOCK3: partners {self.partners[i]['segid']} and {self.partners[j]['segid']} come from "
                         f"different structures but touch; moving them apart.")
                return True
        return False

    def _check_interface(self) -> None:
        """Log a warning when no two partners touch, as 'alascan' then finds no interface to scan."""
        coords = [self._partner_coords(index) for index in range(len(self.partners))]
        if not any(ChainMixin._coords_touch(coords[i], coords[j], self.clearance_a)
                   for i, j in itertools.combinations(range(len(coords)), 2)):
            warnings.warn(f"No two partners have heavy atoms within {self.clearance_a} A; alascan will find no "
                          "interface to scan.", stacklevel=3)

    def _place_partners(self) -> dict[int, Placement]:
        """Place every partner after the first apart from the others, keeping the first in place.

        Each one is rotated randomly about its centroid and moved along a random direction until its bounding
        sphere clears all placed spheres by clearance_a; the result is then checked atom by atom.

        Returns:
            {partner index: rotation, origin, and new center}.

        Raises:
            RuntimeError: If two partners still touch after placement.
        """
        rng = np.random.default_rng(self.seed)
        coords = [self._partner_coords(index) for index in range(len(self.partners))]
        centers = [xyz.mean(axis=0) for xyz in coords]
        radii = [float(np.linalg.norm(xyz - center, axis=1).max()) for xyz, center in zip(coords, centers)]
        placed = [(centers[0], radii[0])]
        placements: dict[int, Placement] = {}
        for index in range(1, len(self.partners)):
            rotation = _random_rotation(rng)
            distance = radii[0] + radii[index] + self.clearance_a
            while True:
                for _ in range(HADDOCK_PLACEMENT_TRIALS):
                    center = centers[0] + _random_direction(rng) * distance
                    if all(np.linalg.norm(center - other) >= radius + radii[index] + self.clearance_a
                           for other, radius in placed):
                        break
                else:
                    distance *= HADDOCK_PLACEMENT_GROWTH
                    continue
                break
            placements[index] = {"rotation": rotation, "origin": centers[index], "center": center}
            placed.append((center, radii[index]))
        self.placements = placements
        moved = [self._place_coords(index, xyz) for index, xyz in enumerate(coords)]
        for i, j in itertools.combinations(range(len(moved)), 2):
            if ChainMixin._coords_touch(moved[i], moved[j], self.clearance_a):
                raise RuntimeError(f"Partners {i + 1} and {j + 1} are still within {self.clearance_a} A "
                                   "after placement.")
        log_info(f"HADDOCK3: moved {len(placements)} partner(s) apart with random rotations (seed {self.seed}).")
        return placements

    def residue_table(self) -> pd.DataFrame:
        """Return one row per written residue: partner, segment ID, HADDOCK3 number, and the original numbering.

        Returns:
            One row per residue written to the partner PDB files.
        """
        rows = []
        for (index, key), entries in self._chains.items():
            partner = self.partners[index]
            record = partner["structure"].chains[key]
            for entry in entries:
                residue = entry["residue"]
                reference_index = entry["reference_index"]
                rows.append({
                    "partner": index + 1, "segid": partner["segid"], "haddock_resid": entry["resid"],
                    "structure_id": partner["structure"].structure_id, "chain_key": key,
                    "chain_id": record.chain_id, "auth_chain_id": record.auth_chain_id,
                    "pdb_resnum": residue.id[1], "icode": residue.id[2].strip(),
                    "resname": residue.get_resname().strip().upper(),
                    "sequence_index": entry["sequence_index"] + 1,
                    "reference_index": None if reference_index is None else reference_index + 1,
                })
        frame = pd.DataFrame(rows, columns=HADDOCK_RESIDUE_MAP_COLUMNS)
        return frame.astype({"reference_index": "Int64"})

###################################### RESTRAINT RESOLUTION

    def _locate_chain(self, chain: ChainSpec) -> ChainTarget:
        """Find the partner holding a chain.

        A chain key matches before a label or author ID, so 'A' means the chain keyed 'A' even when another
        partner holds a chain whose label ID is 'A'.

        Args:
            chain: Chain key, label ID, or author ID unique among the partners, or (partner number from 1, chain).

        Returns:
            (partner index, chain key).

        Raises:
            ValueError: If the partner number is invalid, the chain is not in that partner, or the name is ambiguous.
            KeyError: If no partner holds the chain.
        """
        if isinstance(chain, (tuple, list)):
            number, name = chain
            index = int(number) - 1
            if not 0 <= index < len(self.partners):
                raise ValueError(f"No partner number {number}.")
            key = self.partners[index]["structure"]._chain_key(name)
            if key not in self.partners[index]["keys"]:
                raise ValueError(f"Chain {name} is not part of partner {number}.")
            return index, key
        matches = [(index, str(chain)) for index, partner in enumerate(self.partners) if str(chain) in partner["keys"]]
        if not matches:
            for index, partner in enumerate(self.partners):
                try:
                    key = partner["structure"]._chain_key(chain)
                except KeyError:
                    continue
                if key in partner["keys"]:
                    matches.append((index, key))
        if not matches:
            raise KeyError(f"Chain {chain} is not part of any docking partner.")
        if len(matches) > 1:
            raise ValueError(f"Chain {chain} matches several partners; pass (partner number, chain) instead.")
        return matches[0]

    def _haddock_numbers(self, chain: ChainSpec, positions: Iterable[int], numbering: str) -> tuple[int, list[int]]:
        """Translate 1-based positions on a chain into HADDOCK3 residue numbers.

        Positions that match no observed residue are logged and dropped.

        Args:
            chain: Chain identifier, see _locate_chain.
            positions: 1-based positions.
            numbering: 'reference' (reference sequence), 'sequence' (observed sequence), or 'author' (PDB number).

        Returns:
            (partner index, sorted HADDOCK3 residue numbers).

        Raises:
            ValueError: If numbering is invalid or no position matches an observed residue.
        """
        if numbering not in HADDOCK_NUMBERINGS:
            raise ValueError(f"numbering must be one of {HADDOCK_NUMBERINGS}.")
        index, key = self._locate_chain(chain)
        lookup: dict[int, list[int]] = {}
        for entry in self._chains[(index, key)]:
            if numbering == "author":
                position: int | None = entry["residue"].id[1]
            elif numbering == "sequence":
                position = entry["sequence_index"] + 1
            else:
                reference_index = entry["reference_index"]
                position = None if reference_index is None else reference_index + 1
            if position is not None:
                lookup.setdefault(position, []).append(entry["resid"])
        wanted = [int(position) for position in positions]
        resids: set[int] = set()
        missing: list[int] = []
        for position in wanted:
            if position in lookup:
                resids.update(lookup[position])
            else:
                missing.append(position)
        if missing:
            log_info(f"HADDOCK3: {len(missing)} {numbering} position(s) of chain {key} have no observed residue "
                     f"and were dropped: {missing[:20]}{' ...' if len(missing) > 20 else ''}")
        if not resids and wanted:
            raise ValueError(f"None of the {numbering} positions given for chain {key} match an observed residue.")
        return index, sorted(resids)

    def _check_partners_differ(self, index_a: int, index_b: int) -> None:
        """Check that two restrained selections belong to different partners.

        Args:
            index_a: Partner index of the first selection.
            index_b: Partner index of the second selection.

        Raises:
            ValueError: If both selections are in the same partner.
        """
        if index_a == index_b:
            raise ValueError("Restraints must join two different partners; both selections are in partner "
                             f"{self.partners[index_a]['segid']}.")

###################################### RESTRAINTS

    def add_ambiguous_restraints(self, chain_a: ChainSpec, active_a: Iterable[int], chain_b: ChainSpec,
                                 active_b: Iterable[int], passive_a: Iterable[int] = (),
                                 passive_b: Iterable[int] = (), numbering: str = "reference",
                                 distance: Sequence[float] = HADDOCK_RESTRAINT_DISTANCE,
                                 origin: str = "custom") -> Self:
        """Add ambiguous interaction restraints between active and passive residues of two chains.

        Each active residue of one side must lie near any active or passive residue of the other side.

        Args:
            chain_a: First chain, see _locate_chain.
            active_a: 1-based positions of the active residues of chain_a.
            chain_b: Second chain, in another partner.
            active_b: 1-based positions of the active residues of chain_b.
            passive_a: 1-based positions of passive residues of chain_a.
            passive_b: 1-based positions of passive residues of chain_b.
            numbering: 'reference', 'sequence', or 'author'.
            distance: (distance, lower correction, upper correction) in A.
            origin: Label stored with the restraint in restraint_table.

        Returns:
            self.

        Raises:
            ValueError: If neither side has active residues.
        """
        act_in_a, act_in_b, pas_in_a, pas_in_b = (list(values) for values in (active_a, active_b, passive_a, passive_b))
        index_a, act_a = self._haddock_numbers(chain_a, act_in_a, numbering) if act_in_a else (self._locate_chain(chain_a)[0], [])
        index_b, act_b = self._haddock_numbers(chain_b, act_in_b, numbering) if act_in_b else (self._locate_chain(chain_b)[0], [])
        pas_a = self._haddock_numbers(chain_a, pas_in_a, numbering)[1] if pas_in_a else []
        pas_b = self._haddock_numbers(chain_b, pas_in_b, numbering)[1] if pas_in_b else []
        self._check_partners_differ(index_a, index_b)
        if not act_a and not act_b:
            raise ValueError("Ambiguous restraints need active residues on at least one side.")
        self._restraints.append({
            "kind": "ambiguous", "origin": origin, "partner_a": index_a, "residues_a": act_a, "passive_a": pas_a,
            "partner_b": index_b, "residues_b": act_b, "passive_b": pas_b, "distance": tuple(distance),
        })
        return self

    def add_unambiguous_restraint(self, chain_a: ChainSpec, positions_a: Iterable[int], chain_b: ChainSpec,
                                  positions_b: Iterable[int], numbering: str = "reference",
                                  distance: Sequence[float] = HADDOCK_RESTRAINT_DISTANCE,
                                  origin: str = "custom") -> Self:
        """Add one unambiguous restraint joining a patch of chain_a to a patch of chain_b.

        Args:
            chain_a: First chain, see _locate_chain.
            positions_a: 1-based positions of the residues of the first patch.
            chain_b: Second chain, in another partner.
            positions_b: 1-based positions of the residues of the second patch.
            numbering: 'reference', 'sequence', or 'author'.
            distance: (distance, lower correction, upper correction) in A.
            origin: Label stored with the restraint in restraint_table.

        Returns:
            self.

        Raises:
            ValueError: If either patch has no residues.
        """
        index_a, resids_a = self._haddock_numbers(chain_a, positions_a, numbering)
        index_b, resids_b = self._haddock_numbers(chain_b, positions_b, numbering)
        self._check_partners_differ(index_a, index_b)
        if not resids_a or not resids_b:
            raise ValueError("Unambiguous restraints need residues on both sides.")
        self._restraints.append({
            "kind": "unambiguous", "origin": origin, "partner_a": index_a, "residues_a": resids_a, "passive_a": [],
            "partner_b": index_b, "residues_b": resids_b, "passive_b": [], "distance": tuple(distance),
        })
        return self

    def add_patch_pair_restraints(self, pairs: Iterable[PatchPair], kind: str = "ambiguous",
                                  numbering: str = "reference",
                                  distance: Sequence[float] = HADDOCK_RESTRAINT_DISTANCE) -> Self:
        """Add restraints joining each couple of patches.

        With 'ambiguous', every residue of each patch is active; with 'unambiguous', each couple is one restraint.

        Args:
            pairs: ((chain_a, positions_a), (chain_b, positions_b)) couples.
            kind: 'ambiguous' or 'unambiguous'.
            numbering: 'reference', 'sequence', or 'author'.
            distance: (distance, lower correction, upper correction) in A.

        Returns:
            self.

        Raises:
            ValueError: If kind is invalid or no pairs are given.
        """
        if kind not in HADDOCK_RESTRAINT_KINDS:
            raise ValueError(f"kind must be one of {HADDOCK_RESTRAINT_KINDS}.")
        add = self.add_ambiguous_restraints if kind == "ambiguous" else self.add_unambiguous_restraint
        count = 0
        for (chain_a, positions_a), (chain_b, positions_b) in pairs:
            add(chain_a, list(positions_a), chain_b, list(positions_b), numbering=numbering, distance=distance,
                origin="patch_pair")
            count += 1
        if not count:
            raise ValueError("No patch pairs given.")
        return self

    def _contact_source(self, source: SyncedStructure | None,
                        chain_map: Mapping[ChainSpec, str] | None) -> ReferenceSource:
        """Resolve the structure to read contacts from and how its chains map onto partner chains.

        Args:
            source: Structure to read contacts from; None for the reference set with set_reference(), or else the
                partners' own structure.
            chain_map: {partner chain: source chain}, or None to match chain keys.

        Returns:
            (source structure, {source chain key: (partner index, chain key)}).

        Raises:
            ValueError: If the default source is ambiguous, chain keys match ambiguously, or the source holds fewer
                than two partners.
        """
        if source is None and self._reference_source is not None and chain_map is None:
            return self._reference_source
        if source is None:
            structures = {id(partner["structure"]) for partner in self.partners}
            if len(structures) > 1:
                raise ValueError("Partners come from different structures; call set_reference() with the complex, "
                                 "or pass source= (and chain_map= if its chains differ).")
            source = self.partners[0]["structure"]
        if chain_map is None:
            matched = [(key, (index, key)) for index, partner in enumerate(self.partners)
                       for key in partner["keys"] if key in source]
            if len({key for key, _ in matched}) < len(matched):
                raise ValueError("Several partners have a chain with the same key; pass chain_map= with "
                                 "(partner number, chain) keys.")
            resolved = dict(matched)
        else:
            resolved = {source._chain_key(source_chain): self._locate_chain(chain)
                        for chain, source_chain in chain_map.items()}
        if len({index for index, _ in resolved.values()}) < 2:
            raise ValueError("The source structure must hold chains of at least two partners.")
        return source, resolved

    def _position_map(self, source: SyncedStructure, source_key: str, target: ChainTarget) -> dict[int, int]:
        """Map the residues of a source chain onto a partner chain.

        A chain of a partner's own structure maps onto itself. Any other chain maps through reference positions; when
        the two reference sequences differ (tags, truncations, mutations), through a global alignment of them.

        Args:
            source: Structure holding the source chain.
            source_key: Chain key in the source.
            target: (partner index, chain key) to map onto.

        Returns:
            {0-based source sequence index: 0-based target sequence index}, for residues with a counterpart.
        """
        index, key = target
        target_structure = self.partners[index]["structure"]
        source_record, target_record = source.chains[source_key], target_structure.chains[key]
        if source is target_structure and source_key == key:
            return {seq: seq for seq in range(len(source_record.residues))}
        ref_to_ref = reference_index_map(source_record.reference_sequence, target_record.reference_sequence)
        if source_record.reference_sequence != target_record.reference_sequence:
            source_ref, target_ref = source_record.reference_sequence, target_record.reference_sequence
            same = sum(source_ref[a] == target_ref[b] for a, b in ref_to_ref.items())
            pair = (str(source.structure_id), source_key, target)
            if pair not in self._aligned_pairs:
                self._aligned_pairs.add(pair)
                log_info(f"HADDOCK3: reference sequences of {source.structure_id} chain {source_key} and "
                         f"{target_structure.structure_id} chain {key} differ; matched by alignment "
                         f"({same} identical of {len(target_ref)} positions).")
        return sequence_index_map(source_record, target_record, ref_to_ref)

    def partner_contact_map(self, cutoff: float = CONTACTS_CUTOFF_NM, scheme: str = "ca",
                              source: SyncedStructure | None = None, chain_map: Mapping[ChainSpec, str] | None = None,
                              ) -> dict[tuple[ChainTarget, ChainTarget], list[tuple[int, int]]]:
        """Return the contacts between partners, per chain pair.

        Contacts without a counterpart in a partner chain are dropped.

        Args:
            cutoff: Contact cutoff in nm.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            source: Structure to read contacts from; None for the reference set with set_reference(), or else the
                partners' own structure.
            chain_map: {partner chain: source chain} when the source chains differ.

        Returns:
            {((partner, chain), (partner, chain)): sorted (sequence_index_a, sequence_index_b) pairs}, 0-based.
        """
        source, resolved = self._contact_source(source, chain_map)
        contacts = source.contact_table(cutoff=cutoff, contact_type="inter", chain_ids=list(resolved),
                                           scheme=scheme, moieties=POLYMER_MOIETIES)
        maps = {key: self._position_map(source, key, target) for key, target in resolved.items()}
        pairs: dict[tuple[ChainTarget, ChainTarget], list[tuple[int, int]]] = {}
        dropped = 0
        columns = (contacts[name].tolist() for name in ("chain_key_a", "sequence_index_a", "chain_key_b",
                                                         "sequence_index_b"))
        for key_1, index_1, key_2, index_2 in zip(*columns):
            seq_1, seq_2 = int(index_1) - 1, int(index_2) - 1
            target_1, target_2 = resolved[key_1], resolved[key_2]
            if target_1[0] == target_2[0]:
                continue
            mapped_1, mapped_2 = maps[key_1].get(seq_1), maps[key_2].get(seq_2)
            if mapped_1 is None or mapped_2 is None:
                dropped += 1
                continue
            if target_1 > target_2:
                target_1, target_2, mapped_1, mapped_2 = target_2, target_1, mapped_2, mapped_1
            pairs.setdefault((target_1, target_2), []).append((mapped_1, mapped_2))
        if dropped:
            log_info(f"HADDOCK3: {dropped} contact(s) have no counterpart in the partner chains and were dropped.")
        return {chains: sorted(set(values)) for chains, values in pairs.items()}

    def add_contact_restraints(self, kind: str = "ambiguous", cutoff: float = CONTACTS_CUTOFF_NM, scheme: str = "ca",
                               source: SyncedStructure | None = None, chain_map: Mapping[ChainSpec, str] | None = None,
                               patch_gap: float = HADDOCK_PATCH_GAP,
                               distance: Sequence[float] = HADDOCK_RESTRAINT_DISTANCE) -> Self:
        """Add restraints from inter-partner contacts.

        'ambiguous' makes the contacted residues of each chain active. 'unambiguous' clusters the contacts of each
        chain pair into patches, joining contacts whose positions differ by at most patch_gap on both sequences,
        and adds one restraint per patch pair. Positions are reference positions when every contacted residue of
        the chain has one, otherwise observed sequence positions.

        Args:
            kind: 'ambiguous' or 'unambiguous'.
            cutoff: Contact cutoff in nm.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            source: Structure to read contacts from; None for the reference set with set_reference(), or else the
                partners' own structure. Another structure's contacts are transferred through reference positions.
            chain_map: {partner chain: source chain} when the source chains differ.
            patch_gap: Largest position difference that joins two contacts into one patch.
            distance: (distance, lower correction, upper correction) in A.

        Returns:
            self.

        Raises:
            ValueError: If kind is invalid or the partners have no contacts.
        """
        if kind not in HADDOCK_RESTRAINT_KINDS:
            raise ValueError(f"kind must be one of {HADDOCK_RESTRAINT_KINDS}.")
        pairs = self.partner_contact_map(cutoff=cutoff, scheme=scheme, source=source, chain_map=chain_map)
        if not pairs:
            raise ValueError("No contacts between partners at this cutoff.")
        for (target_a, target_b), contacts in pairs.items():
            chain_a, chain_b = (target_a[0] + 1, target_a[1]), (target_b[0] + 1, target_b[1])
            if kind == "ambiguous":
                self.add_ambiguous_restraints(chain_a, sorted({a + 1 for a, _ in contacts}), chain_b,
                                              sorted({b + 1 for _, b in contacts}), numbering="sequence",
                                              distance=distance, origin="contacts")
                continue
            axes = [self._patch_positions(target, [pair[side] for pair in contacts])
                    for side, target in enumerate((target_a, target_b))]
            points = [(axes[0][a], axes[1][b]) for a, b in contacts]
            for members in cluster_contact_patches(points, patch_gap):
                side_a = sorted({contacts[m][0] + 1 for m in members})
                side_b = sorted({contacts[m][1] + 1 for m in members})
                self.add_unambiguous_restraint(chain_a, side_a, chain_b, side_b, numbering="sequence",
                                               distance=distance, origin="contact_patch")
                self._patches.append({
                    "patch": len(self._patches) + 1,
                    "segid_a": self.partners[target_a[0]]["segid"], "chain_a": target_a[1],
                    "residues_a": self._restraints[-1]["residues_a"],
                    "segid_b": self.partners[target_b[0]]["segid"], "chain_b": target_b[1],
                    "residues_b": self._restraints[-1]["residues_b"], "n_pairs": len(members),
                })
        log_info(f"HADDOCK3: added {kind} restraints from contacts between {len(pairs)} chain pair(s).")
        return self

    def _patch_positions(self, target: ChainTarget, seq_indices: Sequence[int]) -> dict[int, int]:
        """Return the clustering position of contacted residues: reference positions when all are mapped.

        Args:
            target: (partner index, chain key).
            seq_indices: 0-based sequence indices of contacted residues.

        Returns:
            {0-based sequence index: clustering position}.
        """
        mapping = self.partners[target[0]]["structure"].chains[target[1]].sequence_index_to_ref_index
        if all(index in mapping for index in seq_indices):
            return {index: mapping[index] for index in seq_indices}
        return {index: index for index in seq_indices}

    def add_scan_residues(self, chain: ChainSpec, positions: Iterable[int], numbering: str = "reference") -> Self:
        """Restrict the alanine scan of 'alascan' to some residues; without any, the whole interface is scanned.

        alascan mutates only residues within its int_cutoff of the other partner in each model. haddock3 2026.7 fails
        when a model ends up with no residue to mutate, so pick residues that stay at the interface.

        Args:
            chain: Chain, see _locate_chain.
            positions: 1-based positions of the residues to scan.
            numbering: 'reference', 'sequence', or 'author'.

        Returns:
            self.

        Raises:
            ValueError: If the method is not 'alascan'.
        """
        if self.method != "alascan":
            raise ValueError("Scan residues need method='alascan'.")
        index, resids = self._haddock_numbers(chain, positions, numbering)
        params = self.module_params.setdefault("alascan", {})
        name = f"resdic_{self.partners[index]['segid']}"
        params[name] = sorted(set(params.get(name, [])) | set(resids))
        return self

    def clear_restraints(self) -> Self:
        """Remove every restraint and patch.

        Returns:
            self.
        """
        self._restraints.clear()
        self._patches.clear()
        return self

    def restraint_table(self) -> pd.DataFrame:
        """Return one row per restraint set with segment IDs, HADDOCK3 residue numbers, and distance.

        Returns:
            One row per restraint set.
        """
        rows = [{
            "kind": item["kind"], "origin": item["origin"],
            "segid_a": self.partners[item["partner_a"]]["segid"], "residues_a": item["residues_a"],
            "passive_a": item["passive_a"],
            "segid_b": self.partners[item["partner_b"]]["segid"], "residues_b": item["residues_b"],
            "passive_b": item["passive_b"], "distance": item["distance"],
        } for item in self._restraints]
        return pd.DataFrame(rows, columns=HADDOCK_RESTRAINT_COLUMNS)

    def contact_patch_table(self) -> pd.DataFrame:
        """Return contact patches behind unambiguous contact restraints, in HADDOCK3 residue numbers.

        Returns:
            One row per contact patch.
        """
        return pd.DataFrame(self._patches, columns=HADDOCK_PATCH_COLUMNS)

    def restraint_text(self, kind: str) -> str:
        """Return the CNS restraint file text of one kind.

        Args:
            kind: 'ambiguous' or 'unambiguous'.

        Returns:
            The .tbl text, or an empty string when there are no such restraints.

        Raises:
            ValueError: If kind is invalid.
        """
        if kind not in HADDOCK_RESTRAINT_KINDS:
            raise ValueError(f"kind must be one of {HADDOCK_RESTRAINT_KINDS}.")
        blocks = []
        for item in self._restraints:
            if item["kind"] != kind:
                continue
            segid_a = self.partners[item["partner_a"]]["segid"]
            segid_b = self.partners[item["partner_b"]]["segid"]
            bounds = " ".join(f"{value:.1f}" for value in item["distance"])
            if kind == "unambiguous":
                blocks.append(f"assign {_tbl_selection(segid_a, item['residues_a'])}\n"
                              f"       {_tbl_selection(segid_b, item['residues_b'])} {bounds}\n")
                continue
            for active, own, targets, other in (
                    (item["residues_a"], segid_a, set(item["residues_b"]) | set(item["passive_b"]), segid_b),
                    (item["residues_b"], segid_b, set(item["residues_a"]) | set(item["passive_a"]), segid_a)):
                if not targets:
                    continue
                for resid in active:
                    choices = "\n        or\n".join(f"       (resi {target} and segid {other})"
                                                   for target in sorted(targets))
                    blocks.append(f"assign (resi {resid} and segid {own})\n(\n{choices}\n) {bounds}\n")
        return "\n".join(blocks)

###################################### CONFIGURATION

    @property
    def run_dir(self) -> Path:
        """Directory haddock3 runs in."""
        return self.workdir / "run"

    @property
    def input_paths(self) -> InputPaths:
        """Paths of the files export_inputs() creates."""
        inputs = self.workdir / "inputs"
        if self.method == "alascan":
            structure_ids = dict.fromkeys(self._safe_name(partner["structure"].structure_id) for partner in self.partners)
            segids = "".join(partner["segid"] for partner in self.partners)
            molecules = [inputs / f"complex_{'-'.join(structure_ids)}_{segids}.pdb"]
        else:
            molecules = [inputs / f"partner_{partner['segid']}_{self._safe_name(partner['structure'].structure_id)}"
                                  f"_{'-'.join(self._safe_name(key) for key in partner['keys'])}.pdb"
                         for partner in self.partners]
        return {
            "molecules": molecules,
            "reference": inputs / "reference.pdb",
            "ambiguous": self.workdir / "restraints" / "ambig.tbl",
            "unambiguous": self.workdir / "restraints" / "unambig.tbl",
            "haddock_residue_map": self.workdir / "haddock_residue_map.tsv",
            "config": self.workdir / "config.toml",
            "log": self.workdir / "haddock3.log",
        }

    def _restraint_path(self, kind: str) -> Path:
        """Return the restraint file path of one kind.

        Args:
            kind: 'ambiguous' or 'unambiguous'.

        Returns:
            Path of ambig.tbl or unambig.tbl.
        """
        paths = self.input_paths
        return paths["ambiguous"] if kind == "ambiguous" else paths["unambiguous"]

    @staticmethod
    def _safe_name(text: object) -> str:
        """Return text usable in a file name.

        Args:
            text: Any value; it is converted to a string.

        Returns:
            The text with characters outside letters, digits, dot, and dash replaced by underscores.
        """
        return re.sub(r"[^A-Za-z0-9.-]", "_", str(text))

    def _check_overrides(self) -> None:
        """Check module and global parameter overrides.

        Raises:
            ValueError: For overrides of unknown modules or protected globals, or an iniseed different from seed.
        """
        modules = {name for name, _ in HADDOCK_PRESETS[self.method]} | {"caprieval"}
        unknown = set(self.module_params) - modules
        if unknown:
            raise ValueError(f"module_params names modules not in the '{self.method}' workflow: {sorted(unknown)}.")
        for name, params in self.module_params.items():
            if "iniseed" in params and int(params["iniseed"]) != self.seed:
                raise ValueError(f"iniseed of {name} must equal seed ({self.seed}); set seed= instead.")
        protected = HADDOCK_PROTECTED_GLOBALS & set(self.global_params)
        if protected:
            raise ValueError(f"global_params cannot set {sorted(protected)}; they are managed by HaddockDocking.")

    def _ensure_method_restraints(self) -> None:
        """Add the restraints a method needs by default: contact patches for 'oracle', a check for 'targeted'.

        Raises:
            ValueError: If 'targeted' has no patch couples.
        """
        if self.method == "oracle" and not self._restraints:
            log_info("HADDOCK3: oracle without restraints; adding unambiguous contact-patch restraints.")
            self.add_contact_restraints(kind="unambiguous")
        if self.method == "targeted" and not any(item["origin"] == "patch_pair" for item in self._restraints):
            raise ValueError("The 'targeted' method needs patch couples; add them with add_patch_pair_restraints().")

    def set_reference(self, reference: PathLike | SyncedStructure | None,
                      chain_map: Mapping[ChainSpec, str] | None = None) -> Self:
        """Set the reference complex that caprieval compares models with.

        A SyncedStructure is written to inputs/reference.pdb in partner segment IDs and HADDOCK3 numbering, its
        residues matched to the partner chains through reference positions, so it can be another entry than the
        partners (for example the holo complex of two apo partners). It also becomes the default source of
        add_contact_restraints().

        Args:
            reference: 'auto' (the input complex when all partners share a structure ID), a SyncedStructure, a PDB
                path already in partner segment IDs and HADDOCK3 numbering, or None for no reference.
            chain_map: {partner chain: reference chain} for a SyncedStructure, or None to match chain keys. Partner
                chains are named as in restraints, such as 'A' or (2, 'A').

        Returns:
            self.
        """
        if isinstance(reference, _structure.SyncedStructure):
            structure = reference
            self._reference_source = self._contact_source(structure, chain_map)
            self._reference = self.input_paths["reference"]
            mapped = len(self._reference_rows())
            total = sum(len(entries) for entries in self._chains.values())
            log_info(f"HADDOCK3: reference {structure.structure_id} covers {mapped} of {total} partner residues.")
        else:
            self._reference_source = None
            self._reference = reference
        return self

    def _reference_rows(self) -> list[tuple[int, ResidueEntry, Residue]]:
        """Return the reference residue matched to each partner residue, for a reference set as a structure.

        Returns:
            (partner index, partner residue entry, reference residue) triples, in partner chain order.
        """
        if self._reference_source is None:
            return []
        source, resolved = self._reference_source
        order = list(self._chains)
        rows = []
        for source_key, target in sorted(resolved.items(), key=lambda item: order.index(item[1])):
            inverse = {seq: source_seq for source_seq, seq in self._position_map(source, source_key, target).items()}
            residues = source.chains[source_key].residues
            rows += [(target[0], entry, residues[inverse[entry["sequence_index"]]])
                     for entry in self._chains[target] if entry["sequence_index"] in inverse]
        return rows

    def _reference_path(self) -> Path | None:
        """Return the reference complex for caprieval and run_caprieval().

        Returns:
            The reference path, or None when none is available.

        Raises:
            FileNotFoundError: If a given reference file does not exist.
        """
        reference = self._reference
        if reference == "auto":
            structure_ids = {partner["structure"].structure_id for partner in self.partners}
            return self.input_paths["reference"] if len(structure_ids) == 1 else None
        if reference is None:
            return None
        path = Path(reference).expanduser().resolve()
        if self._reference_source is None and not path.is_file():
            raise FileNotFoundError(f"Reference structure not found: {path}")
        return path

    def _caprieval_params(self) -> dict[str, Any] | None:
        """Return caprieval parameters.

        Without a reference, caprieval compares every model with the best-scoring one, so scores are always
        tabulated.

        Returns:
            The parameters, or None when caprieval is off.
        """
        if not self.caprieval:
            return None
        reference = self._reference_path()
        params: dict[str, Any] = {"receptor_chain": self.partners[0]["segid"],
                                  "ligand_chains": [partner["segid"] for partner in self.partners[1:]]}
        if reference is None:
            log_info("HADDOCK3: no reference structure; caprieval compares models with the best-scoring model.")
        else:
            params = {"reference_fname": str(reference), **params}
        return params

    def workflow_modules(self) -> list[tuple[str, dict[str, Any]]]:
        """Return the workflow with seed, restraints, and overrides applied.

        seletopclusts gets its top_clusters from the clustfcc step of the run, once run_haddock() has counted them.

        Returns:
            (module, parameters) pairs in workflow order.
        """
        self._ensure_method_restraints()
        has = {kind: any(item["kind"] == kind for item in self._restraints) for kind in HADDOCK_RESTRAINT_KINDS}
        workflow: list[tuple[str, dict[str, Any]]] = []
        for name, preset in HADDOCK_PRESETS[self.method]:
            params = dict(preset)
            if name in HADDOCK_RESTRAINT_MODULES:
                if has["ambiguous"]:
                    params["ambig_fname"] = str(self.input_paths["ambiguous"])
                if has["unambiguous"]:
                    params["unambig_fname"] = str(self.input_paths["unambiguous"])
            if name == "rigidbody" and params.get("ranair") and has["ambiguous"]:
                log_info("HADDOCK3: ambiguous restraints given; rigidbody ranair is turned off to use them.")
                params["ranair"] = False
            params.update(self.module_params.get(name, {}))
            if name == "seletopclusts" and "top_clusters" not in params and self._top_clusters is not None:
                params["top_clusters"] = self._top_clusters
            if name in HADDOCK_SEEDED_MODULES:
                params["iniseed"] = self.seed
            workflow.append((name, params))
        capri = self._caprieval_params()
        if capri is not None:
            capri.update(self.module_params.get("caprieval", {}))
            workflow.append(("caprieval", capri))
        return workflow

    def config_toml(self, run_dir: PathLike | None = None, stop: int | None = None) -> str:
        """Return the haddock3 TOML configuration.

        Args:
            run_dir: Run directory to write in the config, or None for workdir/run.
            stop: Number of workflow modules to write, or None for all.

        Returns:
            The TOML text.
        """
        globals_ = {"run_dir": str(run_dir or self.run_dir),
                    "molecules": [str(path) for path in self.input_paths["molecules"]],
                    "ncores": self.ncores, "mode": "local", **self.global_params}
        lines = [f"# HADDOCK3 configuration written by synced_structure (method: {self.method}, seed: {self.seed})"]
        lines += [f"{name} = {_toml_literal(value)}" for name, value in globals_.items()]
        for name, params in self.workflow_modules()[:stop]:
            lines += ["", f"[{name}]"] + [f"{key} = {_toml_literal(value)}" for key, value in params.items()]
        return "\n".join(lines) + "\n"

###################################### WRITING

    def _partner_rows(self, partner_indices: Iterable[int]) -> list[tuple[int, ResidueEntry, Residue]]:
        """Return the residues of some partners as written to their own PDB files.

        Args:
            partner_indices: Partners to include.

        Returns:
            (partner index, residue entry, residue) triples, in partner chain order.
        """
        return [(index, entry, entry["residue"]) for index in partner_indices
                for key in self.partners[index]["keys"] for entry in self._chains[(index, key)]]

    def _write_pdb(self, path: Path, rows: Sequence[tuple[int, ResidueEntry, Residue]], placed: bool) -> Path:
        """Write the heavy atoms of residues to a PDB file, each under its partner's segment ID and HADDOCK3 number.

        Args:
            path: Output file.
            rows: (partner index, residue entry, residue to write) triples; a TER record ends each run of one chain.
            placed: Whether to apply each partner's placement.

        Returns:
            The path written.
        """
        lines = []
        serial = 1
        previous: tuple[int, int] | None = None
        for index, entry, residue in rows:
            chain = (index, id(entry["residue"].get_parent()))
            if previous is not None and chain != previous:
                lines.append("TER")
            previous = chain
            segid = self.partners[index]["segid"]
            resname = residue.get_resname().strip().upper()
            for atom in residue.get_atoms():
                if not is_heavy_atom(atom):
                    continue
                coord = atom.coord.astype(np.float64)
                if placed:
                    coord = self._place_coords(index, coord[None])[0]
                lines.append(_pdb_atom_line(serial, atom, resname, segid, entry["resid"], coord))
                serial += 1
        if previous is not None:
            lines.append("TER")
        lines.append("END")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n")
        return path

    @staticmethod
    def _validate_tbl(text: str, path: Path) -> None:
        """Check a restraint file with haddock3's own validator, when haddock3 is importable.

        Args:
            text: Restraint file text.
            path: File the text is written to, for the error message.

        Raises:
            ValueError: If haddock3 rejects the restraints.
        """
        try:
            from haddock.libs.librestraints import validate_tbldata
        except ImportError:
            return
        try:
            validate_tbldata(text)
        except Exception as error:
            raise ValueError(f"HADDOCK3 rejected the restraints in {path}: {error}") from error

    def export_inputs(self) -> Path:
        """Write partner PDBs (one complex PDB for 'alascan'), the reference, restraints, residue map, and config.toml.

        Returns:
            Path of config.toml.
        """
        start = time.perf_counter()
        paths = self.input_paths
        config = self.config_toml()
        if self.method == "alascan":
            self._write_pdb(paths["molecules"][0], self._partner_rows(range(len(self.partners))), placed=False)
        else:
            for index, path in enumerate(paths["molecules"]):
                self._write_pdb(path, self._partner_rows([index]), placed=True)
        if self._reference_source is not None:
            self._write_pdb(paths["reference"], self._reference_rows(), placed=False)
        elif self._reference == "auto" and self._reference_path() == paths["reference"]:
            self._write_pdb(paths["reference"], self._partner_rows(range(len(self.partners))), placed=False)
        for kind in HADDOCK_RESTRAINT_KINDS:
            text = self.restraint_text(kind)
            restraint_path = self._restraint_path(kind)
            if text:
                restraint_path.parent.mkdir(parents=True, exist_ok=True)
                restraint_path.write_text(text)
                self._validate_tbl(text, restraint_path)
            elif restraint_path.exists():
                restraint_path.unlink()
        self.residue_table().to_csv(paths["haddock_residue_map"], sep=EXPORT_SEP, index=False)
        paths["config"].write_text(config)
        log_completion(start, f"Wrote HADDOCK3 '{self.method}' set up in {self.workdir}")
        return paths["config"]

###################################### RUNNING

    def _find_executable(self, name: str = HADDOCK_EXECUTABLE) -> str:
        """Find a haddock3 command.

        Args:
            name: Command name; with haddock3_exe set, other commands are looked up next to it.

        Returns:
            Path of the executable.

        Raises:
            FileNotFoundError: If it cannot be found.
        """
        candidates: list[PathLike | None]
        if self.haddock3_exe:
            candidates = [self.haddock3_exe if name == HADDOCK_EXECUTABLE else Path(self.haddock3_exe).parent / name]
        else:
            candidates = [shutil.which(name), Path(sys.executable).parent / name]
        for candidate in candidates:
            if candidate and Path(candidate).is_file():
                return str(candidate)
        raise FileNotFoundError(f"{name} executable not found; install haddock3 with 'pip install haddock3' "
                                "or pass haddock3_exe=.")

    def _run_command(self, args: Sequence[str], log_path: Path, name: str = HADDOCK_EXECUTABLE,
                     progress: bool = False, append: bool = False) -> None:
        """Run a haddock3 command, writing its output to a log file.

        The command is stopped if Python is interrupted while it runs.

        Args:
            args: Command-line arguments.
            log_path: Log file for stdout and stderr.
            name: Command name.
            progress: Whether to log module starts, finishes, and model counts while the command runs.
            append: Whether to append to the log file instead of replacing it.

        Raises:
            RuntimeError: If the command exits with a nonzero code.
        """
        offset = log_path.stat().st_size if append and log_path.is_file() else 0
        with open(log_path, "a" if append else "w") as handle:
            process = subprocess.Popen([self._find_executable(name), *args], cwd=self.workdir, stdout=handle,
                                       stderr=subprocess.STDOUT)
            try:
                if progress:
                    self._watch_run(process, log_path, offset)
                process.wait()
            finally:
                if process.poll() is None:          # Python was interrupted: stop haddock3 too
                    process.terminate()
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
        if process.returncode != 0:
            raise RuntimeError(f"{name} exited with code {process.returncode}; see {log_path}\n{_log_tail(log_path)}")

    def _model_counts(self) -> list[int | None]:
        """Return how many models each step of the workflow writes, as far as the configuration tells.

        Returns:
            One count per module of workflow_modules(), 0 for steps that write no models, None when unknown. After
            seletopclusts the count is an upper bound: clusters can have fewer than top_models members.
        """
        counts: list[int | None] = []
        flowing: int | None = len(self.input_paths["molecules"])
        for name, params in self.workflow_modules():
            if name == "topoaa":
                counts.append(flowing)
            elif name in HADDOCK_DEFAULT_SAMPLING:
                flowing = int(params.get("sampling", HADDOCK_DEFAULT_SAMPLING[name]))
                counts.append(flowing)
            elif name in HADDOCK_REFINE_MODULES:
                flowing = None if flowing is None else flowing * int(params.get("sampling_factor", 1))
                counts.append(flowing)
            else:
                if name == "seletop" and flowing is not None:
                    flowing = min(flowing, int(params.get("select", HADDOCK_DEFAULT_SELECT)))
                elif name == "seletopclusts":
                    top = params.get("top_clusters")
                    flowing = None if top is None else int(top) * int(params.get("top_models", HADDOCK_DEFAULT_TOP_MODELS))
                elif name not in HADDOCK_PASSTHROUGH_MODULES:
                    flowing = None
                counts.append(0)
        return counts

    def _step_model_count(self, step: int, params: Mapping[str, Any]) -> int | None:
        """Return how many models a refinement step writes, from the models the previous step passed on.

        Args:
            step: Number of the previous step.
            params: Parameters of the refinement module.

        Returns:
            The count, or None when the previous step's output cannot be read.
        """
        try:
            return len(self.score_table(step)) * int(params.get("sampling_factor", 1))
        except (FileNotFoundError, KeyError, ValueError):
            return None

    def _watch_run(self, process: subprocess.Popen[bytes], log_path: Path, offset: int = 0) -> None:
        """Log progress of a running haddock3 workflow until it exits.

        Module starts and finishes come from the haddock3 log; model counts from the PDB files in the step folder.

        Args:
            process: The running haddock3 process.
            log_path: Log file the process writes to.
            offset: Position in the log file where this process's output starts.
        """
        workflow = self.workflow_modules()
        modules, params = [name for name, _ in workflow], [values for _, values in workflow]
        expected = self._model_counts()
        width = len(str(len(modules)))
        started = re.compile(r"Running \[(\w+)\] module")
        finished = re.compile(r"Module \[(\w+)\] finished")
        run_start = step_start = time.perf_counter()
        step, last_count, last_report = -1, -1, run_start

        def label(index: int) -> str:
            """Return the log prefix of a step, such as 'HADDOCK3 [2/6] rigidbody'.

            Args:
                index: Step number.

            Returns:
                The prefix.
            """
            return f"HADDOCK3 [{index + 1:>{width}}/{len(modules)}] {modules[index]}"

        def report() -> None:
            """Log the models written by the current step, when the count changed since the last line."""
            nonlocal last_count, last_report
            last_report = time.perf_counter()
            if step < 0 or not expected[step]:
                return
            folder = self.run_dir / f"{step}_{modules[step]}"
            count = len(list(folder.glob("*.pdb"))) + len(list(folder.glob("*.pdb.gz"))) if folder.is_dir() else 0
            if count == last_count:
                return
            last_count = count
            total = expected[step]
            elapsed = last_report - step_start
            line = f"{label(step)}: {count}/{total or '?'} models | {_duration(elapsed)}"
            if total and 0 < count < total:
                line += f" | ~{_duration(elapsed / count * (total - count))} left"
            log_info(line)

        with open(log_path) as log:
            log.seek(offset)
            while True:
                done = process.poll() is not None
                for text in log.readlines():
                    if (match := started.search(text)) and match.group(1) in modules[step + 1:]:
                        step = modules.index(match.group(1), step + 1)
                        step_start, last_count = time.perf_counter(), -1
                        if modules[step] in HADDOCK_REFINE_MODULES and step > 0:
                            expected[step] = self._step_model_count(step - 1, params[step]) or expected[step]
                        log_info(f"{label(step)} started | run {_duration(step_start - run_start)}")
                    elif (match := finished.search(text)) and 0 <= step and match.group(1) == modules[step]:
                        report()
                        log_info(f"{label(step)} finished in {_duration(time.perf_counter() - step_start)}")
                if done:
                    return
                if time.perf_counter() - last_report >= HADDOCK_PROGRESS_SECONDS:
                    report()
                time.sleep(HADDOCK_POLL_SECONDS)

    def validate_config(self) -> Self:
        """Check the set up with 'haddock3 --setup' in a scratch run directory, removed afterwards.

        Returns:
            self.
        """
        self.export_inputs()
        scratch = self.workdir / "setup_check"
        config = self.workdir / "config_setup_check.toml"
        if scratch.exists():
            shutil.rmtree(scratch)
        config.write_text(self.config_toml(run_dir=scratch))
        try:
            self._run_command([str(config), "--setup"], self.workdir / "haddock3_setup_check.log")
        finally:
            if scratch.exists():
                shutil.rmtree(scratch)
            config.unlink(missing_ok=True)
        log_info(f"HADDOCK3: set up in {self.workdir} passed haddock3 validation.")
        return self

    def run_haddock(self, overwrite: bool = False, progress: bool = True) -> Self:
        """Write the set up and run haddock3.

        With progress, each module start and finish is logged, and while a module writes models the count written so
        far, with an estimate of the time left, is logged every few seconds when it changes. Interrupting Python (for
        example Ctrl+C or a notebook stop) stops haddock3.

        Args:
            overwrite: Whether to delete a previous run directory made by this class; others are never deleted.
            progress: Whether to log progress while haddock3 runs.

        When clustfcc is followed by seletopclusts without a top_clusters override, haddock3 runs up to clustfcc,
        the top fifth of the clusters (at least five, or all when fewer) is chosen, and the run restarts from
        seletopclusts.

        Returns:
            self.

        Raises:
            FileExistsError: If the run directory exists and may not be deleted.
        """
        start = time.perf_counter()
        self._top_clusters = None
        config = self.export_inputs()
        sentinel = self.workdir / RUN_SENTINEL
        if self.run_dir.exists() and any(self.run_dir.iterdir()):
            if not overwrite:
                raise FileExistsError(f"{self.run_dir} exists; pass overwrite=True to replace it.")
            if not sentinel.is_file():
                raise FileExistsError(f"{self.run_dir} was not created by HaddockDocking and is not deleted.")
            shutil.rmtree(self.run_dir)
        sentinel.write_text(f"{self.run_dir}\n")
        split = self._cluster_split()
        if split is None:
            self._run_command([str(config)], self.input_paths["log"], progress=progress)
        else:
            # top_clusters depends on how many clusters clustfcc finds: run up to clustfcc, count, then restart
            first = self.workdir / "config_clustering.toml"
            first.write_text(self.config_toml(stop=split))
            self._run_command([str(first)], self.input_paths["log"], progress=progress)
            self._top_clusters = self._top_cluster_count(split - 1)
            config.write_text(self.config_toml())
            self._run_command([str(config), "--restart", str(split)], self.input_paths["log"], progress=progress,
                              append=True)
        log_completion(start, f"HADDOCK3 '{self.method}' run finished in {self.run_dir}")
        return self

    def _cluster_split(self) -> int | None:
        """Return the index of a seletopclusts step whose top_clusters is chosen from the preceding clustfcc.

        Returns:
            The step index, or None when the workflow has no such step or top_clusters is given.
        """
        names = [name for name, _ in self.workflow_modules()]
        for index in range(1, len(names)):
            if names[index - 1:index + 1] == ["clustfcc", "seletopclusts"]:
                return None if "top_clusters" in self.module_params.get("seletopclusts", {}) else index
        return None

    def _top_cluster_count(self, step: int) -> int:
        """Return how many clusters to refine: a fixed fraction of the clusters found, at least a minimum.

        Clusters are read from the step's clustfcc.tsv, which needs no haddock3 import; clustfcc marks unclustered
        models with '-' and writes no table when it finds no clusters.

        Args:
            step: Number of the clustfcc step.

        Returns:
            Number of clusters for seletopclusts.

        Raises:
            RuntimeError: If clustfcc found no clusters.
        """
        path = self._step_path(step) / "clustfcc.tsv"
        clusters = 0
        if path.is_file():
            clusters = pd.read_csv(path, sep="\t", na_values="-")["cluster_id"].dropna().nunique()
        if clusters == 0:
            raise RuntimeError("clustfcc found no clusters; lower min_population with "
                               "module_params={'clustfcc': {'min_population': 1}}.")
        top = min(clusters, max(HADDOCK_MIN_TOP_CLUSTERS, math.ceil(HADDOCK_TOP_CLUSTER_FRACTION * clusters)))
        log_info(f"HADDOCK3: clustfcc found {clusters} clusters; refining the top {top}.")
        return top

###################################### RESULTS

    def step_paths(self) -> dict[str, Path]:
        """Return the step folders of the run directory.

        Returns:
            {step folder name: path}, in workflow order.

        Raises:
            FileNotFoundError: If there is no run directory.
        """
        if not self.run_dir.is_dir():
            raise FileNotFoundError(f"No run directory at {self.run_dir}; call run_haddock() first.")
        found = [(int(match.group(1)), path) for path in self.run_dir.iterdir()
                 if path.is_dir() and (match := STEP_PATTERN.match(path.name))]
        return {path.name: path for _, path in sorted(found)}

    def _step_path(self, step: StepRef) -> Path:
        """Find a step folder.

        Args:
            step: Folder name ('3_emref'), module name (its last step), or step number.

        Returns:
            The step folder.

        Raises:
            KeyError: If the run has no such step.
        """
        steps = self.step_paths()
        if str(step) in steps:
            return steps[str(step)]
        part = 0 if isinstance(step, (int, np.integer)) else 1
        matches = [path for name, path in steps.items() if _split_step_name(name)[part] == str(step)]
        if not matches:
            raise KeyError(f"No step {step!r} in {self.run_dir}; steps are {list(steps)}.")
        return matches[-1]

    def _caprieval_step_paths(self) -> list[Path]:
        """Return the caprieval step folders of the run, in workflow order.

        Returns:
            Caprieval step folders.
        """
        return [path for name, path in self.step_paths().items() if _split_step_name(name)[1] == "caprieval"]

    @staticmethod
    def _model_path(folder: PathLike, model: PathLike) -> Path:
        """Resolve a model path given relative to a folder.

        Args:
            folder: Folder the path is relative to.
            model: Model path as written by haddock3.

        Returns:
            Absolute path, of the .gz file when the model was compressed.
        """
        path = (Path(folder) / model).resolve()
        if not path.is_file() and Path(f"{path}.gz").is_file():
            path = Path(f"{path}.gz")
        return path

    def _read_capri_tables(self, folder: Path) -> dict[str, pd.DataFrame]:
        """Read the caprieval tables of a folder.

        Args:
            folder: Caprieval step or analysis folder.

        Returns:
            {'models': capri_ss} with an absolute model path, plus 'clusters' (capri_clt) when present.

        Raises:
            FileNotFoundError: If the folder has no capri_ss.tsv.
        """
        if not (folder / "capri_ss.tsv").is_file():
            raise FileNotFoundError(f"No capri_ss.tsv in {folder}.")
        models = pd.read_csv(folder / "capri_ss.tsv", sep="\t")
        models.insert(1, "path", [self._model_path(folder, model) for model in models["model"]])
        models.insert(2, "reference", _caprieval_reference(folder))
        tables = {"models": models.sort_values("caprieval_rank", ignore_index=True)}
        if (folder / "capri_clt.tsv").is_file():
            tables["clusters"] = pd.read_csv(folder / "capri_clt.tsv", sep="\t", comment="#")
        return tables

    def _read_io_rows(self, folder: Path) -> list[Record] | None:
        """Read one row per output model of a step from its io.json.

        Decoding needs jsonpickle and haddock3 importable, as io.json stores haddock3 objects.

        Args:
            folder: Step folder.

        Returns:
            Rows with model, path, score, clusters, origin, seed, and energies; None when io.json cannot be decoded.
        """
        if importlib.util.find_spec("jsonpickle") is None or importlib.util.find_spec("haddock") is None:
            return None
        import jsonpickle

        if not (folder / "io.json").is_file():
            return None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            data = jsonpickle.decode((folder / "io.json").read_text())
        rows: list[Record] = []
        for model in _flatten(data.get("output", [])):
            if not (hasattr(model, "file_name") and hasattr(model, "score")):
                continue
            rows.append({
                "step": folder.name, "model": model.file_name,
                "path": self._model_path(folder, model.rel_path), "score": model.score,
                "cluster_id": model.clt_id, "cluster_rank": model.clt_rank,
                "model_cluster_rank": model.clt_model_rank, "origin": model.ori_name, "seed": model.seed,
                **dict(getattr(model, "unw_energies", None) or {}),
            })
        return rows

    def _read_remark_rows(self, folder: Path) -> list[Record]:
        """Read one row per model file of a step from the HADDOCK REMARK header of each PDB.

        Args:
            folder: Step folder.

        Returns:
            Rows with model, path, score, and energy terms.
        """
        files = sorted(folder.glob("*.pdb")) + sorted(folder.glob("*.pdb.gz"))
        rows: list[Record] = []
        for path in files:
            row: Record = {"step": folder.name, "model": path.name.removesuffix(".gz"), "path": path.resolve(),
                           **_read_remark_scores(path)}
            if "score" in row:
                rows.append(row)
        return rows

    def score_table(self, step: StepRef | None = None) -> pd.DataFrame:
        """Return the HADDOCK scores of every model of a step, best first, with a model path and energy terms.

        Rows come from the step's io.json, which also gives clusters, origin, and seed; without haddock3 importable,
        from the REMARK header of the step's own PDB files, which gives the score and energy terms only.

        Args:
            step: Step folder name, module name (its last step), or step number; None for the last step with models.

        Returns:
            One row per model, with a 'rank' column after 'score'.

        Raises:
            FileNotFoundError: If no step has an io.json, or the step has no scored models.
        """
        if step is None:
            folders = [path for path in self.step_paths().values() if (path / "io.json").is_file()]
            if not folders:
                raise FileNotFoundError(f"No step of {self.run_dir} has an io.json.")
            folders = folders[::-1]
        else:
            folders = [self._step_path(step)]
        for folder in folders:
            rows = self._read_io_rows(folder)
            if rows is None:
                rows = self._read_remark_rows(folder)
            if rows:
                break
        else:
            raise FileNotFoundError(f"No scored models in {folders[0]}; pick a step that writes models.")
        frame = pd.DataFrame(rows)
        for column in HADDOCK_SCORE_COLUMNS:
            if column not in frame:
                frame[column] = None
        energies = [column for column in frame.columns if column not in HADDOCK_SCORE_COLUMNS]
        frame = frame[HADDOCK_SCORE_COLUMNS + energies].sort_values("score", ignore_index=True)
        frame.insert(list(frame.columns).index("score") + 1, "rank", range(1, len(frame) + 1))
        return frame

    def caprieval_tables(self, step: StepRef | None = None) -> dict[str, pd.DataFrame]:
        """Return the caprieval tables of a step.

        Args:
            step: None for the last caprieval step; a caprieval step; or a step analysed with run_caprieval().

        Returns:
            {'models': capri_ss} with absolute model paths, plus 'clusters' (capri_clt) when present.

        Raises:
            FileNotFoundError: If the run or step has no caprieval tables.
        """
        if step is None:
            folders = self._caprieval_step_paths()
            if not folders:
                raise FileNotFoundError("The run has no caprieval step; run with caprieval=True "
                                        "or use run_caprieval().")
            return self._read_capri_tables(folders[-1])
        folder = self._step_path(step)
        if _split_step_name(folder.name)[1] == "caprieval":
            return self._read_capri_tables(folder)
        analysed = self.run_dir / "analysis" / f"{folder.name}_analysis"
        if not (analysed / "capri_ss.tsv").is_file():
            raise FileNotFoundError(f"{folder.name} has no caprieval tables; "
                                    f"call run_caprieval(['{folder.name}']) first.")
        return self._read_capri_tables(analysed)

    def alascan_tables(self, step: StepRef | None = None) -> dict[str, pd.DataFrame]:
        """Return the alanine-scan results of an alascan step, mapped back to the original chains.

        delta_score is the wild-type HADDOCK score minus the mutant's. Lower scores are better, so negative values
        mark hot spots: residues whose mutation worsens the interface score.

        Args:
            step: Alascan step folder name or number; None for the last alascan step.

        Returns:
            {'models': one row per scanned residue of each model, and 'clusters': one row per scanned residue of each
            cluster, averaged over its models; most negative delta_score first within each model or cluster}.

        Raises:
            FileNotFoundError: If the run has no alascan step or the step wrote no scan tables.
        """
        if step is None:
            folders = [path for name, path in self.step_paths().items() if _split_step_name(name)[1] == "alascan"]
            if not folders:
                raise FileNotFoundError("The run has no alascan step; use method='alascan'.")
            folder = folders[-1]
        else:
            folder = self._step_path(step)
        models, clusters = [], []
        for path in sorted(folder.glob("scan_*.tsv")):
            frame = pd.read_csv(path, sep="\t", comment="#")
            if match := SCAN_CLUSTER_PATTERN.match(path.name):
                label = match.group(1)
                frame.insert(0, "cluster_id", int(label) if label.isdigit() else label)
                clusters.append(frame.drop(columns=["full_resname"], errors="ignore"))
            elif match := SCAN_MODEL_PATTERN.match(path.name):
                frame.insert(0, "model", match.group(1))
                models.append(frame.rename(columns={"res": "resid"}))
        if not models and not clusters:
            raise FileNotFoundError(f"No scan tables in {folder}; alascan found no interface residues to mutate.")
        return {"models": self._scan_frame(models, "model"), "clusters": self._scan_frame(clusters, "cluster_id")}

    def _scan_frame(self, frames: list[pd.DataFrame], group: str) -> pd.DataFrame:
        """Join alascan rows with the residue map and sort them.

        Args:
            frames: Scan tables, with 'chain' and 'resid' columns in segment IDs and HADDOCK3 numbers.
            group: Column the rows are grouped by, 'model' or 'cluster_id'.

        Returns:
            One row per scanned residue, most negative delta_score first within each group.
        """
        residues = self.residue_table()[HADDOCK_ALASCAN_COLUMNS]
        if not frames:
            return pd.DataFrame(columns=[group, *HADDOCK_ALASCAN_COLUMNS])
        scan = pd.concat(frames, ignore_index=True).rename(columns={"chain": "segid", "resid": "haddock_resid"})
        scan = scan.astype({"segid": str, "haddock_resid": int})
        merged = scan.merge(residues, on=["segid", "haddock_resid"], how="left")
        leading = [group, *HADDOCK_ALASCAN_COLUMNS]
        merged = merged[leading + [column for column in merged.columns if column not in leading]]
        return merged.sort_values([group, "delta_score"], ignore_index=True)

    def plot_caprieval(self, x: str = "irmsd", y: str = "score", color: str | None = "dockq",
                       step: StepRef | None = None, title: str | None = None,
                       ax: Axes | None = None) -> tuple[Figure | SubFigure, Axes]:
        """Plot two caprieval columns of the models against each other, colored by a third.

        The default, HADDOCK score against i-RMSD colored by DockQ, shows whether the best-scoring models are also
        the closest to the reference. Without a reference, caprieval measures models against the best-scoring
        model, and the title says so.

        Args:
            x: Column on the x-axis, such as 'irmsd', 'lrmsd', 'fnat', or 'dockq'.
            y: Column on the y-axis.
            color: Column mapped to the colorbar, or None for one color.
            step: Caprieval or analysed step, as in caprieval_tables(); None for the last caprieval step.
            title: Plot title, or None for one naming the reference.
            ax: Axes to draw into, or None for a new figure.

        Returns:
            (figure, axes).

        Raises:
            ValueError: If a column is missing or has no values.
        """
        models = self.caprieval_tables(step)["models"]
        columns = [x, y] + ([color] if color else [])
        for column in columns:
            if column not in models or pd.to_numeric(models[column], errors="coerce").isna().all():
                numeric = [name for name in models.columns
                           if pd.to_numeric(models[name], errors="coerce").notna().any()]
                raise ValueError(f"No values in caprieval column {column!r}; choose from {numeric}.")
        values = models[columns].apply(pd.to_numeric, errors="coerce").dropna()
        fig: Figure | SubFigure
        if ax is None:
            fig, ax = plt.subplots(figsize=(CAPRI_SIZE_IN, CAPRI_SIZE_IN), facecolor="white")
        else:
            fig = ax.figure
        ax.set_facecolor("white")
        ax.grid(True, linestyle="--", alpha=0.3, color=CAPRI_GRID_COLOR)
        ax.set_axisbelow(True)
        if color:
            points = ax.scatter(values[x], values[y], c=values[color], cmap=CAPRI_CMAP, s=CAPRI_MARKER_SIZE,
                                edgecolors="none", alpha=CAPRI_ALPHA)
            colorbar = fig.colorbar(points, ax=ax, pad=0.03, aspect=20)
            colorbar.outline.set_visible(False)
            colorbar.ax.tick_params(labelsize=CAPRI_TICK_PT - 1, colors=CAPRI_AXIS_COLOR)
            colorbar.set_label(CAPRI_LABELS.get(color, color), fontsize=CAPRI_LABEL_PT - 1, color=CAPRI_TEXT_COLOR)
        else:
            ax.scatter(values[x], values[y], color=CAPRI_MARKER_COLOR, s=CAPRI_MARKER_SIZE, edgecolors="none",
                       alpha=CAPRI_ALPHA)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(CAPRI_EDGE_COLOR)
            ax.spines[side].set_linewidth(0.8)
        ax.tick_params(colors=CAPRI_AXIS_COLOR, labelsize=CAPRI_TICK_PT)
        ax.set_xlabel(CAPRI_LABELS.get(x, x), fontsize=CAPRI_LABEL_PT, color=CAPRI_TEXT_COLOR)
        ax.set_ylabel(CAPRI_LABELS.get(y, y), fontsize=CAPRI_LABEL_PT, color=CAPRI_TEXT_COLOR)
        if title is None:
            reference = str(models["reference"].iloc[0]) if "reference" in models and len(models) else "best model"
            title = f"Models vs {Path(reference).stem if reference != 'best model' else reference}"
        ax.set_title(title, fontsize=CAPRI_TITLE_PT, fontweight="bold", pad=12)
        return fig, ax

    def docked_model_table(self, step: StepRef | None = None) -> pd.DataFrame:
        """Return the final models, best first: caprieval metrics, HADDOCK score and energies, and model path.

        The rows are the models that docked_model_paths() and load_docked_model() return. Without a caprieval step
        the table is score_table() of the last step with models.

        Args:
            step: Caprieval or analysed step, as in caprieval_tables(); None for the last caprieval step.

        Returns:
            One row per model.
        """
        if step is None and not self._caprieval_step_paths():
            return self.score_table()
        models = self.caprieval_tables(step)["models"]
        extra_rows = self._read_io_rows(self._step_path(step) if step is not None else self._caprieval_step_paths()[-1])
        if extra_rows:
            extra = pd.DataFrame(extra_rows)
            keep = ["path"] + [column for column in extra.columns
                               if column not in models.columns and column not in ("step", "model")]
            models = models.merge(extra[keep], on="path", how="left")
        return models

    def docked_model_paths(self, step: StepRef | None = None) -> list[Path]:
        """Return the final model paths, best first, as listed by docked_model_table().

        Args:
            step: Passed to docked_model_table().

        Returns:
            Model file paths.
        """
        return list(self.docked_model_table(step)["path"])

    def load_docked_models(self, top: int | None = None, step: StepRef | None = None) -> dict[str, SyncedStructure]:
        """Load the best models, in docked_model_table() order.

        Args:
            top: Number of models to load, or None for all.
            step: Passed to docked_model_table().

        Returns:
            {model name: SyncedStructure}.
        """
        paths = self.docked_model_paths(step)[:top]
        return {Path(path).name.split(".")[0]: self.load_docked_model(path) for path in paths}

    def run_caprieval(self, steps: StepRef | Iterable[StepRef],
                               ncores: int | None = None) -> dict[str, dict[str, pd.DataFrame]]:
        """Run caprieval on finished steps with haddock3-analyse.

        Models are compared with the reference structure when one is available, else with the best-scoring model.
        Compressed steps are unpacked for the analysis and compressed again by haddock3-analyse.

        Args:
            steps: Step folder names, module names, or step numbers, such as ['rigidbody', 'emref'].
            ncores: Number of cores, or None for the set up's ncores.

        Returns:
            {step folder name: caprieval_tables()}.

        Raises:
            RuntimeError: If haddock3-analyse produced no tables for a step.
        """
        selected = [steps] if isinstance(steps, (str, int)) else list(steps)
        folders = [self._step_path(step) for step in selected]
        cleaned = any(any(folder.glob("*.pdb.gz")) for folder in folders)
        for folder in folders:
            previous = self.run_dir / "analysis" / f"{folder.name}_analysis"
            if previous.is_dir():
                shutil.rmtree(previous)
        args = ["-r", str(self.run_dir), "-m", *[folder.name.split("_")[0] for folder in folders],
                "--ncores", str(ncores or self.ncores), "--is_cleaned", str(cleaned)]
        reference = self._reference_path()
        if reference is not None:
            args += ["-p", "reference_fname", str(reference)]
        log_path = self.workdir / "haddock3_analyse.log"
        self._run_command(args, log_path, name=HADDOCK_ANALYSE_EXECUTABLE)
        missing = [folder.name for folder in folders
                   if not (self.run_dir / "analysis" / f"{folder.name}_analysis" / "capri_ss.tsv").is_file()]
        if missing:
            raise RuntimeError(f"haddock3-analyse produced no tables for {missing}; see {log_path}\n"
                               f"{_log_tail(log_path)}")
        return {folder.name: self.caprieval_tables(folder.name) for folder in folders}

    def load_docked_model(self, model: int | PathLike = 0) -> SyncedStructure:
        """Load a docked model with the original chain keys, residue IDs, and reference sequences.

        Args:
            model: Index into docked_model_paths(), or a model path.

        Returns:
            The model as a SyncedStructure.
        """
        path = self.docked_model_paths()[model] if isinstance(model, int) else Path(model)
        parsed = next(_read_model(path).get_models())
        structure = Structure(path.name.split(".")[0])
        rebuilt = Model(0)
        structure.add(rebuilt)
        keys: dict[ChainTarget, str] = {}
        missing = 0
        for (index, key), entries in self._chains.items():
            segid = self.partners[index]["segid"]
            new_key = key if key not in keys.values() else f"{key}{segid}"
            keys[(index, key)] = new_key
            chain = Chain(new_key)
            for entry in entries:
                residue_id = (" ", entry["resid"], " ")
                if segid not in parsed or residue_id not in parsed[segid]:
                    missing += 1
                    continue
                residue = parsed[segid][residue_id]
                parsed[segid].detach_child(residue_id)
                residue.id = entry["residue"].id
                residue.resname = entry["residue"].get_resname()
                chain.add(residue)
            rebuilt.add(chain)
        if missing:
            log_info(f"HADDOCK3: {missing} residue(s) of the input are missing from {path.name}.")
        first = self.partners[0]["structure"]
        loaded = type(first).from_biopython(structure, structure_id=structure.id, file_format="pdb")
        for (index, key), new_key in keys.items():
            source = self.partners[index]["structure"]
            record = source.chains[key]
            loaded.set_reference_sequence(new_key, record.reference_sequence)
            if record.entity_id is not None:
                loaded.set_entity_id(new_key, record.entity_id)
            loaded.set_accession(new_key, source.get_accession(key))
        return loaded

###################################### STRUCTURE MIXIN


class DockingMixin(StructureState):
    """Set up HADDOCK3 docking of this structure's chains."""

    def prepare_haddock_docking(self, *partners: str | Sequence[str] | PartnerSpec, workdir: PathLike,
                                **kwargs: Any) -> HaddockDocking:
        """Return a HaddockDocking of the given partners.

        Args:
            *partners: Chains of this structure (one chain or a list per partner), or (SyncedStructure, chains).
            workdir: Directory for the HADDOCK3 inputs and run.
            **kwargs: Passed to HaddockDocking, such as method, seed, ncores, or module_params.

        Returns:
            The docking set up.
        """
        host = cast("SyncedStructure", self)
        normalized: list[PartnerSpec] = []
        for item in partners:
            if isinstance(item, (tuple, list)) and len(item) == 2 and hasattr(item[0], "chains"):
                normalized.append(cast(PartnerSpec, tuple(item)))
            else:
                normalized.append((host, cast("str | Sequence[str]", item)))
        return HaddockDocking(normalized, workdir, **kwargs)


# imported last because structure.py imports this module; type checks resolve it when functions run
from . import structure as _structure  # noqa: E402
