#!/usr/bin/env python
"""Helpers for configuration, logging, sequence cleaning, aligners, and contact arguments."""

###################################### IMPORTS AND CONSTANTS

import datetime
import re
import resource
import time
from pathlib import Path

import yaml
from Bio.Align import PairwiseAligner
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1_extended

from constants import GLOBAL_GAP_EXTEND_SCORE, GLOBAL_GAP_OPEN_SCORE, GLOBAL_MATCH_SCORE, GLOBAL_MISMATCH_SCORE

###################################### CONFIG AND LOGGING


def resolve_config(yaml_arg, defaults, section):
    """Return a section of a YAML config, filling missing keys from defaults and saving the result.

    Args:
        yaml_arg: Path to the YAML file, or None to use defaults only.
        defaults: Dict of default values for the section.
        section: Top-level YAML key that holds this section.

    Returns:
        Dict with the merged section values.
    """
    if yaml_arg is None:
        return dict(defaults)
    path = Path(yaml_arg)
    data = (yaml.safe_load(path.read_text()) if path.exists() else None) or {}
    config = {**defaults, **(data.get(section) or {})}
    if data.get(section) != config:
        data[section] = config
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(data, sort_keys=False))
    return config


def finish_log(start_time, message):
    """Print a timestamped log line with elapsed time and peak memory.

    Args:
        start_time: Value of time.perf_counter() taken at the start of the task.
        message: Human-readable description of the completed task.
    """
    elapsed = time.perf_counter() - start_time
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    mem_mb = usage / 1024 if usage < 10**7 else usage / 1024**2
    timestamp = datetime.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[PYTHON-INFO] {timestamp} | {message} | execution time: {elapsed:.2f}s | max memory: {mem_mb:.1f} MB")

###################################### SEQUENCE HELPERS


def clean_sequence(sequence):
    """Return an uppercase one-letter sequence with non-alphabetic characters replaced by X."""
    sequence = "".join(str(sequence).split()).upper()
    return "".join(char if char.isalpha() else "X" for char in sequence)


def cif_sequence(sequence):
    """Return a cleaned mmCIF sequence with parenthesized modified residues replaced by X."""
    sequence = re.sub(r"\([^)]*\)", "X", str(sequence))
    return clean_sequence(sequence)


def to_list(value):
    """Return an mmCIF field as a list of strings.

    Args:
        value: None, a single string, or an iterable of values.

    Returns:
        List of string values; empty if value is None.
    """
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def residue_to_letter(residue):
    """Return the one-letter code of a residue, or X if unknown."""
    return protein_letters_3to1_extended.get(residue.resname.upper(), "X")


def is_protein_residue(residue):
    """Return True if the residue is recognized as an amino acid."""
    return is_aa(residue, standard=False)


def make_global_aligner():
    """Return a global PairwiseAligner with the fixed scoring scheme."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = GLOBAL_MATCH_SCORE
    aligner.mismatch_score = GLOBAL_MISMATCH_SCORE
    aligner.open_gap_score = GLOBAL_GAP_OPEN_SCORE
    aligner.extend_gap_score = GLOBAL_GAP_EXTEND_SCORE
    return aligner


def make_local_aligner(config):
    """Return a local PairwiseAligner scored from the config dict."""
    return PairwiseAligner(
        mode="local",
        match_score=config["aln_match_score"],
        mismatch_score=config["aln_mismatch_score"],
        open_gap_score=config["aln_open_gap_score"],
        extend_gap_score=config["aln_extend_gap_score"],
    )

###################################### CONTACT ARGUMENTS


def _check_contact_args(cutoff, contact_type):
    """Raise ValueError if the cutoff is not positive or the contact type is unknown."""
    if cutoff <= 0:
        raise ValueError("cutoff must be greater than zero.")
    if contact_type not in {"inter", "intra", "all"}:
        raise ValueError("contact_type must be 'inter', 'intra', or 'all'.")


def _keep_pair(chain_a, chain_b, contact_type):
    """Return True if a residue pair on these chains matches the requested contact type."""
    if contact_type == "all":
        return True
    same_chain = chain_a == chain_b
    return (contact_type == "intra") == same_chain
