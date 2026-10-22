"""Sequence cleaning, residue letters, and mmCIF field helpers."""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import Any

from Bio.PDB.Residue import Residue

from .constants import MOIETY_NUCLEIC, MOIETY_PROTEIN, NUCLEIC_TO1, PROTEIN_3TO1

###################################### SEQUENCE HELPERS


def clean_sequence(sequence: Any) -> str:
    """Return an uppercase one-letter sequence with non-alphabetic characters replaced by X.

    Args:
        sequence: Sequence text; whitespace is removed.

    Returns:
        The cleaned sequence.
    """
    sequence = "".join(str(sequence).split()).upper()
    return "".join(char if char.isalpha() else "X" for char in sequence)


def clean_cif_sequence(sequence: Any) -> str:
    """Return a cleaned mmCIF sequence with parenthesized modified residues replaced by X.

    Args:
        sequence: mmCIF one-letter sequence.

    Returns:
        The cleaned sequence.
    """
    return clean_sequence(re.sub(r"\([^)]*\)", "X", str(sequence)))


def one_letter_code(residue: Residue, moiety: str) -> str:
    """Return the one-letter code of a protein or nucleic residue.

    Args:
        residue: Biopython Residue.
        moiety: Biomolecule type from residue_class.

    Returns:
        The one-letter code, or X if unknown.
    """
    table = {MOIETY_PROTEIN: PROTEIN_3TO1, MOIETY_NUCLEIC: NUCLEIC_TO1}.get(moiety, {})
    return table.get(residue.get_resname().strip().upper(), "X")

###################################### MMCIF FIELDS


def to_list(value: str | Iterable[Any] | None) -> list[str]:
    """Return an mmCIF field as a list of strings.

    Args:
        value: Field value: a string, a list of values, or None.

    Returns:
        The values as strings; empty if value is None.
    """
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def clean_cif_value(value: str | None) -> str | None:
    """Return None for mmCIF missing-value markers ('?' or '.'), otherwise the value unchanged.

    Args:
        value: mmCIF field value.

    Returns:
        The value, or None for a missing-value marker.
    """
    return None if value in ("?", ".", None) else value
