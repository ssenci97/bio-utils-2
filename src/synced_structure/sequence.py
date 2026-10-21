"""Sequence cleaning, residue letters, and mmCIF field helpers."""

###################################### IMPORTS AND CONSTANTS

import re

from .constants import MOIETY_NUCLEIC, MOIETY_PROTEIN, NUCLEIC_TO1, TO1

###################################### SEQUENCE HELPERS


def clean_sequence(sequence):
    """Return an uppercase one-letter sequence with non-alphabetic characters replaced by X."""
    sequence = "".join(str(sequence).split()).upper()
    return "".join(char if char.isalpha() else "X" for char in sequence)


def cif_sequence(sequence):
    """Return a cleaned mmCIF sequence with parenthesized modified residues replaced by X."""
    return clean_sequence(re.sub(r"\([^)]*\)", "X", str(sequence)))


def residue_to_letter(residue, moiety):
    """Return the one-letter code of a protein or nucleic residue, or X if unknown.

    Args:
        residue: Biopython Residue.
        moiety: Biomolecule type from residue_moiety.
    """
    table = {MOIETY_PROTEIN: TO1, MOIETY_NUCLEIC: NUCLEIC_TO1}.get(moiety, {})
    return table.get(residue.get_resname().strip().upper(), "X")

###################################### MMCIF FIELDS


def to_list(value):
    """Return an mmCIF field as a list of strings, empty if value is None."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def clean_val(value):
    """Return None for mmCIF missing-value markers ('?' or '.'), otherwise the value unchanged."""
    return None if value in ("?", ".", None) else value
