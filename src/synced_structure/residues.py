"""Residue classification, numbering, reference-index mapping, and per-residue tables with standard and author chain IDs."""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

import copy

import numpy as np
import pandas as pd
from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from .constants import (
    MAPPING_COLUMNS, MOIETIES, MOIETY_NUCLEIC, MOIETY_OTHER, MOIETY_PROTEIN, MOIETY_SUMMARY_COLUMNS,
    NUCLEIC_RESNAMES, POLYMER_MOIETIES, PROTEIN_RESNAMES, REPRESENTATIVE_ATOMS, RESIDUE_COLUMNS,
)
from ._types import ChainRefs, Self
from .base import StructureState

ResidueMoiety = tuple[Residue, str]

###################################### RESIDUE CLASSIFICATION


def residue_class(residue: Residue) -> str:
    """Return the biomolecule type of a residue.

    Single-atom residues are 'other' before the nucleic-acid names are checked, so a uranium ion ('U') is not
    read as uracil.

    Args:
        residue: Biopython Residue.

    Returns:
        'protein', 'nucleic', or 'other' (ligands, ions, and water).
    """
    name = residue.get_resname().strip().upper()
    if residue.id[0] == "W":
        return MOIETY_OTHER
    if name in PROTEIN_RESNAMES:
        return MOIETY_PROTEIN
    if len(residue) == 1:
        return MOIETY_OTHER
    if name in NUCLEIC_RESNAMES:
        return MOIETY_NUCLEIC
    return MOIETY_OTHER


def is_heavy_atom(atom: Atom | DisorderedAtom) -> bool:
    """Return True if an atom is not hydrogen.

    Args:
        atom: Biopython Atom, or a DisorderedAtom (alternate locations).

    Returns:
        Whether the atom is a heavy atom.
    """
    return (atom.element or atom.get_name()[:1]).strip().upper() != "H"


def contact_point_coord(residue: Residue, moiety: str) -> np.ndarray | None:
    """Return a residue's contact point in A.

    Protein uses C-alpha, nucleic acids C4' (or C1', then P), other residues the heavy-atom centroid.

    Args:
        residue: Biopython Residue.
        moiety: Biomolecule type from residue_class.

    Returns:
        The contact point, or None if the residue has none.
    """
    if moiety == MOIETY_OTHER:
        coords = [atom.coord for atom in residue.get_atoms() if is_heavy_atom(atom)]
        return np.mean(coords, axis=0) if coords else None
    for name in REPRESENTATIVE_ATOMS.get(moiety, ()):
        if name in residue:
            return residue[name].coord
    return None


def split_residues(chain: Chain) -> tuple[list[ResidueMoiety], list[ResidueMoiety]]:
    """Split a chain's residues into polymer and hetero residues, in chain order.

    Args:
        chain: Biopython Chain.

    Returns:
        (polymer, hetero) lists of (residue, moiety) pairs.
    """
    polymer: list[ResidueMoiety] = []
    hetero: list[ResidueMoiety] = []
    for residue in chain.get_residues():
        moiety = residue_class(residue)
        (polymer if moiety in POLYMER_MOIETIES else hetero).append((residue, moiety))
    return polymer, hetero


def _rebuild_residue_index(chain: Chain) -> None:
    """Rebuild a chain's ID index after its residue IDs change, keeping residue order.

    Args:
        chain: Biopython Chain whose residue IDs changed.
    """
    chain.child_dict = {child.id: child for child in chain.child_list}


###################################### RESIDUE MAPPING


class ResidueMixin(StructureState):
    """Residue tables, reference numbering, and trimming of unmapped residues."""

    def _chain_residue_table(self, chain_id: str) -> pd.DataFrame:
        """Return one row per polymer residue of a chain with 1-based sequence and reference indices.

        Args:
            chain_id: Chain key as held by the structure.

        Returns:
            Residue table with RESIDUE_COLUMNS.
        """
        data = self.chains[chain_id]
        mapping = data.sequence_index_to_ref_index
        rows = [
            {
                "structure_id": self.structure_id,
                "chain_id": data.chain_id,
                "auth_chain_id": data.auth_chain_id,
                "pdb_resnum": residue.id[1],
                "resname": residue.resname.strip().upper(),
                "sequence_index": index + 1,
                "reference_index": mapping[index] + 1 if index in mapping else None,
                "moiety": moiety,
            }
            for index, (residue, moiety) in enumerate(zip(data.residues, data.moieties))
        ]
        frame = pd.DataFrame(rows, columns=RESIDUE_COLUMNS)
        frame["reference_index"] = frame["reference_index"].astype("Int64")
        return frame

    def chain_residue_table(self, chain_id: str) -> pd.DataFrame:
        """Return the author number, sequence index, and reference index of each residue of a chain.

        Args:
            chain_id: Chain key, label ID, or author ID.

        Returns:
            Table with standard and author chain IDs, pdb_resnum, sequence_index, and reference_index.
        """
        return self._chain_residue_table(self._chain_key(chain_id))[MAPPING_COLUMNS].reset_index(drop=True)

    @property
    def residue_table(self) -> pd.DataFrame:
        """Return all polymer residues across chains with their indices.

        Returns:
            Residue table with RESIDUE_COLUMNS.
        """
        if not self.chains:
            return pd.DataFrame(columns=RESIDUE_COLUMNS)
        return pd.concat([self._chain_residue_table(chain_id) for chain_id in self.chains], ignore_index=True)

    @property
    def moiety_table(self) -> pd.DataFrame:
        """Return residue and atom counts per chain and biomolecule type.

        Returns:
            One row per chain and moiety.
        """
        rows = []
        for data in self.chains.values():
            counts: dict[str, tuple[int, int]] = {}
            pairs = list(zip(data.residues, data.moieties)) + data.hetero
            for residue, moiety in pairs:
                n_res, n_atoms = counts.get(moiety, (0, 0))
                counts[moiety] = (n_res + 1, n_atoms + len(residue))
            for moiety, (n_res, n_atoms) in counts.items():
                rows.append({
                    "structure_id": self.structure_id,
                    "chain_id": data.chain_id,
                    "auth_chain_id": data.auth_chain_id,
                    "moiety": moiety,
                    "residues": n_res,
                    "atoms": n_atoms,
                })
        return pd.DataFrame(rows, columns=MOIETY_SUMMARY_COLUMNS)

    def chain_ids_by_moiety(self, moiety: str) -> list[str]:
        """Return the chain IDs that contain at least one residue of a biomolecule type.

        Args:
            moiety: One of 'protein', 'nucleic', or 'other'.

        Returns:
            Label chain IDs, in chain order.

        Raises:
            ValueError: If moiety is not a known biomolecule type.
        """
        if moiety not in MOIETIES:
            raise ValueError(f"moiety must be one of {MOIETIES}.")
        frame = self.moiety_table
        return frame.loc[frame["moiety"] == moiety, "chain_id"].drop_duplicates().tolist()

###################################### RENUMBERING AND TRIMMING

    @staticmethod
    def _set_residue_ids(chain: Chain, pairs: list[tuple[Residue, tuple[str, int, str]]]) -> None:
        """Set residue IDs from (residue, new_id) pairs in two phases to avoid transient collisions.

        Args:
            chain: Biopython Chain whose residues are renumbered.
            pairs: (residue, new_id) tuples.
        """
        for offset, (residue, _) in enumerate(pairs):
            residue.id = ("#", offset, " ")
        _rebuild_residue_index(chain)
        for residue, new_id in pairs:
            residue.id = new_id
        _rebuild_residue_index(chain)

    @staticmethod
    def _renumber_polymers(chain: Chain) -> None:
        """Renumber a chain's polymer residues 1..N in place, leaving hetero residues untouched.

        Args:
            chain: Biopython Chain to renumber.
        """
        polymer, _ = split_residues(chain)
        ResidueMixin._set_residue_ids(chain, [(residue, (" ", number, " "))
                                       for number, (residue, _) in enumerate(polymer, start=1)])

    @classmethod
    def _trim_to_mapped(cls, chain: Chain, index_mapping: dict[int, int]) -> None:
        """Remove unmapped polymer residues and renumber the kept ones by reference position; hetero residues stay.

        Args:
            chain: Biopython Chain to trim in place.
            index_mapping: 0-based observed index to 0-based reference index.
        """
        polymer = [residue for residue, _ in split_residues(chain)[0]]
        pairs = [(polymer[index], (" ", ref_index + 1, " "))
                 for index, ref_index in index_mapping.items() if index < len(polymer)]
        kept = {id(residue) for residue, _ in pairs}
        for residue in polymer:
            if id(residue) not in kept:
                chain.detach_child(residue.id)
        cls._set_residue_ids(chain, pairs)

    def apply_reference_numbering(self, chain_ids: ChainRefs | None = None) -> Self:
        """Trim unmapped polymer residues and renumber mapped ones by reference position, in place.

        Args:
            chain_ids: Chains to process, or None for all chains.

        Returns:
            This structure.
        """
        selected = self.chain_ids if chain_ids is None else [self._chain_key(c) for c in chain_ids]
        trimmed = []
        for chain_id in selected:
            mapping = self.chains[chain_id].sequence_index_to_ref_index
            if mapping:
                self._trim_to_mapped(self.chains[chain_id].chain, mapping)
                trimmed.append(chain_id)
        self._build_chains()
        for chain_id in trimmed:
            data = self.chains[chain_id]
            data.sequence_index_to_ref_index = {index: residue.id[1] - 1 for index, residue in enumerate(data.residues)}
        return self

    def renumber_chain(self, chain_id: str, in_place: bool = False) -> Chain:
        """Renumber a chain's polymer residues 1..N.

        Args:
            chain_id: Chain key, label ID, or author ID.
            in_place: Whether to renumber the stored chain itself instead of a copy.

        Returns:
            The renumbered chain.
        """
        chain_id = self._chain_key(chain_id)
        chain = self.chains[chain_id].chain if in_place else copy.deepcopy(self.chains[chain_id].chain)
        self._renumber_polymers(chain)
        if in_place:
            self._build_chains()
        return chain
