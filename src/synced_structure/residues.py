"""Residue classification, numbering, reference-index mapping, and per-residue tables with standard and author chain IDs."""

###################################### IMPORTS AND CONSTANTS

import copy

import numpy as np
import pandas as pd

from .constants import (
    ION_RESNAMES, MAPPING_COLUMNS, MOIETIES, MOIETY_ION, MOIETY_LIGAND, MOIETY_NUCLEIC, MOIETY_PROTEIN,
    MOIETY_SUMMARY_COLUMNS, MOIETY_WATER, NUCLEIC_RESNAMES, POLYMER_MOIETIES, PROTEIN_RESNAMES,
    REPRESENTATIVE_ATOMS, RESIDUE_COLUMNS, WATER_RESNAMES,
)

###################################### RESIDUE CLASSIFICATION


def residue_class(residue):
    """Return the biomolecule type of a residue: 'protein', 'nucleic', 'ligand', 'ion', or 'water'.

    Single-atom residues are ions before the nucleic-acid names are checked, so a uranium ion ('U') is not
    read as uracil.
    """
    name = residue.get_resname().strip().upper()
    if residue.id[0] == "W" or name in WATER_RESNAMES:
        return MOIETY_WATER
    if name in PROTEIN_RESNAMES:
        return MOIETY_PROTEIN
    if len(residue) == 1:
        return MOIETY_ION
    if name in NUCLEIC_RESNAMES:
        return MOIETY_NUCLEIC
    if name in ION_RESNAMES:
        return MOIETY_ION
    return MOIETY_LIGAND


residue_moiety = residue_class  # deprecated name, kept for existing callers


def is_heavy(atom):
    """Return True if an atom is not hydrogen."""
    return (atom.element or atom.get_name()[:1]).strip().upper() != "H"


def representative_coord(residue, moiety):
    """Return a residue's contact point in A, or None if it has none.

    Protein uses C-alpha, nucleic acids C4' (or C1', then P), ligands and ions the heavy-atom centroid.

    Args:
        residue: Biopython Residue.
        moiety: Biomolecule type from residue_class.
    """
    if moiety in (MOIETY_LIGAND, MOIETY_ION):
        coords = [atom.coord for atom in residue.get_atoms() if is_heavy(atom)]
        return np.mean(coords, axis=0) if coords else None
    for name in REPRESENTATIVE_ATOMS.get(moiety, ()):
        if name in residue:
            return residue[name].coord
    return None


def split_residues(chain):
    """Return (polymer, hetero) lists of (residue, moiety) pairs for a chain, in chain order."""
    polymer, hetero = [], []
    for residue in chain.get_residues():
        moiety = residue_class(residue)
        (polymer if moiety in POLYMER_MOIETIES else hetero).append((residue, moiety))
    return polymer, hetero

def _reindex_children(chain):
    """Rebuild a Biopython entity's ID index after its residue IDs change, keeping residue order."""
    chain.child_dict = {child.id: child for child in chain.child_list}


###################################### RESIDUE MAPPING


class ResidueMixin:
    """Residue tables, reference numbering, and trimming of unmapped residues."""

    def _residue_frame(self, chain_id):
        """Return one row per polymer residue with 1-based sequence and reference indices."""
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

    def get_residue_mapping(self, chain_id):
        """Return a DataFrame with standard and author chain IDs, pdb_resnum, sequence_index, and reference_index."""
        return self._residue_frame(self._require_chain(chain_id))[MAPPING_COLUMNS].reset_index(drop=True)

    @property
    def residue_map(self):
        """Return a DataFrame of all polymer residues across chains with their indices."""
        if not self.chains:
            return pd.DataFrame(columns=RESIDUE_COLUMNS)
        return pd.concat([self._residue_frame(chain_id) for chain_id in self.chains], ignore_index=True)

    @property
    def moiety_summary(self):
        """Return one row per chain and biomolecule type with residue and atom counts."""
        rows = []
        for data in self.chains.values():
            counts = {}
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

    def chains_by_moiety(self, moiety):
        """Return the chain IDs that contain at least one residue of a biomolecule type.

        Args:
            moiety: One of 'protein', 'nucleic', 'ligand', 'ion', or 'water'.
        """
        if moiety not in MOIETIES:
            raise ValueError(f"moiety must be one of {MOIETIES}.")
        frame = self.moiety_summary
        return frame.loc[frame["moiety"] == moiety, "chain_id"].drop_duplicates().tolist()

###################################### RENUMBERING AND TRIMMING

    @staticmethod
    def _set_residue_ids(chain, pairs):
        """Set residue IDs from (residue, new_id) pairs in two phases to avoid transient collisions.

        Args:
            chain: Biopython Chain whose residues are renumbered.
            pairs: List of (residue, new_id) tuples.
        """
        for offset, (residue, _) in enumerate(pairs):
            residue.id = ("#", offset, " ")
        _reindex_children(chain)
        for residue, new_id in pairs:
            residue.id = new_id
        _reindex_children(chain)

    @staticmethod
    def _renumber(chain):
        """Renumber a chain's polymer residues 1..N in place, leaving hetero residues untouched."""
        polymer, _ = split_residues(chain)
        ResidueMixin._set_residue_ids(chain, [(residue, (" ", number, " "))
                                       for number, (residue, _) in enumerate(polymer, start=1)])

    @classmethod
    def _keep_mapped(cls, chain, index_mapping):
        """Remove unmapped polymer residues and renumber the kept ones by reference position; hetero residues stay.

        Args:
            chain: Biopython Chain to trim in place.
            index_mapping: Dict from 0-based observed to 0-based reference indices.
        """
        polymer = [residue for residue, _ in split_residues(chain)[0]]
        pairs = [(polymer[index], (" ", ref_index + 1, " "))
                 for index, ref_index in index_mapping.items() if index < len(polymer)]
        kept = {id(residue) for residue, _ in pairs}
        for residue in polymer:
            if id(residue) not in kept:
                chain.detach_child(residue.id)
        cls._set_residue_ids(chain, pairs)

    def apply_reference_numbering(self, chain_ids=None):
        """Trim unmapped polymer residues and renumber mapped ones by reference position in place, and return self.

        Args:
            chain_ids: Chains to process, or None for all chains.
        """
        selected = self.chain_ids if chain_ids is None else [self._require_chain(c) for c in chain_ids]
        trimmed = []
        for chain_id in selected:
            mapping = self.chains[chain_id].sequence_index_to_ref_index
            if mapping:
                self._keep_mapped(self.chains[chain_id].chain, mapping)
                trimmed.append(chain_id)
        self._build_chains()
        for chain_id in trimmed:
            data = self.chains[chain_id]
            data.sequence_index_to_ref_index = {index: residue.id[1] - 1 for index, residue in enumerate(data.residues)}
        return self

    def reindex_chain(self, chain_id, in_place=False):
        """Return a chain with polymer residues renumbered 1..N, modifying the structure only if in_place is True.

        Args:
            chain_id: Chain identifier.
            in_place: Whether to renumber the stored chain itself.
        """
        chain_id = self._require_chain(chain_id)
        chain = self.chains[chain_id].chain if in_place else copy.deepcopy(self.chains[chain_id].chain)
        self._renumber(chain)
        if in_place:
            self._build_chains()
        return chain
