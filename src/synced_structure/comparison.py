"""Structure comparison against a reference: residue matching, Kabsch superposition, RMSD, and CAPRI/DockQ metrics."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd

from ._types import ChainRefs, PathLike, Record, Self
from .alignment import sequence_index_map
from .base import StructureState
from .constants import (
    CAPRI_CLASSES, CAPRI_COLUMNS, CAPRI_CONTACT_A, CAPRI_INCORRECT, CAPRI_INTERFACE_A, DOCKQ_IRMSD_SCALE_A,
    DOCKQ_LRMSD_SCALE_A, RESIDUE_MATCH_COLUMNS, SUPERPOSITION_ATOMS,
)
from .interface_metrics import InterfaceMixin
from .residues import is_heavy_atom
from .utils import log_info

if TYPE_CHECKING:
    from Bio.PDB.Residue import Residue

    from .structure import SyncedStructure

# (model chain key, model sequence index, reference chain key, reference sequence index), 0-based
ResiduePair = tuple[str, int, str, int]

###################################### GEOMETRY


def _kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the rotation and translation that best fit mobile onto target, in least squares.

    Args:
        mobile: (n, 3) coordinates to move.
        target: (n, 3) coordinates to fit onto.

    Returns:
        (rotation, translation) such that mobile @ rotation + translation fits target, the Biopython convention.

    Raises:
        ValueError: If the sets differ in shape or hold fewer than three points.
    """
    if mobile.shape != target.shape or len(mobile) < 3:
        raise ValueError(f"Superposition needs two equal sets of at least 3 atoms, got {len(mobile)} and "
                         f"{len(target)}.")
    mobile_center, target_center = mobile.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - mobile_center).T @ (target - target_center))
    flip = np.diag([1.0, 1.0, float(np.sign(np.linalg.det(u @ vt)) or 1.0)])
    rotation = u @ flip @ vt
    return rotation, target_center - mobile_center @ rotation


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Return the root mean square deviation of two equal coordinate sets.

    Args:
        a: (n, 3) coordinates.
        b: (n, 3) coordinates.

    Returns:
        The RMSD, or NaN for empty sets.
    """
    if not len(a):
        return float("nan")
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def _atom_names(residue: Residue, moiety: str, atoms: str) -> list[str]:
    """Return the names of the atoms of a residue used for an atom selection.

    Args:
        residue: Biopython residue.
        moiety: Moiety of the residue.
        atoms: Atom selection, a key of SUPERPOSITION_ATOMS.

    Returns:
        Atom names in residue order.
    """
    if atoms == "heavy":
        return [atom.get_name() for atom in residue if is_heavy_atom(atom)]
    wanted = SUPERPOSITION_ATOMS[atoms].get(moiety, ())
    return [name for name in wanted if name in residue]


def _chain_list(chains: str | ChainRefs) -> list[str]:
    """Return chain references as a list, treating a single string as one chain.

    Args:
        chains: One chain reference or several.

    Returns:
        The chain references.
    """
    return [chains] if isinstance(chains, str) else list(chains)


def capri_class(fnat: float, irmsd: float, lrmsd: float) -> str:
    """Return the CAPRI quality class of a model.

    Args:
        fnat: Fraction of native contacts.
        irmsd: Interface RMSD in Å.
        lrmsd: Ligand RMSD in Å.

    Returns:
        'high', 'medium', 'acceptable', or 'incorrect'.
    """
    for name, min_fnat, max_lrmsd, max_irmsd in CAPRI_CLASSES:
        if fnat >= min_fnat and (lrmsd <= max_lrmsd or irmsd <= max_irmsd):
            return name
    return CAPRI_INCORRECT


def dockq_score(fnat: float, irmsd: float, lrmsd: float) -> float:
    """Return the DockQ score of a model, from 0 (incorrect) to 1 (identical to the reference).

    Args:
        fnat: Fraction of native contacts.
        irmsd: Interface RMSD in Å.
        lrmsd: Ligand RMSD in Å.

    Returns:
        The DockQ score.
    """
    return (fnat + 1 / (1 + (irmsd / DOCKQ_IRMSD_SCALE_A) ** 2) + 1 / (1 + (lrmsd / DOCKQ_LRMSD_SCALE_A) ** 2)) / 3

###################################### COMPARISON MIXIN


class ComparisonMixin(StructureState):
    """Comparison of this structure (the model) with a reference: residue matching, superposition, and CAPRI metrics."""

    def _chain_pairs(self, reference: SyncedStructure, chain_map: Mapping[str, str] | None,
                     chains: str | ChainRefs | None = None) -> list[tuple[str, str]]:
        """Resolve which model chain is compared with which reference chain.

        Args:
            reference: Reference structure.
            chain_map: {model chain: reference chain}, or None to pair chains with the same key.
            chains: Model chains to keep, or None for every paired chain.

        Returns:
            (model chain key, reference chain key) pairs in model chain order.

        Raises:
            ValueError: If no chain pairs, or a requested chain has no counterpart.
        """
        if chain_map is None:
            pairs = {key: key for key in self.chains if key in reference.chains}
        else:
            pairs = {self._chain_key(model): reference._chain_key(ref) for model, ref in chain_map.items()}
        if chains is not None:
            keys = [self._chain_key(chain) for chain in _chain_list(chains)]
            missing = [key for key in keys if key not in pairs]
            if missing:
                raise ValueError(f"Chains {missing} of {self.structure_id} have no counterpart in "
                                 f"{reference.structure_id}; pass chain_map=.")
            pairs = {key: pairs[key] for key in keys}
        if not pairs:
            raise ValueError(f"No chain of {self.structure_id} matches a chain of {reference.structure_id}; "
                             "pass chain_map={model chain: reference chain}.")
        order = {key: index for index, key in enumerate(self.chains)}
        return sorted(pairs.items(), key=lambda pair: order[pair[0]])

    def _residue_pairs(self, reference: SyncedStructure, chain_map: Mapping[str, str] | None,
                       chains: str | ChainRefs | None = None) -> list[ResiduePair]:
        """Return the residues of this structure matched to residues of the reference.

        Args:
            reference: Reference structure.
            chain_map: {model chain: reference chain}, or None to pair chains with the same key.
            chains: Model chains to keep, or None for every paired chain.

        Returns:
            Matched residues, 0-based, in model chain and sequence order.
        """
        pairs: list[ResiduePair] = []
        for key, ref_key in self._chain_pairs(reference, chain_map, chains):
            mapping = sequence_index_map(self.chains[key], reference.chains[ref_key])
            pairs.extend((key, seq, ref_key, mapping[seq]) for seq in sorted(mapping))
        return pairs

    def _paired_coords(self, reference: SyncedStructure, pairs: Sequence[ResiduePair],
                       atoms: str) -> tuple[np.ndarray, np.ndarray]:
        """Return the coordinates of the atoms shared by matched residue pairs.

        Args:
            reference: Reference structure.
            pairs: Matched residues.
            atoms: Atom selection: 'ca', 'backbone', or 'heavy'.

        Returns:
            (model coordinates, reference coordinates) in Å, row by row the same atom.

        Raises:
            ValueError: If atoms is unknown.
        """
        if atoms not in SUPERPOSITION_ATOMS:
            raise ValueError(f"atoms must be one of {tuple(SUPERPOSITION_ATOMS)}.")
        model_xyz, ref_xyz = [], []
        for key, seq, ref_key, ref_seq in pairs:
            record = self.chains[key]
            residue, ref_residue = record.residues[seq], reference.chains[ref_key].residues[ref_seq]
            for name in _atom_names(residue, record.moieties[seq], atoms):
                if name in ref_residue:
                    model_xyz.append(residue[name].coord)
                    ref_xyz.append(ref_residue[name].coord)
        return np.array(model_xyz, float).reshape(-1, 3), np.array(ref_xyz, float).reshape(-1, 3)

    def residue_match_table(self, reference: SyncedStructure, chain_map: Mapping[str, str] | None = None
                            ) -> pd.DataFrame:
        """Return the residues of this structure matched to residues of a reference structure.

        Residues match through their reference positions; when the reference sequences of two chains differ (tags,
        truncations, mutations), through a global alignment of them.

        Args:
            reference: Reference structure.
            chain_map: {model chain: reference chain}, e.g. {'A': 'H'}, or None to pair chains with the same key.

        Returns:
            One row per matched residue with RESIDUE_MATCH_COLUMNS; sequence indices are 1-based.
        """
        rows = []
        for key, seq, ref_key, ref_seq in self._residue_pairs(reference, chain_map):
            record, ref_record = self.chains[key], reference.chains[ref_key]
            residue, ref_residue = record.residues[seq], ref_record.residues[ref_seq]
            rows.append({
                "structure_id": self.structure_id, "chain_key": key, "chain_id": record.chain_id,
                "auth_chain_id": record.auth_chain_id, "pdb_resnum": residue.id[1],
                "resname": residue.resname.strip().upper(), "sequence_index": seq + 1,
                "reference_structure_id": reference.structure_id, "reference_chain_key": ref_key,
                "reference_chain_id": ref_record.chain_id, "reference_auth_chain_id": ref_record.auth_chain_id,
                "reference_pdb_resnum": ref_residue.id[1], "reference_resname": ref_residue.resname.strip().upper(),
                "reference_sequence_index": ref_seq + 1,
            })
        return pd.DataFrame(rows, columns=RESIDUE_MATCH_COLUMNS)

    def rmsd_metrics(self, reference: SyncedStructure, chains: str | ChainRefs | None = None,
                     chain_map: Mapping[str, str] | None = None, atoms: str = "ca",
                     fit_chains: str | ChainRefs | None = None, superpose: bool = True) -> Record:
        """Return the RMSD of this structure to a reference, after an optimal (Kabsch) superposition.

        Neither structure moves; apply_superposition() moves this one. fit_chains gives RMSDs such as the ligand RMSD
        of docking: fit on the receptor, measure the ligand.

        Args:
            reference: Reference structure.
            chains: Model chains to measure, or None for every chain with a counterpart.
            chain_map: {model chain: reference chain}, or None to pair chains with the same key.
            atoms: Atoms compared: 'ca' (C-alpha, P for nucleic acids), 'backbone', or 'heavy' (shared heavy atoms).
            fit_chains: Model chains to superpose on, or None for the measured chains.
            superpose: Whether to superpose first; False measures the coordinates as they are.

        Returns:
            {'rmsd_A', 'n_residues', 'n_atoms', 'n_fit_atoms', 'rotation', 'translation'}; rotation and translation
            move this structure onto the reference as coord @ rotation + translation.
        """
        pairs = self._residue_pairs(reference, chain_map, chains)
        model_xyz, ref_xyz = self._paired_coords(reference, pairs, atoms)
        rotation, translation = np.eye(3), np.zeros(3)
        n_fit = 0
        if superpose:
            fit_pairs = pairs if fit_chains is None else self._residue_pairs(reference, chain_map, fit_chains)
            fit_model, fit_ref = self._paired_coords(reference, fit_pairs, atoms)
            rotation, translation = _kabsch(fit_model, fit_ref)
            n_fit = len(fit_model)
        return {
            "rmsd_A": round(_rmsd(model_xyz @ rotation + translation, ref_xyz), 4),
            "n_residues": len(pairs), "n_atoms": len(model_xyz), "n_fit_atoms": n_fit,
            "rotation": rotation, "translation": translation,
        }

    def apply_superposition(self, reference: SyncedStructure, chains: str | ChainRefs | None = None,
                            chain_map: Mapping[str, str] | None = None, atoms: str = "ca") -> Self:
        """Move this structure onto a reference by an optimal (Kabsch) superposition of matched atoms.

        Every atom of every model moves, ligands and water included; cached interface metrics are cleared.

        Args:
            reference: Reference structure.
            chains: Model chains to superpose on, or None for every chain with a counterpart.
            chain_map: {model chain: reference chain}, or None to pair chains with the same key.
            atoms: Atoms fitted: 'ca', 'backbone', or 'heavy'.

        Returns:
            This structure.
        """
        metrics = self.rmsd_metrics(reference, chains, chain_map, atoms)
        rotation, translation = metrics["rotation"], metrics["translation"]
        for atom in self.structure.get_atoms():
            for child in atom.disordered_get_list() if atom.is_disordered() else [atom]:
                child.transform(rotation, translation)
        self._interface_cache.clear()
        log_info(f"{self.structure_id}: superposed onto {reference.structure_id} over {metrics['n_atoms']} atoms, "
                 f"RMSD {metrics['rmsd_A']:.2f} Å")
        return self

    def _interface_contacts(self, receptor: Sequence[str], ligand: Sequence[str], cutoff: float
                            ) -> tuple[set[tuple[str, int, str, int]], set[tuple[str, int]]]:
        """Return the receptor-ligand residue contacts and the residues near the other partner.

        Distances are the minimum heavy-atom distance between two residues.

        Args:
            receptor: Receptor chain keys.
            ligand: Ligand chain keys.
            cutoff: Distance in Å below which two residues are in contact.

        Returns:
            ({(receptor key, index, ligand key, index)}, {(chain key, index)} of residues within CAPRI_INTERFACE_A
            of the other partner); indices are 0-based.
        """
        block = int(self.config["interface_block_size"])
        contacts: set[tuple[str, int, str, int]] = set()
        interface: set[tuple[str, int]] = set()
        coords = {key: InterfaceMixin._heavy_atom_coords(self.chains[key].residues) for key in [*receptor, *ligand]}
        for key_a in receptor:
            xa, oa = coords[key_a]
            for key_b in ligand:
                xb, ob = coords[key_b]
                if not len(xa) or not len(xb):
                    continue
                distances = InterfaceMixin._residue_min_distances(
                    xa, oa, xb, ob, len(self.chains[key_a].residues), len(self.chains[key_b].residues), block)
                contacts.update((key_a, int(i), key_b, int(j)) for i, j in zip(*np.nonzero(distances < cutoff)))
                near = distances < CAPRI_INTERFACE_A
                interface.update((key_a, int(i)) for i in np.nonzero(near.any(axis=1))[0])
                interface.update((key_b, int(j)) for j in np.nonzero(near.any(axis=0))[0])
        return contacts, interface

    def capri_metrics(self, reference: SyncedStructure, receptor: str | ChainRefs, ligand: str | ChainRefs,
                      chain_map: Mapping[str, str] | None = None) -> Record:
        """Return the CAPRI quality metrics of this structure, a model of a complex, against the reference complex.

        Definitions follow CAPRI and DockQ: native contacts are receptor-ligand residue pairs of the reference within
        CAPRI_CONTACT_A (5 Å, heavy atoms); interface residues are reference residues within CAPRI_INTERFACE_A (10 Å)
        of the other partner. i-RMSD is the backbone RMSD of the interface residues after superposing them; l-RMSD
        the ligand backbone RMSD after superposing the receptor backbone; rmsd_A the backbone RMSD of both partners.
        Backbone atoms are N, CA, C, O for proteins and P, O5', C5', C4', C3', O3' for nucleic acids.

        Args:
            reference: Reference (native) complex.
            receptor: Model chains of the receptor; a key, label ID, or author ID, or several.
            ligand: Model chains of the ligand.
            chain_map: {model chain: reference chain}, or None to pair chains with the same key.

        Returns:
            The CAPRI_COLUMNS values: rmsd_A, irmsd_A, lrmsd_A, fnat, fnonnat, dockq, capri_class, and counts.

        Raises:
            ValueError: If the receptor and ligand share a chain or the reference has no receptor-ligand contacts.
        """
        receptor_keys = [self._chain_key(chain) for chain in _chain_list(receptor)]
        ligand_keys = [self._chain_key(chain) for chain in _chain_list(ligand)]
        if set(receptor_keys) & set(ligand_keys):
            raise ValueError("receptor and ligand must not share a chain.")
        chain_pairs = dict(self._chain_pairs(reference, chain_map, [*receptor_keys, *ligand_keys]))
        ref_receptor = [chain_pairs[key] for key in receptor_keys]
        ref_ligand = [chain_pairs[key] for key in ligand_keys]
        pairs = self._residue_pairs(reference, chain_map, [*receptor_keys, *ligand_keys])
        to_ref = {(key, seq): (ref_key, ref_seq) for key, seq, ref_key, ref_seq in pairs}

        native, ref_interface = reference._interface_contacts(ref_receptor, ref_ligand, CAPRI_CONTACT_A)
        if not native:
            raise ValueError(f"{reference.structure_id} has no contact between receptor {ref_receptor} and ligand "
                             f"{ref_ligand} within {CAPRI_CONTACT_A} Å.")
        model_contacts, _ = self._interface_contacts(receptor_keys, ligand_keys, CAPRI_CONTACT_A)
        mapped = {(*to_ref.get((ka, i), (None, -1)), *to_ref.get((kb, j), (None, -1)))
                  for ka, i, kb, j in model_contacts}
        fnat = len(native & mapped) / len(native)
        fnonnat = len(mapped - native) / len(mapped) if mapped else 0.0

        interface_pairs = [pair for pair in pairs if (pair[2], pair[3]) in ref_interface]
        irmsd = self._fitted_rmsd(reference, interface_pairs, interface_pairs)
        receptor_pairs = [pair for pair in pairs if pair[0] in receptor_keys]
        ligand_pairs = [pair for pair in pairs if pair[0] in ligand_keys]
        lrmsd = self._fitted_rmsd(reference, receptor_pairs, ligand_pairs)
        rmsd = self._fitted_rmsd(reference, pairs, pairs)
        return {
            "structure_id": self.structure_id, "reference_structure_id": reference.structure_id,
            "receptor": ",".join(receptor_keys), "ligand": ",".join(ligand_keys),
            "rmsd_A": round(rmsd, 3), "irmsd_A": round(irmsd, 3), "lrmsd_A": round(lrmsd, 3),
            "fnat": round(fnat, 4), "fnonnat": round(fnonnat, 4), "dockq": round(dockq_score(fnat, irmsd, lrmsd), 4),
            "capri_class": capri_class(fnat, irmsd, lrmsd), "n_native_contacts": len(native),
            "n_model_contacts": len(model_contacts), "n_interface_residues": len(ref_interface),
            "n_matched_residues": len(pairs),
        }

    def _fitted_rmsd(self, reference: SyncedStructure, fit_pairs: Sequence[ResiduePair],
                     measure_pairs: Sequence[ResiduePair], atoms: str = "backbone") -> float:
        """Return the RMSD of some matched residues after superposing others.

        Args:
            reference: Reference structure.
            fit_pairs: Matched residues to superpose on.
            measure_pairs: Matched residues to measure.
            atoms: Atom selection.

        Returns:
            The RMSD in Å.
        """
        fit_model, fit_ref = self._paired_coords(reference, fit_pairs, atoms)
        rotation, translation = _kabsch(fit_model, fit_ref)
        model_xyz, ref_xyz = self._paired_coords(reference, measure_pairs, atoms)
        return _rmsd(model_xyz @ rotation + translation, ref_xyz)

    def capri_table(self, models: Mapping[str, SyncedStructure | PathLike] | Iterable[SyncedStructure | PathLike],
                    receptor: str | ChainRefs, ligand: str | ChainRefs,
                    chain_map: Mapping[str, str] | None = None) -> pd.DataFrame:
        """Return the CAPRI quality metrics of several models against this structure, the reference complex.

        Args:
            models: Models as structures or structure file paths, or {name: model}; names replace structure IDs.
            receptor: Model chains of the receptor; a key, label ID, or author ID, or several.
            ligand: Model chains of the ligand.
            chain_map: {model chain: reference chain}, or None to pair chains with the same key.

        Returns:
            One row per model with CAPRI_COLUMNS, best DockQ first.
        """
        from .structure import SyncedStructure

        reference = cast("SyncedStructure", self)
        named = models.items() if isinstance(models, Mapping) else ((None, model) for model in models)
        rows = []
        for name, model in named:
            structure = model if isinstance(model, SyncedStructure) else SyncedStructure.from_path(Path(model))
            row = structure.capri_metrics(reference, receptor, ligand, chain_map)
            if name is not None:
                row["structure_id"] = str(name)
            rows.append(row)
        table = pd.DataFrame(rows, columns=CAPRI_COLUMNS)
        return table.sort_values("dockq", ascending=False, kind="stable").reset_index(drop=True)
