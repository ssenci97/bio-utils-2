"""C-alpha and sidechain contact detection, contact strings, and contact summaries with standard and author chain IDs."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any

import mdtraj as md
import numpy as np
import pandas as pd

from ._types import ChainRefs, Moieties, PathLike, Record
from .base import StructureState
from .constants import (
    ANGSTROM_PER_NM, BACKBONE_ATOMS, CONTACT_COLUMNS, CONTACT_SCHEMES, CONTACT_SEQUENCE_COLUMNS, CONTACT_TYPES,
    CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES, DEFAULT_CONTACT_TYPE, EXPORT_SEP, MD_SCHEMES, MOIETIES,
    MOIETY_NUCLEIC, MOIETY_PROTEIN, NUCLEIC_PHOSPHATE_ATOMS,
)
from .residues import is_heavy, representative_coord

if TYPE_CHECKING:
    from collections.abc import Callable

    from Bio.PDB.Atom import Atom
    from Bio.PDB.Residue import Residue

###################################### CONTACT ARGUMENTS


def _sidechain_atoms(residue: Residue, moiety: str) -> list[Atom]:
    """Return the heavy sidechain atoms of a protein residue or the heavy base atoms of a nucleotide.

    Args:
        residue: Biopython residue.
        moiety: Biomolecule type of the residue.

    Returns:
        The atoms; empty for other moieties and for residues without such atoms, such as glycine.
    """
    backbone: Callable[[str], bool]
    if moiety == MOIETY_PROTEIN:
        backbone = lambda name: name in BACKBONE_ATOMS
    elif moiety == MOIETY_NUCLEIC:
        backbone = lambda name: name in NUCLEIC_PHOSPHATE_ATOMS or "'" in name or "*" in name
    else:
        return []
    return [atom for atom in residue.get_atoms() if not backbone(atom.get_id().strip()) and is_heavy(atom)]


def _check_contact_args(cutoff: float, contact_type: str, scheme: str, moieties: Moieties) -> None:
    """Validate contact arguments.

    Args:
        cutoff: Distance cutoff in nm.
        contact_type: 'inter', 'intra', or 'all'.
        scheme: Contact scheme.
        moieties: Biomolecule types.

    Raises:
        ValueError: If any argument is invalid.
    """
    if cutoff <= 0:
        raise ValueError("cutoff must be greater than zero.")
    if contact_type not in CONTACT_TYPES:
        raise ValueError(f"contact_type must be one of {CONTACT_TYPES}.")
    if scheme not in CONTACT_SCHEMES:
        raise ValueError(f"scheme must be one of {CONTACT_SCHEMES}.")
    if not set(moieties) <= set(MOIETIES):
        raise ValueError(f"moieties must be a subset of {MOIETIES}.")


def _pair_mask(same_chain: np.ndarray, contact_type: str) -> np.ndarray:
    """Select residue pairs of the requested contact type.

    Args:
        same_chain: Boolean array, True where both residues of a pair are in the same chain.
        contact_type: 'inter', 'intra', or 'all'.

    Returns:
        Boolean mask over the pairs.
    """
    if contact_type == "all":
        return np.ones_like(same_chain, dtype=bool)
    return same_chain if contact_type == "intra" else ~same_chain

###################################### CONTACT POINTS AND TRAJECTORY


class ContactMixin(StructureState):
    """Contact points, contact detection with MDTraj, and contact strings."""

    def _contact_records(self, chain_ids: ChainRefs | None = None,
                         moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> list[Record]:
        """Collect one record per residue of the given moieties that has a contact point.

        Args:
            chain_ids: Chain keys to include, or None for all.
            moieties: Biomolecule types to include, from MOIETIES.

        Returns:
            Records with the chain key, both chain IDs, 0-based sequence index (None for hetero residues),
            residue number and name, moiety, residue, and contact point in A.
        """
        selected = None if chain_ids is None else {str(c) for c in chain_ids}
        records = []
        for key, data in self.chains.items():
            if selected is not None and key not in selected:
                continue
            residues: list[tuple[int | None, Residue, str]] = [
                (index, residue, moiety) for index, (residue, moiety) in enumerate(zip(data.residues, data.moieties))]
            residues += [(None, residue, moiety) for residue, moiety in data.hetero]
            for index, residue, moiety in residues:
                if moiety not in moieties:
                    continue
                point = representative_coord(residue, moiety)
                if point is not None:
                    records.append({"key": key, "chain_id": data.chain_id, "auth_chain_id": data.auth_chain_id,
                                    "sequence_index": index, "resnum": int(residue.id[1]),
                                    "resname": residue.resname.strip().upper(), "moiety": moiety,
                                    "residue": residue, "point": point})
        return records

    @staticmethod
    def _make_trajectory(records: list[Record], scheme: str = "ca") -> md.Trajectory:
        """Build a one-frame MDTraj trajectory in nm with one MDTraj residue per record.

        With 'ca', each residue is one atom named CA at its contact point. With 'sidechain', protein residues
        contribute their heavy sidechain atoms and nucleotides their heavy base atoms; residues without such atoms
        (such as Gly) and other moieties contribute one CP atom at their contact point.

        Args:
            records: Records from _contact_records.
            scheme: 'ca' or 'sidechain'.

        Returns:
            The trajectory.

        Raises:
            ValueError: If there are no records.
        """
        if not records:
            raise ValueError("No contact points found in the selected chains for the requested moieties.")
        topology = md.Topology()
        md_chains: dict[str, Any] = {}
        coordinates = []
        for record in records:
            if record["key"] not in md_chains:
                md_chains[record["key"]] = topology.add_chain(chain_id=record["key"])
            md_residue = topology.add_residue(record["resname"], md_chains[record["key"]], resSeq=record["resnum"])
            residue = record["residue"]
            atoms = _sidechain_atoms(residue, record["moiety"]) if scheme == "sidechain" else []
            if atoms:
                for atom in atoms:
                    letters = "".join(char for char in atom.get_id() if char.isalpha())
                    symbol = (atom.element or letters[:1] or "C").strip().capitalize()
                    topology.add_atom(atom.get_id().strip(), md.element.get_by_symbol(symbol), md_residue)
                    coordinates.append(atom.coord)
            else:
                topology.add_atom("CA" if scheme == "ca" else "CP", md.element.carbon, md_residue)
                coordinates.append(record["point"])
        xyz = np.asarray(coordinates, dtype=np.float32)[None, :, :] / ANGSTROM_PER_NM
        return md.Trajectory(xyz=xyz, topology=topology)

    def contact_point_trajectory(self, chain_ids: ChainRefs | None = None,
                                 moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> md.Trajectory:
        """Return the contact-point MDTraj trajectory in nm: Cα, C4′, or centroid per residue.

        Args:
            chain_ids: Chains to include, or None for all.
            moieties: Biomolecule types to include.

        Returns:
            The trajectory.
        """
        if chain_ids is not None:
            chain_ids = [self._require_chain(chain_id) for chain_id in chain_ids]
        return self._make_trajectory(self._contact_records(chain_ids, moieties), "ca")

###################################### CONTACTS

    def compute_contacts(self, cutoff: float = CONTACTS_CUTOFF_NM, contact_type: str = DEFAULT_CONTACT_TYPE,
                         chain_ids: ChainRefs | None = None, scheme: str = "ca",
                         moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> pd.DataFrame:
        """Find residue pairs closer than the cutoff.

        Args:
            cutoff: Distance cutoff in nm.
            contact_type: 'inter', 'intra', or 'all'.
            chain_ids: Chains to include, or None for all.
            scheme: 'ca' (one contact point per residue) or 'sidechain' (protein sidechain and nucleotide base
                atoms, plus one contact point per other residue).
            moieties: Biomolecule types to include, from MOIETIES. Hetero residues have no sequence index.

        Returns:
            One row per contact, with 1-based sequence indices.
        """
        _check_contact_args(cutoff, contact_type, scheme, moieties)
        if chain_ids is not None:
            chain_ids = [self._require_chain(chain_id) for chain_id in chain_ids]
        records = self._contact_records(chain_ids, tuple(moieties))
        trajectory = self._make_trajectory(records, scheme)
        keys = np.array([record["key"] for record in records])
        first, second = np.triu_indices(len(records), k=1)
        keep = _pair_mask(keys[first] == keys[second], contact_type)
        pairs = np.column_stack([first[keep], second[keep]])
        if not len(pairs):
            return pd.DataFrame(columns=CONTACT_COLUMNS)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, message=".*No alpha carbon.*")
            distances, residue_pairs = md.compute_contacts(trajectory, contacts=pairs, scheme=MD_SCHEMES[scheme],
                                                           ignore_nonprotein=False, periodic=False)
        rows = []
        for (i, j), distance in zip(residue_pairs, distances[0]):
            if distance > cutoff:
                continue
            a, b = records[int(i)], records[int(j)]
            rows.append({
                "chain_1": a["chain_id"], "auth_chain_1": a["auth_chain_id"], "resnum_1": a["resnum"],
                "resname_1": a["resname"],
                "sequence_index_1": None if a["sequence_index"] is None else a["sequence_index"] + 1,
                "chain_2": b["chain_id"], "auth_chain_2": b["auth_chain_id"], "resnum_2": b["resnum"],
                "resname_2": b["resname"],
                "sequence_index_2": None if b["sequence_index"] is None else b["sequence_index"] + 1,
                "min_distance_nm": round(float(distance), 4),
                "contact_type": "intra" if a["key"] == b["key"] else "inter",
                "moiety_1": a["moiety"],
                "moiety_2": b["moiety"],
            })
        frame = pd.DataFrame(rows, columns=CONTACT_COLUMNS)
        return frame.astype({"sequence_index_1": "Int64", "sequence_index_2": "Int64"})

###################################### CONTACT STRINGS

    def _string_indices(self, chain_id: str, moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> list[int]:
        """List the polymer residues of a chain whose moiety is in moieties.

        Args:
            chain_id: Chain key.
            moieties: Biomolecule types to keep.

        Returns:
            Their 0-based sequence indices.
        """
        return [index for index, moiety in enumerate(self.chains[chain_id].moieties) if moiety in moieties]

    def _contact_bits(self, chain_id: str, cutoff: float, contact_type: str, contacts: pd.DataFrame | None,
                      scheme: str, moieties: Moieties) -> list[str]:
        """Mark each polymer residue of a chain as X (no contact point), 0, or 1.

        Args:
            chain_id: Chain key.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts table, or None.
            scheme: Contact scheme, used if contacts is None.
            moieties: Biomolecule types, used if contacts is None.

        Returns:
            One character per polymer residue.
        """
        label = self.chains[chain_id].chain_id
        if contacts is None:
            contacts = self.compute_contacts(cutoff=cutoff, contact_type=contact_type, scheme=scheme, moieties=moieties)
        bits = ["X"] * self.chains[chain_id].sequence_length
        for record in self._contact_records([chain_id], moieties):
            if record["sequence_index"] is not None:
                bits[record["sequence_index"]] = "0"
        for side in ("1", "2"):
            indices = contacts.loc[contacts[f"chain_{side}"] == label, f"sequence_index_{side}"].dropna()
            for index in indices.astype(int):
                bits[index - 1] = "1"
        return bits

    def contact_string(self, chain_id: str, cutoff: float = CONTACTS_CUTOFF_NM,
                       contact_type: str = DEFAULT_CONTACT_TYPE, contacts: pd.DataFrame | None = None,
                       scheme: str = "ca", moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> str:
        """Return a per-residue string of X (no contact point), 0 (no contact), or 1 (contact).

        Only polymer residues whose moiety is in moieties appear in the string.

        Args:
            chain_id: Chain identifier.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts table, or None.
            scheme: Contact scheme, used if contacts is None.
            moieties: Biomolecule types included in the string.

        Returns:
            The contact string along the observed sequence.
        """
        chain_id = self._require_chain(chain_id)
        bits = self._contact_bits(chain_id, cutoff, contact_type, contacts, scheme, moieties)
        return "".join(bits[index] for index in self._string_indices(chain_id, moieties))

    def reference_contact_string(self, chain_id: str, cutoff: float = CONTACTS_CUTOFF_NM,
                                 contact_type: str = DEFAULT_CONTACT_TYPE, contacts: pd.DataFrame | None = None,
                                 scheme: str = "ca", moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> str:
        """Return the contact string laid out along the full reference sequence.

        Reference positions not mapped to an included observed residue are 'X'; mapped positions are '1' in
        contact and '0' otherwise, including residues without a contact point.

        Args:
            chain_id: Chain identifier.
            cutoff: Distance cutoff in nm, used if contacts is None.
            contact_type: Contact type, used if contacts is None.
            contacts: Precomputed contacts table, or None.
            scheme: Contact scheme, used if contacts is None.
            moieties: Biomolecule types included in the string.

        Returns:
            The contact string along the reference sequence.

        Raises:
            ValueError: If the chain has no reference sequence.
        """
        chain_id = self._require_chain(chain_id)
        data = self.chains[chain_id]
        reference = data.reference_sequence
        if not reference:
            raise ValueError(f"Chain {data.chain_id} has no reference sequence. Set one with set_reference_sequence() first.")
        observed = self._contact_bits(chain_id, cutoff, contact_type, contacts, scheme, moieties)
        kept = set(self._string_indices(chain_id, moieties))
        bits = ["X"] * len(reference)
        for seq_index, ref_index in data.sequence_index_to_ref_index.items():
            if seq_index in kept and 0 <= ref_index < len(bits):
                bits[ref_index] = "1" if observed[seq_index] == "1" else "0"
        return "".join(bits)

    def export_contact_strings(self, output_path: PathLike, cutoff: float = CONTACTS_CUTOFF_NM,
                               contact_type: str = DEFAULT_CONTACT_TYPE, reference: bool = False,
                               scheme: str = "ca", moieties: Moieties = DEFAULT_CONTACT_MOIETIES) -> Path:
        """Write each chain's sequence and contact string to a TSV file.

        Args:
            output_path: Destination TSV path.
            cutoff: Distance cutoff in nm.
            contact_type: Contact type.
            reference: Whether to use reference sequences, skipping chains without one.
            scheme: Contact scheme ('ca' or 'sidechain').
            moieties: Biomolecule types included in the output.

        Returns:
            The written path.
        """
        contacts = self.compute_contacts(cutoff=cutoff, contact_type=contact_type, scheme=scheme, moieties=moieties)
        options: dict[str, Any] = {"cutoff": cutoff, "contact_type": contact_type, "contacts": contacts,
                                   "scheme": scheme, "moieties": moieties}
        rows = []
        for chain_key, data in self.chains.items():
            if reference:
                if not data.reference_sequence:
                    continue
                sequence = data.reference_sequence
                text = self.reference_contact_string(chain_key, **options)
            else:
                sequence = "".join(data.sequence[index] for index in self._string_indices(chain_key, moieties))
                text = self.contact_string(chain_key, **options)
            rows.append({
                "structure_id": self.structure_id,
                "chain_id": data.chain_id,
                "auth_chain_id": data.auth_chain_id,
                "sequence": sequence,
                "contacts": text,
                "contact_type": contact_type,
                "mdtraj_scheme": scheme,
                "cutoff_nm": cutoff,
            })
        path = self._prepare_path(output_path)
        pd.DataFrame(rows, columns=CONTACT_SEQUENCE_COLUMNS).to_csv(path, sep=EXPORT_SEP, index=False)
        return path
