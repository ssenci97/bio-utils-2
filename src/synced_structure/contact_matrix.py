"""Residue-pair contact matrices for one chain, a chain pair, or a block of several chains."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

from collections.abc import Sequence

import numpy as np
import pandas as pd

from ._types import Moieties, Record
from .base import StructureState
from .constants import COMPLETE_MIN_SEPARATION, CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES

INDEX_COLUMNS = ["sequence_index_1", "sequence_index_2"]

###################################### CONTACT MATRIX


class ContactMatrixMixin(StructureState):
    """Residue-pair and block contact matrices."""

    def _resolve_pair(self, chain_a: str, chain_b: str | None) -> tuple[str, str, str, str]:
        """Resolve two chains, with chain_b defaulting to chain_a for intra-chain maps.

        Args:
            chain_a: First chain identifier, standard or author.
            chain_b: Second chain identifier, or None for an intra-chain map.

        Returns:
            (key_a, key_b, label_a, label_b).
        """
        key_a = self._require_chain(chain_a)
        key_b = key_a if chain_b is None else self._require_chain(chain_b)
        return key_a, key_b, self.chains[key_a].chain_id, self.chains[key_b].chain_id

    @staticmethod
    def _pair_rows(contacts: pd.DataFrame, label_a: str, label_b: str) -> pd.DataFrame:
        """Select the contact rows that join two chains, in either order.

        Args:
            contacts: Contacts table.
            label_a: Standard ID of the first chain.
            label_b: Standard ID of the second chain.

        Returns:
            The matching rows.
        """
        forward = (contacts["chain_1"] == label_a) & (contacts["chain_2"] == label_b)
        backward = (contacts["chain_1"] == label_b) & (contacts["chain_2"] == label_a)
        return contacts[forward | backward]

    @staticmethod
    def _multi_block_entries(contacts: pd.DataFrame, labels: Sequence[str],
                             offsets: Sequence[int] | np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Place polymer contacts among the labeled chains in a block matrix; hetero contacts are dropped.

        Chain k occupies indices offsets[k] onward, in the order of labels.

        Args:
            contacts: Contacts table covering all chains in labels.
            labels: Standard IDs of the chains, in block order.
            offsets: Start index of each chain in the block matrix.

        Returns:
            0-based rows, columns, and distances in nm.
        """
        position = {label: k for k, label in enumerate(labels)}
        contacts = contacts.dropna(subset=INDEX_COLUMNS)
        contacts = contacts[contacts["chain_1"].isin(position) & contacts["chain_2"].isin(position)]
        starts = np.asarray(offsets, dtype=int)
        block_rows = starts[contacts["chain_1"].map(position).to_numpy(dtype=int)] + contacts["sequence_index_1"].to_numpy(dtype=int) - 1
        block_cols = starts[contacts["chain_2"].map(position).to_numpy(dtype=int)] + contacts["sequence_index_2"].to_numpy(dtype=int) - 1
        return block_rows, block_cols, contacts["min_distance_nm"].to_numpy(dtype=float)

    def _block_matrix(self, keys: Sequence[str], cutoff: float, scheme: str, moieties: Moieties) -> Record:
        """Build a symmetric block matrix of minimum distances for several chains.

        Args:
            keys: Chain keys in block order.
            cutoff: Distance cutoff in nm.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            moieties: Biomolecule types to include.

        Returns:
            Dict with 'matrix' (nm, NaN for non-contacts), 'sizes', 'offsets', 'labels', 'keys', and 'contacts'.
        """
        labels = [self.chains[k].chain_id for k in keys]
        sizes = [self.chains[k].sequence_length for k in keys]
        offsets = np.cumsum([0] + sizes[:-1]).astype(int)
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all", chain_ids=list(dict.fromkeys(keys)),
                                         scheme=scheme, moieties=moieties)
        matrix = np.full((sum(sizes), sum(sizes)), np.nan)
        if not contacts.empty:
            gi, gj, distance = self._multi_block_entries(contacts, labels, offsets)
            matrix[gi, gj] = distance
            matrix[gj, gi] = distance
        return {"matrix": matrix, "sizes": sizes, "offsets": offsets, "labels": labels, "keys": list(keys),
                "contacts": contacts}

    def _matrix_from_contacts(self, contacts: pd.DataFrame, chain_a: str, chain_b: str) -> np.ndarray:
        """Build the chain_a x chain_b minimum-distance matrix, symmetric if intra-chain.

        Args:
            contacts: Contacts table from compute_contacts, covering both chains.
            chain_a: Chain key on the rows.
            chain_b: Chain key on the columns.

        Returns:
            Distances in nm, NaN for non-contacts.
        """
        label_a, label_b = self.chains[chain_a].chain_id, self.chains[chain_b].chain_id
        matrix = np.full((self.chains[chain_a].sequence_length, self.chains[chain_b].sequence_length), np.nan)
        rows = self._pair_rows(contacts, label_a, label_b).dropna(subset=INDEX_COLUMNS)
        forward = (rows["chain_1"] == label_a).to_numpy()
        first = rows["sequence_index_1"].to_numpy(dtype=int) - 1
        second = rows["sequence_index_2"].to_numpy(dtype=int) - 1
        i, j = np.where(forward, first, second), np.where(forward, second, first)
        distance = rows["min_distance_nm"].to_numpy(dtype=float)
        matrix[i, j] = distance
        if chain_a == chain_b:
            matrix[j, i] = distance
        return matrix

    def contact_matrix(self, chain_a: str = "A", chain_b: str | None = None, cutoff: float = CONTACTS_CUTOFF_NM,
                       moieties: Moieties = DEFAULT_CONTACT_MOIETIES, scheme: str = "ca") -> Record:
        """Return the residue-pair minimum-distance matrix of one chain or a chain pair.

        Args:
            chain_a: Chain on the rows, standard or author ID.
            chain_b: Chain on the columns, or None for an intra-chain map.
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.

        Returns:
            Dict with 'matrix' (nm, NaN where not in contact), 'labels' (standard IDs), 'keys', and 'contacts'.
        """
        key_a, key_b, label_a, label_b = self._resolve_pair(chain_a, chain_b)
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all", chain_ids=sorted({key_a, key_b}),
                                         scheme=scheme, moieties=moieties)
        return {"matrix": self._matrix_from_contacts(contacts, key_a, key_b), "labels": (label_a, label_b),
                "keys": (key_a, key_b), "contacts": contacts}

    def block_contact_matrix(self, chain_a: str = "A", chain_b: str = "B", cutoff: float = CONTACTS_CUTOFF_NM,
                             moieties: Moieties = DEFAULT_CONTACT_MOIETIES, scheme: str = "ca",
                             min_separation: int = COMPLETE_MIN_SEPARATION) -> Record:
        """Return the block contact matrix of two chains with intra- and inter-chain masks.

        Args:
            chain_a: First chain, standard or author ID.
            chain_b: Second chain, standard or author ID.
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            min_separation: Smallest sequence separation kept in the intra-chain mask.

        Returns:
            Dict with 'matrix', 'intra', 'inter' (boolean masks), 'sizes', 'labels', 'keys', and 'contacts'.

        Raises:
            ValueError: If both arguments name the same chain.
        """
        key_a, key_b, _, _ = self._resolve_pair(chain_a, chain_b)
        if key_a == key_b:
            raise ValueError("chain_a and chain_b must be different chains for a block map.")
        block = self._block_matrix([key_a, key_b], cutoff, scheme, moieties)
        n1 = block["sizes"][0]
        index = np.arange(len(block["matrix"]))
        same_block = (index < n1)[:, None] == (index < n1)[None, :]
        separation = np.abs(index[:, None] - index[None, :])
        contact = ~np.isnan(block["matrix"])
        return {
            "matrix": block["matrix"],
            "intra": contact & same_block & (separation >= min_separation),
            "inter": contact & ~same_block,
            "sizes": tuple(block["sizes"]),
            "labels": tuple(block["labels"]),
            "keys": (key_a, key_b),
            "contacts": block["contacts"],
        }
