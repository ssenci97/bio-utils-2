"""Residue-pair contact matrices for one chain, a chain pair, or a block of several chains."""

###################################### IMPORTS AND CONSTANTS

import numpy as np

from .constants import COMPLETE_MIN_SEPARATION, CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES

INDEX_COLUMNS = ["sequence_index_1", "sequence_index_2"]

###################################### CONTACT MATRIX


class ContactMatrixMixin:
    def _pair_keys(self, chain_a, chain_b):
        """Return (key_a, key_b, label_a, label_b), with chain_b defaulting to chain_a for intra-chain maps.

        Args:
            chain_a: First chain identifier, standard or author.
            chain_b: Second chain identifier, or None for an intra-chain map.
        """
        chain_a = self._require_chain(chain_a)
        chain_b = chain_a if chain_b is None else self._require_chain(chain_b)
        return chain_a, chain_b, self.chains[chain_a]["chain_id"], self.chains[chain_b]["chain_id"]

    @staticmethod
    def _pair_rows(contacts, label_a, label_b):
        """Return the contact rows that join two chains, in either order."""
        forward = (contacts["chain_1"] == label_a) & (contacts["chain_2"] == label_b)
        backward = (contacts["chain_1"] == label_b) & (contacts["chain_2"] == label_a)
        return contacts[forward | backward]

    @staticmethod
    def _multi_block_entries(contacts, labels, offsets):
        """Return 0-based block-matrix rows, columns, and distances for polymer contacts among the labeled chains.

        Chain k occupies indices offsets[k] onward, in the order of labels. Hetero contacts are dropped.

        Args:
            contacts: Contacts DataFrame covering all chains in labels.
            labels: Standard IDs of the chains, in block order.
            offsets: Start index of each chain in the block matrix.
        """
        position = {label: k for k, label in enumerate(labels)}
        contacts = contacts.dropna(subset=INDEX_COLUMNS)
        contacts = contacts[contacts["chain_1"].isin(position) & contacts["chain_2"].isin(position)]
        offsets = np.asarray(offsets, dtype=int)
        gi = offsets[contacts["chain_1"].map(position).to_numpy(dtype=int)] + contacts["sequence_index_1"].to_numpy(dtype=int) - 1
        gj = offsets[contacts["chain_2"].map(position).to_numpy(dtype=int)] + contacts["sequence_index_2"].to_numpy(dtype=int) - 1
        return gi, gj, contacts["min_distance_nm"].to_numpy(dtype=float)

    def _block_matrix(self, keys, cutoff, scheme, moieties):
        """Return a symmetric block matrix of minimum distances in nm for several chains, NaN for non-contacts.

        Args:
            keys: Chain keys in block order.
            cutoff: Distance cutoff in nm.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            moieties: Biomolecule types to include.

        Returns:
            Dict with 'matrix', 'sizes', 'offsets', 'labels', 'keys', and 'contacts'.
        """
        labels = [self.chains[k]["chain_id"] for k in keys]
        sizes = [self.chains[k]["sequence_length"] for k in keys]
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

    def _matrix_from_contacts(self, contacts, chain_a, chain_b):
        """Return the chain_a x chain_b minimum-distance matrix in nm, NaN for non-contacts; symmetric if intra-chain.

        Args:
            contacts: Contacts DataFrame from compute_contacts, covering both chains.
            chain_a: Chain key on the rows.
            chain_b: Chain key on the columns.
        """
        label_a, label_b = self.chains[chain_a]["chain_id"], self.chains[chain_b]["chain_id"]
        matrix = np.full((self.chains[chain_a]["sequence_length"], self.chains[chain_b]["sequence_length"]), np.nan)
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

    def contact_matrix(self, chain_a="A", chain_b=None, cutoff=CONTACTS_CUTOFF_NM, moieties=DEFAULT_CONTACT_MOIETIES,
                       scheme="ca"):
        """Return the residue-pair minimum-distance matrix in nm, NaN where residues are not in contact.

        Args:
            chain_a: Chain on the rows, standard or author ID.
            chain_b: Chain on the columns, or None for an intra-chain map.
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.

        Returns:
            Dict with 'matrix', 'labels' (standard IDs), 'keys', and 'contacts'.
        """
        chain_a, chain_b, label_a, label_b = self._pair_keys(chain_a, chain_b)
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all", chain_ids=sorted({chain_a, chain_b}),
                                         scheme=scheme, moieties=moieties)
        return {"matrix": self._matrix_from_contacts(contacts, chain_a, chain_b), "labels": (label_a, label_b),
                "keys": (chain_a, chain_b), "contacts": contacts}

    def complete_contacts(self, chain_a="A", chain_b="B", cutoff=CONTACTS_CUTOFF_NM,
                          moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca", min_separation=COMPLETE_MIN_SEPARATION):
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
        """
        k1, k2, _, _ = self._pair_keys(chain_a, chain_b)
        if k1 == k2:
            raise ValueError("chain_a and chain_b must be different chains for a block map.")
        block = self._block_matrix([k1, k2], cutoff, scheme, moieties)
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
            "keys": (k1, k2),
            "contacts": block["contacts"],
        }
