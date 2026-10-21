"""Contact statistics: per chain pair coverage and density, per-residue profiles, sequence separation, and moiety pairs."""

###################################### IMPORTS AND CONSTANTS

from collections import Counter

import numpy as np
import pandas as pd

from .constants import (
    CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES, MOIETY_PAIR_COLUMNS, PAIR_STAT_COLUMNS, PROFILE_COLUMNS,
    SEPARATION_COLUMNS,
)


def _residue_token(index, resnum, resname):
    """Return a residue identifier: the sequence index for polymer residues, (resnum, resname) for hetero ones."""
    return ("hetero", int(resnum), resname) if pd.isna(index) else int(index)


def _contact_sides(contacts):
    """Yield (label, residue_id, distance) for both residues of every contact row."""
    for row in contacts.itertuples(index=False):
        yield row.chain_1, _residue_token(row.sequence_index_1, row.resnum_1, row.resname_1), row.min_distance_nm
        yield row.chain_2, _residue_token(row.sequence_index_2, row.resnum_2, row.resname_2), row.min_distance_nm

###################################### STATISTICS


class ContactStatsMixin:
    """Contact statistics per chain pair, per contact-point residue, per sequence separation, and per moiety pair."""

    def _pair_statistics(self, contacts, moieties):
        """Return one row per chain pair (or chain, for intra-chain contacts) with counts, coverage, density, and distances.

        Coverage uses the number of contact-point residues in each chain, so moieties left out do not dilute it.

        Args:
            contacts: Contacts DataFrame over all chains.
            moieties: Biomolecule types used for the contacts.
        """
        if contacts.empty:
            return pd.DataFrame(columns=PAIR_STAT_COLUMNS)
        n_points = Counter(record["chain_id"] for record in self._contact_records(moieties=moieties))
        pairs = sorted({tuple(sorted(p)) for p in zip(contacts["chain_1"], contacts["chain_2"])})
        rows = []
        for label_a, label_b in pairs:
            same = label_a == label_b
            pair_rows = self._pair_rows(contacts, label_a, label_b)
            contacted_by_chain = {}
            for label, residue, _ in _contact_sides(pair_rows):
                contacted_by_chain.setdefault(label, set()).add(residue)
            seen_a, seen_b = contacted_by_chain.get(label_a, set()), contacted_by_chain.get(label_b, set())
            n_a, n_b = n_points[label_a], n_points[label_b]
            possible = n_a * (n_a - 1) / 2 if same else n_a * n_b
            distances = pair_rows["min_distance_nm"]
            rows.append({
                "chain_a": label_a,
                "chain_b": label_b,
                "contact_type": "intra" if same else "inter",
                "contact_pairs": len(pair_rows),
                "contacted_a": len(seen_a),
                "contacted_b": len(seen_b),
                "fraction_a": len(seen_a) / n_a if n_a else np.nan,
                "fraction_b": len(seen_b) / n_b if n_b else np.nan,
                "contact_density": len(pair_rows) / possible if possible else np.nan,
                "min_nm": round(float(distances.min()), 4),
                "median_nm": round(float(distances.median()), 4),
                "mean_nm": round(float(distances.mean()), 4),
            })
        return pd.DataFrame(rows, columns=PAIR_STAT_COLUMNS)

    def _residue_profile(self, contacts, moieties, chain_keys=None):
        """Return one row per contact-point residue with its number of contacts and nearest partner distance.

        Args:
            contacts: Contacts DataFrame.
            moieties: Biomolecule types to include.
            chain_keys: Chain keys to profile, or None for all.
        """
        counts, nearest = Counter(), {}
        for label, residue, distance in _contact_sides(contacts):
            counts[(label, residue)] += 1
            nearest[(label, residue)] = min(nearest.get((label, residue), np.inf), float(distance))
        rows = []
        for record in self._contact_records(chain_keys, moieties):
            index = record["sequence_index"]
            number = None if index is None else index + 1
            point = (record["chain_id"], _residue_token(np.nan if index is None else number, record["resnum"], record["resname"]))
            rows.append({
                "structure_id": self.structure_id,
                "chain_id": record["chain_id"],
                "auth_chain_id": record["auth_chain_id"],
                "sequence_index": number,
                "resnum": record["resnum"],
                "resname": record["resname"],
                "moiety": record["moiety"],
                "n_contacts": counts.get(point, 0),
                "nearest_nm": round(nearest[point], 4) if point in nearest else np.nan,
            })
        return pd.DataFrame(rows, columns=PROFILE_COLUMNS).astype({"sequence_index": "Int64"})

    @staticmethod
    def _separation_frame(contacts):
        """Return intra-chain polymer contacts with their sequence separation |i - j|."""
        intra = contacts[contacts["chain_1"] == contacts["chain_2"]].dropna(subset=["sequence_index_1", "sequence_index_2"])
        return pd.DataFrame({
            "chain_id": intra["chain_1"].to_numpy(),
            "separation": (intra["sequence_index_2"] - intra["sequence_index_1"]).abs().to_numpy(dtype=int),
            "min_distance_nm": intra["min_distance_nm"].to_numpy(),
        }, columns=SEPARATION_COLUMNS)

    @staticmethod
    def _moiety_pair_frame(contacts):
        """Return contact counts and distances grouped by moiety pair, such as protein/nucleic."""
        if contacts.empty:
            return pd.DataFrame(columns=MOIETY_PAIR_COLUMNS)
        frame = pd.DataFrame({
            "moiety_pair": ["/".join(sorted(p)) for p in zip(contacts["moiety_1"], contacts["moiety_2"])],
            "min_distance_nm": contacts["min_distance_nm"].to_numpy(),
        })
        out = frame.groupby("moiety_pair").agg(
            contact_pairs=("min_distance_nm", "size"),
            median_nm=("min_distance_nm", "median"),
            mean_nm=("min_distance_nm", "mean"),
        ).reset_index()
        out[["median_nm", "mean_nm"]] = out[["median_nm", "mean_nm"]].round(4)
        return out[MOIETY_PAIR_COLUMNS]

    def _statistics_from_contacts(self, contacts, moieties, chain_keys=None):
        """Return the contact_statistics dict for precomputed contacts, profiling only chain_keys if given."""
        return {
            "pairs": self._pair_statistics(contacts, moieties),
            "profile": self._residue_profile(contacts, moieties, chain_keys),
            "separation": self._separation_frame(contacts),
            "moiety_pairs": self._moiety_pair_frame(contacts),
        }

    def contact_statistics(self, cutoff=CONTACTS_CUTOFF_NM, moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca"):
        """Return contact statistics as a dict of DataFrames.

        Keys:
            'pairs': per chain pair, contact counts, residue coverage, density, and distance summary.
            'profile': per contact-point residue, contact count and nearest partner distance.
            'separation': intra-chain contacts with sequence separation.
            'moiety_pairs': contact counts and distances per moiety pair.

        Args:
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
        """
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all", scheme=scheme, moieties=moieties)
        return self._statistics_from_contacts(contacts, moieties)
