"""Contact statistics: per chain pair coverage and density, per-residue profiles, sequence separation, and moiety pairs."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

from collections import Counter
from collections.abc import Iterator
from typing import Any

import numpy as np
import pandas as pd

from ._types import ChainRefs, Moieties
from .base import StructureState
from .constants import (
    CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES, MOIETY_PAIR_COLUMNS, PAIR_STAT_COLUMNS, PROFILE_COLUMNS,
    SEPARATION_COLUMNS,
)

ResidueToken = int | tuple[str, int, str]


def _residue_token(index: Any, resnum: int, resname: str) -> ResidueToken:
    """Identify a residue: by sequence index if polymer, by (resnum, resname) if hetero.

    Args:
        index: 1-based sequence index, or NA for hetero residues.
        resnum: Residue number.
        resname: Residue name.

    Returns:
        The sequence index, or ('hetero', resnum, resname).
    """
    return ("hetero", int(resnum), resname) if pd.isna(index) else int(index)


def _contact_sides(contacts: pd.DataFrame) -> Iterator[tuple[str, ResidueToken, float]]:
    """Walk both residues of every contact row.

    Args:
        contacts: Contacts table.

    Yields:
        (chain label, residue token, distance in nm).
    """
    for row in contacts.to_dict("records"):
        distance = float(row["min_distance_nm"])
        for side in ("1", "2"):
            token = _residue_token(row[f"sequence_index_{side}"], int(row[f"resnum_{side}"]), str(row[f"resname_{side}"]))
            yield str(row[f"chain_{side}"]), token, distance

###################################### STATISTICS


class ContactStatsMixin(StructureState):
    """Contact statistics per chain pair, per contact-point residue, per sequence separation, and per moiety pair."""

    def _pair_table(self, contacts: pd.DataFrame, moieties: Moieties) -> pd.DataFrame:
        """Summarize contacts per chain pair (or per chain, for intra-chain contacts).

        Coverage uses the number of contact-point residues in each chain, so moieties left out do not dilute it.

        Args:
            contacts: Contacts table over all chains.
            moieties: Biomolecule types used for the contacts.

        Returns:
            One row per pair with counts, coverage, density, and distances.
        """
        if contacts.empty:
            return pd.DataFrame(columns=PAIR_STAT_COLUMNS)
        n_points = Counter(record["chain_id"] for record in self._contact_records(moieties=moieties))
        pairs = sorted({tuple(sorted(p)) for p in zip(contacts["chain_1"], contacts["chain_2"])})
        rows = []
        for label_a, label_b in pairs:
            same = label_a == label_b
            pair_rows = self._pair_contacts(contacts, label_a, label_b)
            contacted_by_chain: dict[str, set[ResidueToken]] = {}
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

    def _residue_profile_table(self, contacts: pd.DataFrame, moieties: Moieties,
                         chain_keys: ChainRefs | None = None) -> pd.DataFrame:
        """Count contacts and the nearest partner distance of every contact-point residue.

        Args:
            contacts: Contacts table.
            moieties: Biomolecule types to include.
            chain_keys: Chain keys to profile, or None for all.

        Returns:
            One row per contact-point residue.
        """
        counts: Counter[tuple[str, ResidueToken]] = Counter()
        nearest: dict[tuple[str, ResidueToken], float] = {}
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
    def _separation_table(contacts: pd.DataFrame) -> pd.DataFrame:
        """Tabulate intra-chain polymer contacts with their sequence separation.

        Args:
            contacts: Contacts table.

        Returns:
            One row per intra-chain contact with |i - j| and its distance.
        """
        intra = contacts[contacts["chain_1"] == contacts["chain_2"]].dropna(subset=["sequence_index_1", "sequence_index_2"])
        return pd.DataFrame({
            "chain_id": intra["chain_1"].to_numpy(),
            "separation": (intra["sequence_index_2"] - intra["sequence_index_1"]).abs().to_numpy(dtype=int),
            "min_distance_nm": intra["min_distance_nm"].to_numpy(),
        }, columns=SEPARATION_COLUMNS)

    @staticmethod
    def _moiety_pair_table(contacts: pd.DataFrame) -> pd.DataFrame:
        """Group contacts by moiety pair, such as protein/nucleic.

        Args:
            contacts: Contacts table.

        Returns:
            Contact counts and distance summary per moiety pair.
        """
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

    def _stat_tables_from_contacts(self, contacts: pd.DataFrame, moieties: Moieties,
                                  chain_keys: ChainRefs | None = None) -> dict[str, pd.DataFrame]:
        """Compute the contact_stat_tables tables from precomputed contacts.

        Args:
            contacts: Contacts table.
            moieties: Biomolecule types used for the contacts.
            chain_keys: Chain keys to profile, or None for all.

        Returns:
            The same dict as contact_stat_tables.
        """
        return {
            "pairs": self._pair_table(contacts, moieties),
            "profile": self._residue_profile_table(contacts, moieties, chain_keys),
            "separation": self._separation_table(contacts),
            "moiety_pairs": self._moiety_pair_table(contacts),
        }

    def contact_stat_tables(self, cutoff: float = CONTACTS_CUTOFF_NM, moieties: Moieties = DEFAULT_CONTACT_MOIETIES,
                           scheme: str = "ca") -> dict[str, pd.DataFrame]:
        """Compute contact statistics over all chains.

        Args:
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.

        Returns:
            Dict of tables: 'pairs' (per chain pair counts, coverage, density, distances), 'profile' (per
            contact-point residue count and nearest distance), 'separation' (intra-chain contacts with sequence
            separation), and 'moiety_pairs' (counts and distances per moiety pair).
        """
        contacts = self.contact_table(cutoff=cutoff, contact_type="all", scheme=scheme, moieties=moieties)
        return self._stat_tables_from_contacts(contacts, moieties)
