"""Pairwise aligners and local alignment of observed chain sequences against their reference sequences with both chain IDs."""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd
from Bio.Align import PairwiseAligner

from .constants import (
    ALN_SCORE_GAP_WEIGHT, ALN_SCORE_MISMATCH_WEIGHT, GLOBAL_GAP_EXTEND_SCORE, GLOBAL_GAP_OPEN_SCORE,
    GLOBAL_MATCH_SCORE, GLOBAL_MISMATCH_SCORE,
)
from ._types import Record
from .base import StructureState

AlignmentResult = tuple[Record | None, list[Record]]

###################################### ALIGNERS


def make_global_aligner() -> PairwiseAligner:
    """Return a global aligner with the fixed scoring scheme used for reference mapping.

    Returns:
        A global PairwiseAligner.
    """
    return PairwiseAligner(
        mode="global",
        match_score=GLOBAL_MATCH_SCORE,
        mismatch_score=GLOBAL_MISMATCH_SCORE,
        open_gap_score=GLOBAL_GAP_OPEN_SCORE,
        extend_gap_score=GLOBAL_GAP_EXTEND_SCORE,
    )


def make_local_aligner(config: Mapping[str, Any]) -> PairwiseAligner:
    """Return a local aligner scored from the config.

    Args:
        config: Config holding the aln_* scores.

    Returns:
        A local PairwiseAligner.
    """
    return PairwiseAligner(
        mode="local",
        match_score=config["aln_match_score"],
        mismatch_score=config["aln_mismatch_score"],
        open_gap_score=config["aln_open_gap_score"],
        extend_gap_score=config["aln_extend_gap_score"],
    )


def query_coverage(aligned: int, query_length: int) -> float:
    """Return the percent of query residues aligned to a reference residue.

    Args:
        aligned: Number of alignment columns without a gap.
        query_length: Length of the query sequence.

    Returns:
        Coverage in percent, 0.0 for an empty query.
    """
    return round(aligned / query_length * 100, 2) if query_length else 0.0

###################################### LOCAL ALIGNMENT


class AlignmentMixin(StructureState):
    """Local alignments of observed chain sequences against their reference sequences."""

    def _local_alignment(self, chain_key: str) -> AlignmentResult:
        """Return the cached local alignment of a chain, recomputing it if its sequences changed.

        Args:
            chain_key: Chain key as held by the structure.

        Returns:
            (stats, position_rows); stats is None when there is no alignment.
        """
        reference = self.get_reference_sequence(chain_key)
        observed = self.get_sequence(chain_key)
        key = (reference, observed)
        cached = self._alignment_cache.get(chain_key)
        if cached and cached[0] == key:
            return cached[1]
        result = self._compute_local_alignment(chain_key, reference, observed)
        self._alignment_cache[chain_key] = (key, result)
        return result

    def _compute_local_alignment(self, chain_key: str, reference: str, observed: str) -> AlignmentResult:
        """Compute the local alignment of a chain's observed sequence against its reference.

        Args:
            chain_key: Chain key.
            reference: Reference sequence.
            observed: Observed sequence.

        Returns:
            (stats, position_rows) with 1-based positions; (None, []) when there is no alignment.
        """
        if not reference or not observed:
            return None, []
        alignments = self._local_aligner.align(reference, observed)
        if not alignments:
            return None, []
        record = self.chains[chain_key]
        label, auth = record.chain_id, record.auth_chain_id
        ref_aligned, obs_aligned = alignments[0]
        coordinates = alignments[0].coordinates
        ref_start, obs_start = int(coordinates[0, 0]), int(coordinates[1, 0])
        ref_end, obs_end = int(coordinates[0, -1]), int(coordinates[1, -1])
        rows, counts, gap_positions = self._alignment_columns(label, auth, ref_aligned, obs_aligned, ref_start, obs_start)
        stats = self._alignment_stats(label, auth, reference, observed, counts, gap_positions,
                                      (ref_start, ref_end), obs_end, obs_start)
        return stats, rows

    def _alignment_columns(self, label: str, auth: str, ref_aligned: str, obs_aligned: str, ref_start: int,
                           obs_start: int) -> tuple[list[Record], dict[str, int], list[str]]:
        """Tabulate the aligned columns of two sequences.

        Args:
            label: Standard chain ID.
            auth: Author chain ID.
            ref_aligned: Gapped reference string of the alignment.
            obs_aligned: Gapped observed string of the alignment.
            ref_start: 0-based start of the aligned reference region.
            obs_start: 0-based start of the aligned observed region.

        Returns:
            (position_rows, counts, gap_positions).
        """
        counts = {"matches": 0, "mismatches": 0, "insertions": 0, "deletions": 0}
        gap_positions: list[str] = []
        rows: list[Record] = []
        ref_pos, obs_pos = ref_start, obs_start
        for column, (ref_char, obs_char) in enumerate(zip(ref_aligned, obs_aligned), start=1):
            if ref_char != "-":
                ref_pos += 1
            if obs_char != "-":
                obs_pos += 1
            if ref_char == "-":
                counts["insertions"] += 1
            elif obs_char == "-":
                counts["deletions"] += 1
                gap_positions.append(str(obs_pos))
            elif ref_char == obs_char:
                counts["matches"] += 1
            else:
                counts["mismatches"] += 1
            rows.append({
                "structure_id": self.structure_id,
                "chain_id": label,
                "auth_chain_id": auth,
                "alignment_position": column,
                "ref": ref_char,
                "alt": obs_char,
                "posref": ref_pos if ref_char != "-" else None,
                "posalt": obs_pos if obs_char != "-" else None,
                "match": int(ref_char == obs_char and ref_char != "-"),
            })
        return rows, counts, gap_positions

    def _alignment_stats(self, label: str, auth: str, reference: str, observed: str, counts: dict[str, int],
                         gap_positions: list[str], ref_span: tuple[int, int], obs_end: int, obs_start: int) -> Record:
        """Return the one-row summary of a local alignment.

        identity_percent divides matches by all alignment columns, gaps included. query_coverage_percent
        divides matched and mismatched residues by the observed length.

        Args:
            label: Standard chain ID.
            auth: Author chain ID.
            reference: Reference sequence.
            observed: Observed sequence.
            counts: Matches, mismatches, insertions, and deletions.
            gap_positions: 1-based observed positions of deletions in the query.
            ref_span: (start, end) 0-based reference positions of the aligned region.
            obs_end: 0-based end of the aligned observed region.
            obs_start: 0-based start of the aligned observed region.

        Returns:
            The summary row.
        """
        ref_start, ref_end = ref_span
        matches, mismatches = counts["matches"], counts["mismatches"]
        insertions, deletions = counts["insertions"], counts["deletions"]
        aligned_length = sum(counts.values())
        return {
            "structure_id": self.structure_id,
            "chain_id": label,
            "auth_chain_id": auth,
            "query_id": f"{self.structure_id}_{label}",
            "ref_id": f"{self.structure_id}_{label}_reference",
            "query_length": len(observed),
            "ref_length": len(reference),
            "aligned_length": aligned_length,
            "matches": matches,
            "mismatches": mismatches,
            "insertions_in_query": insertions,
            "deletions_in_query": deletions,
            "identity_percent": round(matches / aligned_length * 100, 2) if aligned_length else 0.0,
            "query_coverage_percent": query_coverage(matches + mismatches, len(observed)),
            "alignment_score": round(matches - ALN_SCORE_MISMATCH_WEIGHT * mismatches
                                     - ALN_SCORE_GAP_WEIGHT * (insertions + deletions), 2),
            "range": f"{obs_start + 1}-{obs_end}",
            "ref_start": ref_start,
            "ref_end": ref_end,
            "gappos": ",".join(gap_positions),
            "sequence": reference[ref_start:ref_end],
        }

###################################### EXPOSED ALIGNMENT FRAMES

    @property
    def reference_alignment_stats(self) -> pd.DataFrame:
        """Return local alignment statistics of each chain against its reference.

        Returns:
            One row per aligned chain.
        """
        rows = [self._local_alignment(chain_id)[0] for chain_id in self.chains]
        return pd.DataFrame([row for row in rows if row])

    @property
    def reference_alignment_positions(self) -> pd.DataFrame:
        """Return the column-by-column local alignment of every chain against its reference.

        Returns:
            One row per alignment column.
        """
        rows: list[Record] = []
        for chain_id in self.chains:
            rows.extend(self._local_alignment(chain_id)[1])
        return pd.DataFrame(rows)
