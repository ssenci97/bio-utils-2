"""Pairwise aligners and local alignment of observed chain sequences against their reference sequences with both chain IDs."""

###################################### IMPORTS AND CONSTANTS

import pandas as pd
from Bio.Align import PairwiseAligner

from .constants import (
    ALN_SCORE_GAP_WEIGHT, ALN_SCORE_MISMATCH_WEIGHT, GLOBAL_GAP_EXTEND_SCORE, GLOBAL_GAP_OPEN_SCORE,
    GLOBAL_MATCH_SCORE, GLOBAL_MISMATCH_SCORE,
)

###################################### ALIGNERS


def make_global_aligner():
    """Return a global PairwiseAligner with the fixed scoring scheme used for reference mapping."""
    return PairwiseAligner(
        mode="global",
        match_score=GLOBAL_MATCH_SCORE,
        mismatch_score=GLOBAL_MISMATCH_SCORE,
        open_gap_score=GLOBAL_GAP_OPEN_SCORE,
        extend_gap_score=GLOBAL_GAP_EXTEND_SCORE,
    )


def make_local_aligner(config):
    """Return a local PairwiseAligner scored from the config dict."""
    return PairwiseAligner(
        mode="local",
        match_score=config["aln_match_score"],
        mismatch_score=config["aln_mismatch_score"],
        open_gap_score=config["aln_open_gap_score"],
        extend_gap_score=config["aln_extend_gap_score"],
    )


def query_coverage(aligned, query_length):
    """Return the percent of query residues aligned to a reference residue.

    Args:
        aligned: Number of alignment columns without a gap.
        query_length: Length of the query sequence.
    """
    return round(aligned / query_length * 100, 2) if query_length else 0.0

###################################### LOCAL ALIGNMENT


class AlignmentMixin:
    def _local_alignment(self, chain_id):
        """Return cached (stats, position_rows) for a chain, recomputing if sequences changed."""
        reference = self.get_reference_sequence(chain_id)
        observed = self.get_sequence(chain_id)
        key = (reference, observed)
        cached = self._alignment_cache.get(chain_id)
        if cached and cached[0] == key:
            return cached[1]
        result = self._compute_local_alignment(chain_id, reference, observed)
        self._alignment_cache[chain_id] = (key, result)
        return result

    def _compute_local_alignment(self, chain_id, reference, observed):
        """Compute a local alignment and return (stats, position_rows) with 1-based positions.

        Args:
            chain_id: Chain key.
            reference: Reference sequence.
            observed: Observed sequence.
        """
        if not reference or not observed:
            return None, []
        alignments = self._local_aligner.align(reference, observed)
        if not alignments:
            return None, []
        data = self.chains[chain_id]
        label, auth = data["chain_id"], data["auth_chain_id"]
        best = alignments[0]
        ref_aligned, obs_aligned = best
        ref_start, obs_start = int(best.coordinates[0, 0]), int(best.coordinates[1, 0])
        ref_end, obs_end = int(best.coordinates[0, -1]), int(best.coordinates[1, -1])
        counts = {"matches": 0, "mismatches": 0, "insertions": 0, "deletions": 0}
        gap_positions, rows = [], []
        ref_pos, obs_pos = ref_start, obs_start
        for column, (r, o) in enumerate(zip(ref_aligned, obs_aligned), start=1):
            if r != "-":
                ref_pos += 1
            if o != "-":
                obs_pos += 1
            if r == "-":
                counts["insertions"] += 1
            elif o == "-":
                counts["deletions"] += 1
                gap_positions.append(str(obs_pos))
            elif r == o:
                counts["matches"] += 1
            else:
                counts["mismatches"] += 1
            rows.append({
                "structure_id": self.structure_id,
                "chain_id": label,
                "auth_chain_id": auth,
                "alignment_position": column,
                "ref": r,
                "alt": o,
                "posref": ref_pos if r != "-" else None,
                "posalt": obs_pos if o != "-" else None,
                "match": int(r == o and r != "-"),
            })
        matches, mismatches = counts["matches"], counts["mismatches"]
        insertions, deletions = counts["insertions"], counts["deletions"]
        aligned_length = sum(counts.values())
        stats = {
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
        return stats, rows

###################################### EXPOSED ALIGNMENT FRAMES

    @property
    def alignment_stats(self):
        """Return a DataFrame with one row of local alignment statistics per chain."""
        rows = [self._local_alignment(chain_id)[0] for chain_id in self.chains]
        return pd.DataFrame([row for row in rows if row])

    @property
    def alignment_positions(self):
        """Return a DataFrame with the column-by-column local alignment of every chain."""
        rows = []
        for chain_id in self.chains:
            rows.extend(self._local_alignment(chain_id)[1])
        return pd.DataFrame(rows)
