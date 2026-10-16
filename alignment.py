#!/usr/bin/env python
"""Local alignment of observed chain sequences against their reference sequences with both chain IDs."""

###################################### IMPORTS AND CONSTANTS

import pandas as pd

###################################### LOCAL ALIGNMENT


class AlignmentMixin:
    def _local_alignment(self, chain_id):
        """Return cached (stats, position_rows, sequence_to_reference) for a chain, recomputing if sequences changed."""
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
        """Compute a local alignment and return (stats, position_rows, sequence_to_reference).

        Args:
            chain_id: Chain key.
            reference: Reference sequence.
            observed: Observed sequence.
        """
        empty = (None, [], {})
        if not reference or not observed:
            return empty
        alignments = self._local_aligner.align(reference, observed)
        if not alignments:
            return empty
        data = self.chains[chain_id]
        label, auth = data["chain_id"], data["auth_chain_id"]
        best = alignments[0]
        ref_aligned, obs_aligned = best
        ref_start, obs_start = int(best.coordinates[0, 0]), int(best.coordinates[1, 0])
        ref_end, obs_end = int(best.coordinates[0, 1]), int(best.coordinates[1, 1])
        counts = {"matches": 0, "mismatches": 0, "insertions": 0, "deletions": 0}
        gap_positions, rows, seq_to_ref = [], [], {}
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
            if r != "-" and o != "-":
                seq_to_ref[obs_pos] = ref_pos
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
        matches = counts["matches"]
        mismatches = counts["mismatches"]
        insertions = counts["insertions"]
        deletions = counts["deletions"]
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
            "query_coverage_percent": round((len(observed) - deletions) / len(observed) * 100, 2),
            "alignment_score": round(matches - 0.5 * mismatches - 1.5 * (insertions + deletions), 2),
            "range": f"{obs_start + 1}-{obs_end}",
            "ref_start": ref_start,
            "ref_end": ref_end,
            "gappos": ",".join(gap_positions),
            "sequence": reference[ref_start:ref_end],
        }
        return stats, rows, seq_to_ref

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