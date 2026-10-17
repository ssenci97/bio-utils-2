#!/usr/bin/env python
"""Contact maps and statistics for chains: pair maps with grayscale background contacts, two-chain, three-chain, and four-chain block maps with chain outlines and optional pLDDT tracks, tick and label controls, contact statistics, and a summary panel."""

###################################### IMPORTS AND CONSTANTS

from collections import Counter
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch, Rectangle

from constants import DEFAULT_CONTACT_MOIETIES, DEFAULT_CUTOFF
from helpers import residue_moiety

CONTACT_MAP_BACKGROUND = "white"
CONTACT_MAP_COLOR = "#c0392b"
CONTACT_MAP_CMAP = "viridis_r"
CONTACT_MAP_SIZE_IN = 3.5
CONTACT_FOUR_SIZE_IN = 4.0
CONTACT_MAP_DPI = 150
CONTACT_FONT_PT_PER_IN = 4.0
CONTACT_FONT_MIN_PT = 10.0
CONTACT_FONT_MAX_PT = 24.0
CONTACT_TICK_RATIO = 0.85
CONTACT_TITLE_RATIO = 1.15
CONTACT_SUPTITLE_RATIO = 1.3
CONTACT_LEGEND_ROW_RATIO = 1.4
CONTACT_TICK_BLOCK_RATIO = 1.4
CONTACT_AXIS_TITLE_RATIO = 1.6
CONTACT_MAP_TICK_STEP = 50
CONTACT_MAP_N_TICKS = 4
CONTACT_MAP_MIN_MARKER_PT = 2.5
CONTACT_MAP_STYLES = ("scatter", "heatmap")
CONTACT_MAP_VALUES = ("contact", "distance")
CONTACT_GRAY_FACTOR = 1.5
CONTACT_GRAY_NEAR = "#8c8c8c"
CONTACT_GRAY_FAR = "#f2f2f2"
CONTACT_GRAY_CMAP = LinearSegmentedColormap.from_list("contact_gray", [CONTACT_GRAY_NEAR, CONTACT_GRAY_FAR])
TWO_INTRA_COLOR = "#7f7f7f"
TWO_BOX_COLOR = "black"
TWO_BOX_LW = 0.5
TWO_TRACK_OFFSET = 1.5
THREE_PAIRS = ((0, 1), (0, 2), (1, 2))
THREE_PAIR_COLORS = ("#c0392b", "#2980b9", "#27ae60")
FOUR_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
FOUR_PAIR_COLORS = ("#c0392b", "#2980b9", "#27ae60", "#8e44ad", "#d35400", "#16a085")
PLDDT_BANDS = (("≥90", "#0053D6"), ("70–90", "#65CBF3"), ("50–70", "#FFDB13"), ("<50", "#FF7D45"))
PLDDT_MISSING_COLOR = "#d0d0d0"
COMPLETE_MIN_SEPARATION = 2
CONTACT_PANEL_SIZE_IN = 10.0
CONTACT_PANEL_LINE_COLORS = ("#1f77b4", "#ff7f0e")
CONTACT_HIST_BINS = 30
PAIR_STAT_COLUMNS = [
    "chain_a", "chain_b", "contact_type", "contact_pairs", "contacted_a", "contacted_b",
    "fraction_a", "fraction_b", "contact_density", "min_nm", "median_nm", "mean_nm",
]
PROFILE_COLUMNS = [
    "structure_id", "chain_id", "auth_chain_id", "sequence_index", "resnum", "resname",
    "moiety", "n_contacts", "nearest_nm",
]
SEPARATION_COLUMNS = ["chain_id", "separation", "min_distance_nm"]
MOIETY_PAIR_COLUMNS = ["moiety_pair", "contact_pairs", "median_nm", "mean_nm"]


def _font_size(size_in, font_size=None):
    """Return the base font size in points, proportional to the plot size unless given explicitly.

    Args:
        size_in: Plot width and height in inches.
        font_size: Explicit base font size in points, or None to scale with size_in.
    """
    if font_size is not None:
        return float(font_size)
    return float(np.clip(CONTACT_FONT_PT_PER_IN * size_in, CONTACT_FONT_MIN_PT, CONTACT_FONT_MAX_PT))


def _legend_offset(size_in, font_size, n_entries, ncol):
    """Return the legend anchor in axes coordinates, placed below tick labels and the axis title.

    Args:
        size_in: Map width and height in inches.
        font_size: Base font size in points.
        n_entries: Number of legend entries.
        ncol: Number of legend columns.
    """
    n_rows = math.ceil(n_entries / ncol)
    below_pt = font_size * (CONTACT_TICK_BLOCK_RATIO + CONTACT_AXIS_TITLE_RATIO + CONTACT_LEGEND_ROW_RATIO * n_rows)
    return -below_pt / (size_in * 72 * 0.8)


def _marker_points(map_size_in, n_cols, n_rows, min_marker_pt):
    """Return the marker side in points so that one residue is at least min_marker_pt wide.

    Args:
        map_size_in: Width and height of the map axes in inches.
        n_cols: Number of columns in the matrix.
        n_rows: Number of rows in the matrix.
        min_marker_pt: Smallest allowed marker side in points.
    """
    usable_pt = map_size_in * 72 * 0.8
    cell_pt = usable_pt / max(n_cols, n_rows, 1)
    return max(min_marker_pt, cell_pt)


def _gray_limit(cutoff, gray_cutoff):
    """Return the upper distance in nm for the grayscale background.

    Args:
        cutoff: Colored-contact cutoff in nm.
        gray_cutoff: Explicit background limit in nm, or None to use CONTACT_GRAY_FACTOR times cutoff.
    """
    if gray_cutoff is None:
        return cutoff * CONTACT_GRAY_FACTOR
    return max(gray_cutoff, cutoff)


def _tick_positions(n, tick_step, n_ticks):
    """Return tick positions along an axis of n residues.

    With n_ticks, positions are spread evenly across the axis, including both ends. Otherwise they
    are spaced tick_step residues apart from the first residue.

    Args:
        n: Number of positions along the axis.
        tick_step: Spacing in residues, used when n_ticks is None.
        n_ticks: Number of ticks per axis, or None.
    """
    if n <= 0:
        return np.array([], dtype=int)
    if n_ticks is not None:
        if n_ticks < 1:
            raise ValueError("n_ticks must be at least 1.")
        if n_ticks == 1:
            return np.array([0], dtype=int)
        return np.unique(np.linspace(0, n - 1, int(n_ticks)).round().astype(int))
    if tick_step < 1:
        raise ValueError("tick_step must be at least 1.")
    return np.arange(0, n, tick_step)


def plddt_colors(values):
    """Return AlphaFold pLDDT band colors for per-residue values, with gray where a value is missing.

    Args:
        values: Array of pLDDT values in 0 to 100, NaN where unknown.
    """
    values = np.asarray(values, dtype=float)
    colors = np.where(values >= 90, PLDDT_BANDS[0][1],
                      np.where(values >= 70, PLDDT_BANDS[1][1],
                               np.where(values >= 50, PLDDT_BANDS[2][1], PLDDT_BANDS[3][1])))
    return [PLDDT_MISSING_COLOR if np.isnan(v) else c for v, c in zip(values, colors)]


###################################### CONTACT MATRIX


class ContactMapMixin:
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
        """Return the contact rows that join two chains, in either order.

        Args:
            contacts: Contacts DataFrame from compute_contacts.
            label_a: Standard ID of the first chain.
            label_b: Standard ID of the second chain.
        """
        forward = (contacts["chain_1"] == label_a) & (contacts["chain_2"] == label_b)
        backward = (contacts["chain_1"] == label_b) & (contacts["chain_2"] == label_a)
        return contacts[forward | backward]

    def _matrix_from_contacts(self, contacts, chain_a, chain_b):
        """Return the residue-pair minimum-distance matrix in nm, with NaN for non-contacts.

        Rows follow chain_a and columns follow chain_b. Intra-chain matrices are symmetric.

        Args:
            contacts: Contacts DataFrame from compute_contacts, covering both chains.
            chain_a: First chain key.
            chain_b: Second chain key.
        """
        label_a, label_b = self.chains[chain_a]["chain_id"], self.chains[chain_b]["chain_id"]
        matrix = np.full((self.chains[chain_a]["sequence_length"], self.chains[chain_b]["sequence_length"]), np.nan)
        for row in self._pair_rows(contacts, label_a, label_b).itertuples(index=False):
            if row.chain_1 == label_a:
                i, j = int(row.sequence_index_1), int(row.sequence_index_2)
            else:
                i, j = int(row.sequence_index_2), int(row.sequence_index_1)
            matrix[i, j] = row.min_distance_nm
            if chain_a == chain_b:
                matrix[j, i] = row.min_distance_nm
        return matrix

    def _gray_matrix(self, chain_a, chain_b, gray_cutoff, moieties):
        """Return the pair distance matrix for all contacts up to gray_cutoff, with NaN beyond it.

        Args:
            chain_a: First chain key.
            chain_b: Second chain key.
            gray_cutoff: Upper distance limit in nm for the background.
            moieties: Biomolecule types to include.
        """
        contacts = self.compute_contacts(cutoff=gray_cutoff, contact_type="all",
                                         chain_ids=sorted({chain_a, chain_b}), moieties=moieties)
        return self._matrix_from_contacts(contacts, chain_a, chain_b)

    def contact_matrix(self, chain_a="A", chain_b=None, cutoff=DEFAULT_CUTOFF, moieties=DEFAULT_CONTACT_MOIETIES):
        """Return the residue-pair minimum-distance matrix in nm, with NaN where residues are not in contact.

        Args:
            chain_a: First chain identifier, standard or author.
            chain_b: Second chain identifier, or None for an intra-chain map.
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.

        Returns:
            Tuple of (matrix, label_a, label_b) with standard chain IDs.
        """
        chain_a, chain_b, label_a, label_b = self._pair_keys(chain_a, chain_b)
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all",
                                         chain_ids=sorted({chain_a, chain_b}), moieties=moieties)
        return self._matrix_from_contacts(contacts, chain_a, chain_b), label_a, label_b

    def complete_contacts(self, chain_a="A", chain_b="B", cutoff=DEFAULT_CUTOFF,
                          moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca",
                          min_separation=COMPLETE_MIN_SEPARATION):
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
        k1, k2, label1, label2 = self._pair_keys(chain_a, chain_b)
        if k1 == k2:
            raise ValueError("chain_a and chain_b must be different chains for a block map.")
        n1, n2 = self.chains[k1]["sequence_length"], self.chains[k2]["sequence_length"]
        size = n1 + n2
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all", chain_ids=[k1, k2],
                                         scheme=scheme, moieties=moieties)
        matrix = np.full((size, size), np.nan)
        if not contacts.empty:
            gi, gj, distance = self._block_entries(contacts, label1, n1)
            matrix[gi, gj] = distance
            matrix[gj, gi] = distance
        index = np.arange(size)
        same_block = (index < n1)[:, None] == (index < n1)[None, :]
        separation = np.abs(index[:, None] - index[None, :])
        contact = ~np.isnan(matrix)
        return {
            "matrix": matrix,
            "intra": contact & same_block & (separation >= min_separation),
            "inter": contact & ~same_block,
            "sizes": (n1, n2),
            "labels": (label1, label2),
            "keys": (k1, k2),
            "contacts": contacts,
        }

    @staticmethod
    def _block_entries(contacts, label1, n1):
        """Return block-matrix row indices, column indices, and distances for each contact.

        The first chain occupies indices 0 to n1-1, and the second chain follows it.

        Args:
            contacts: Contacts DataFrame covering both chains.
            label1: Standard ID of the first chain.
            n1: Number of contact-point positions in the first chain.
        """
        first_1 = (contacts["chain_1"] == label1).to_numpy()
        first_2 = (contacts["chain_2"] == label1).to_numpy()
        seq_1 = contacts["sequence_index_1"].to_numpy(dtype=int)
        seq_2 = contacts["sequence_index_2"].to_numpy(dtype=int)
        gi = np.where(first_1, seq_1, n1 + seq_1)
        gj = np.where(first_2, seq_2, n1 + seq_2)
        return gi, gj, contacts["min_distance_nm"].to_numpy(dtype=float)

    @staticmethod
    def _multi_block_entries(contacts, labels, offsets):
        """Return block-matrix row indices, column indices, and distances for contacts spanning any number of chains.

        Chain k occupies indices offsets[k] onward, in the order of labels.

        Args:
            contacts: Contacts DataFrame covering all chains in labels.
            labels: Standard IDs of the chains, in block order.
            offsets: Start index of each chain in the block matrix.
        """
        position = {label: k for k, label in enumerate(labels)}
        offsets = np.asarray(offsets, dtype=int)
        block_1 = contacts["chain_1"].map(position).to_numpy(dtype=int)
        block_2 = contacts["chain_2"].map(position).to_numpy(dtype=int)
        gi = offsets[block_1] + contacts["sequence_index_1"].to_numpy(dtype=int)
        gj = offsets[block_2] + contacts["sequence_index_2"].to_numpy(dtype=int)
        return gi, gj, contacts["min_distance_nm"].to_numpy(dtype=float)

    def _residue_plddt(self, chain_key):
        """Return per-residue pLDDT for a chain, read from the C-alpha B-factor, NaN where missing.

        AlphaFold models store pLDDT in the B-factor column. Experimental structures do not, so their
        values here are not confidence scores.

        Args:
            chain_key: Chain key.
        """
        residues = self.chains[chain_key]["residues"]
        n = self.chains[chain_key]["sequence_length"]
        values = np.full(n, np.nan)
        for position in range(min(n, len(residues))):
            residue = residues[position]
            if "CA" in residue:
                values[position] = residue["CA"].get_bfactor()
            else:
                atoms = list(residue.get_atoms())
                if atoms:
                    values[position] = np.mean([atom.get_bfactor() for atom in atoms])
        return values

###################################### STATISTICS

    def _pair_statistics(self, contacts, moieties):
        """Return one row per chain pair (or chain, for intra-chain contacts) with counts, coverage, density, and distances.

        Coverage uses the number of contact-point residues in each chain, so water and ions
        excluded by moieties do not dilute it.

        Args:
            contacts: Contacts DataFrame over all chains.
            moieties: Biomolecule types used for the contacts.
        """
        if contacts.empty:
            return pd.DataFrame(columns=PAIR_STAT_COLUMNS)
        key_of = {data["chain_id"]: key for key, data in self.chains.items()}
        n_points = {label: len(self._ca_records([key], moieties)) for label, key in key_of.items()}
        pairs = sorted({tuple(sorted(p)) for p in zip(contacts["chain_1"], contacts["chain_2"])})
        rows = []
        for label_a, label_b in pairs:
            same = label_a == label_b
            sub = self._pair_rows(contacts, label_a, label_b)
            seen_a = set(sub.loc[sub["chain_1"] == label_a, "sequence_index_1"]) | \
                set(sub.loc[sub["chain_2"] == label_a, "sequence_index_2"])
            if same:
                seen_b = seen_a
                possible = n_points[label_a] * (n_points[label_a] - 1) / 2
            else:
                seen_b = set(sub.loc[sub["chain_1"] == label_b, "sequence_index_1"]) | \
                    set(sub.loc[sub["chain_2"] == label_b, "sequence_index_2"])
                possible = n_points[label_a] * n_points[label_b]
            distances = sub["min_distance_nm"]
            rows.append({
                "chain_a": label_a,
                "chain_b": label_b,
                "contact_type": "intra" if same else "inter",
                "contact_pairs": len(sub),
                "contacted_a": len(seen_a),
                "contacted_b": len(seen_b),
                "fraction_a": len(seen_a) / n_points[label_a] if n_points[label_a] else np.nan,
                "fraction_b": len(seen_b) / n_points[label_b] if n_points[label_b] else np.nan,
                "contact_density": len(sub) / possible if possible else np.nan,
                "min_nm": round(float(distances.min()), 4),
                "median_nm": round(float(distances.median()), 4),
                "mean_nm": round(float(distances.mean()), 4),
            })
        return pd.DataFrame(rows, columns=PAIR_STAT_COLUMNS)

    def _residue_profile(self, contacts, moieties):
        """Return one row per contact-point residue with its number of contacts and nearest partner distance.

        Args:
            contacts: Contacts DataFrame over all chains.
            moieties: Biomolecule types to include.
        """
        counts, nearest = Counter(), {}
        for row in contacts.itertuples(index=False):
            for label, index in ((row.chain_1, row.sequence_index_1), (row.chain_2, row.sequence_index_2)):
                key = (label, int(index))
                counts[key] += 1
                nearest[key] = min(nearest.get(key, np.inf), float(row.min_distance_nm))
        rows = []
        for key, data in self.chains.items():
            label = data["chain_id"]
            for _, index, residue in self._ca_records([key], moieties):
                point = (label, int(index))
                rows.append({
                    "structure_id": self.structure_id,
                    "chain_id": label,
                    "auth_chain_id": data["auth_chain_id"],
                    "sequence_index": int(index),
                    "resnum": int(residue.id[1]),
                    "resname": residue.resname.strip().upper(),
                    "moiety": residue_moiety(residue),
                    "n_contacts": counts.get(point, 0),
                    "nearest_nm": round(nearest[point], 4) if point in nearest else np.nan,
                })
        return pd.DataFrame(rows, columns=PROFILE_COLUMNS)

    @staticmethod
    def _separation_frame(contacts):
        """Return intra-chain contacts with their sequence separation |i - j|.

        Args:
            contacts: Contacts DataFrame over all chains.
        """
        intra = contacts[contacts["chain_1"] == contacts["chain_2"]]
        return pd.DataFrame({
            "chain_id": intra["chain_1"].to_numpy(),
            "separation": (intra["sequence_index_2"] - intra["sequence_index_1"]).abs().to_numpy(),
            "min_distance_nm": intra["min_distance_nm"].to_numpy(),
        }, columns=SEPARATION_COLUMNS)

    @staticmethod
    def _moiety_pair_frame(contacts):
        """Return contact counts and distances grouped by the moiety pair, such as protein/nucleic.

        Args:
            contacts: Contacts DataFrame over all chains.
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

    def contact_statistics(self, cutoff=DEFAULT_CUTOFF, moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca"):
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
        return {
            "pairs": self._pair_statistics(contacts, moieties),
            "profile": self._residue_profile(contacts, moieties),
            "separation": self._separation_frame(contacts),
            "moiety_pairs": self._moiety_pair_frame(contacts),
        }

###################################### DRAWING HELPERS

    def _chain_name(self, chain_key):
        """Return the quoted author chain ID used in plot text, such as 'A'.

        Args:
            chain_key: Chain key.
        """
        return f"'{self.chains[chain_key]['auth_chain_id']}'"

    def _tick_labels(self, chain_key, n, tick_step, n_ticks=None):
        """Return tick positions and residue-number labels for one axis.

        Args:
            chain_key: Chain key whose residue numbers label the ticks.
            n: Number of positions along the axis.
            tick_step: Spacing between ticks in sequence positions, used when n_ticks is None.
            n_ticks: Number of ticks along the axis, or None.
        """
        residues = self.chains[chain_key]["residues"]
        positions = _tick_positions(n, tick_step, n_ticks)
        labels = [str(residues[p].id[1]) if p < len(residues) else str(p + 1) for p in positions]
        return positions, labels

    @staticmethod
    def _style_axes(ax, xpos, xlabels, ypos, ylabels, *, show_ticks, tick_size, tick_rotation,
                    xlabel, ylabel, font_size):
        """Apply tick positions, tick labels, and axis titles to an Axes.

        Args:
            ax: Matplotlib Axes.
            xpos, xlabels: X tick positions and labels.
            ypos, ylabels: Y tick positions and labels.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel, ylabel: Axis titles, already resolved.
            font_size: Base font size.
        """
        if show_ticks:
            ax.set_xticks(xpos, xlabels)
            ax.set_yticks(ypos, ylabels)
            ax.tick_params(labelsize=tick_size if tick_size is not None else font_size * CONTACT_TICK_RATIO)
            ax.tick_params(axis="x", labelrotation=tick_rotation)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        ax.set_xlabel(xlabel, fontsize=font_size)
        ax.set_ylabel(ylabel, fontsize=font_size)

    @staticmethod
    def _paint(ax, matrix, gray_matrix, value, style, background, color, cmap, cutoff, gray_cutoff, marker_pt):
        """Draw a pair matrix on ax, with an optional grayscale background, and return a mappable or None.

        Contacts within the cutoff are drawn on top of the background. Distance values use a fixed
        scale from 0 to cutoff, so plots with the same cutoff are comparable.

        Args:
            ax: Matplotlib Axes.
            matrix: Contact matrix with NaN for non-contacts.
            gray_matrix: Background distance matrix up to gray_cutoff, or None.
            value: 'contact' or 'distance'.
            style: 'scatter' (one marker per contact) or 'heatmap' (image).
            background: Background color.
            color: Contact color, used when value is 'contact'.
            cmap: Colormap name, used when value is 'distance'.
            cutoff: Distance cutoff in nm.
            gray_cutoff: Upper distance limit of the grayscale background in nm.
            marker_pt: Marker side in points, used by the scatter style.
        """
        if style == "heatmap":
            if gray_matrix is not None:
                gray = CONTACT_GRAY_CMAP.copy()
                gray.set_bad(background)
                ax.imshow(np.ma.masked_invalid(gray_matrix), cmap=gray, vmin=0, vmax=gray_cutoff,
                          interpolation="nearest")
            if value == "contact":
                return ax.imshow(np.ma.masked_invalid(np.where(np.isnan(matrix), np.nan, 1.0)),
                                 cmap=ListedColormap([color]), vmin=0, vmax=1, interpolation="nearest")
            colormap = mpl.colormaps[cmap].copy()
            colormap.set_bad(background)
            return ax.imshow(np.ma.masked_invalid(matrix), cmap=colormap, vmin=0, vmax=cutoff,
                             interpolation="nearest")
        size = marker_pt ** 2
        if gray_matrix is not None:
            rows, cols = np.nonzero(~np.isnan(gray_matrix))
            ax.scatter(cols, rows, c=gray_matrix[rows, cols], cmap=CONTACT_GRAY_CMAP, vmin=0, vmax=gray_cutoff,
                       s=size, marker="s", linewidths=0, zorder=1)
        rows, cols = np.nonzero(~np.isnan(matrix))
        if value == "contact":
            ax.scatter(cols, rows, s=size, marker="s", color=color, linewidths=0, zorder=2)
            return None
        return ax.scatter(cols, rows, c=matrix[rows, cols], cmap=cmap, vmin=0, vmax=cutoff,
                          s=size, marker="s", linewidths=0, zorder=2)

    def _draw_pair(self, ax, matrix, chain_a, chain_b, *, value, style, background, color, cmap,
                   cutoff, font_size, tick_step, min_marker_pt, map_size_in, title,
                   gray_matrix=None, gray_cutoff=None, n_ticks=None, show_ticks=True,
                   tick_size=None, tick_rotation=0, xlabel=None, ylabel=None):
        """Draw a pair contact map on ax with a fixed aspect ratio and quoted chain-ID axis titles.

        Args:
            ax: Matplotlib Axes.
            matrix: Contact matrix from _matrix_from_contacts.
            chain_a: Chain key on the y-axis.
            chain_b: Chain key on the x-axis.
            value, style, background, color, cmap, cutoff: See _paint.
            font_size: Resolved base font size in points.
            tick_step: Spacing of residue-number ticks.
            min_marker_pt: Smallest marker side in points.
            map_size_in: Size of the map axes in inches, used to size markers.
            title: Axes title.
            gray_matrix: Optional background distance matrix.
            gray_cutoff: Upper distance limit of the background in nm.
            n_ticks: Number of ticks per axis, or None to use tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for the quoted chain ID of chain_b.
            ylabel: Y-axis title, or None for the quoted chain ID of chain_a.
        """
        n_rows, n_cols = matrix.shape
        name_a, name_b = self._chain_name(chain_a), self._chain_name(chain_b)
        ax.set_facecolor(background)
        ax.set_xlim(-0.5, n_cols - 0.5)
        ax.set_ylim(n_rows - 0.5, -0.5)
        ax.set_aspect("equal", adjustable="box")
        marker_pt = _marker_points(map_size_in, n_cols, n_rows, min_marker_pt)
        mappable = self._paint(ax, matrix, gray_matrix, value, style, background, color, cmap,
                               cutoff, gray_cutoff, marker_pt)
        if mappable is not None and value == "distance":
            colorbar = ax.figure.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04)
            colorbar.set_label("min distance (nm)", fontsize=font_size)
            colorbar.ax.tick_params(labelsize=font_size * CONTACT_TICK_RATIO)
        xpos, xlabels = self._tick_labels(chain_b, n_cols, tick_step, n_ticks)
        ypos, ylabels = self._tick_labels(chain_a, n_rows, tick_step, n_ticks)
        self._style_axes(ax, xpos, xlabels, ypos, ylabels, show_ticks=show_ticks, tick_size=tick_size,
                         tick_rotation=tick_rotation,
                         xlabel=xlabel or f"{name_b}", ylabel=ylabel or f"{name_a}",
                         font_size=font_size)
        ax.set_title(title, fontsize=font_size * CONTACT_TITLE_RATIO)

    def _checked_matrix(self, matrix, chain_a, chain_b):
        """Return a precomputed matrix as float after checking it matches the two chains.

        Args:
            matrix: Matrix in the contact_matrix convention.
            chain_a: Chain key on the y-axis.
            chain_b: Chain key on the x-axis.
        """
        matrix = np.asarray(matrix, dtype=float)
        expected = (self.chains[chain_a]["sequence_length"], self.chains[chain_b]["sequence_length"])
        if matrix.shape != expected:
            raise ValueError(f"matrix shape {matrix.shape} does not match chains {expected}.")
        return matrix

###################################### PUBLIC PLOTS

    def plot_contact_map(self, chain_a="A", chain_b=None, cutoff=DEFAULT_CUTOFF,
                         moieties=DEFAULT_CONTACT_MOIETIES, value="contact", style="scatter",
                         gray=True, gray_cutoff=None,
                         background=CONTACT_MAP_BACKGROUND, color=CONTACT_MAP_COLOR, cmap=CONTACT_MAP_CMAP,
                         size_in=CONTACT_MAP_SIZE_IN, font_size=None,
                         tick_step=CONTACT_MAP_TICK_STEP, n_ticks=CONTACT_MAP_N_TICKS, show_ticks=True,
                         tick_size=None, tick_rotation=0, xlabel=None, ylabel=None,
                         min_marker_pt=CONTACT_MAP_MIN_MARKER_PT,
                         dpi=CONTACT_MAP_DPI, title=None, ax=None, save_path=None, matrix=None):
        """Plot a fixed-aspect contact map for one chain or a chain pair, with each axis titled by its chain ID.

        With value 'contact' and gray True, every residue pair up to gray_cutoff is drawn in light gray,
        shaded by distance (darker when closer). Contacts within cutoff are drawn on top in color.

        Args:
            chain_a: Chain on the y-axis, standard or author ID.
            chain_b: Chain on the x-axis, or None for an intra-chain map.
            cutoff: Distance cutoff in nm for the colored contacts.
            moieties: Biomolecule types to include.
            value: 'contact' or 'distance' (minimum distances in nm).
            style: 'scatter' or 'heatmap'.
            gray: Whether to draw the grayscale background of all contacts up to gray_cutoff.
            gray_cutoff: Upper distance limit of the background in nm. None uses CONTACT_GRAY_FACTOR times cutoff.
            background: Background color.
            color: Contact color when value is 'contact'.
            cmap: Colormap name when value is 'distance'.
            size_in: Figure width and height in inches.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Spacing of residue-number ticks in sequence positions, used when n_ticks is None.
            n_ticks: Number of ticks per axis, spread evenly along each chain. Overrides tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for the quoted chain ID of chain_b.
            ylabel: Y-axis title, or None for the quoted chain ID of chain_a.
            min_marker_pt: Smallest marker side in points for the scatter style.
            dpi: Figure resolution.
            title: Optional title, or None for a default.
            ax: Existing Matplotlib Axes to draw on, or None to create a figure.
            save_path: Optional file path to save the figure to.
            matrix: Optional precomputed pair matrix from contact_matrix, with NaN for non-contacts.

        Returns:
            Tuple of (figure, axes).
        """
        if value not in CONTACT_MAP_VALUES:
            raise ValueError(f"value must be one of {CONTACT_MAP_VALUES}.")
        if style not in CONTACT_MAP_STYLES:
            raise ValueError(f"style must be one of {CONTACT_MAP_STYLES}.")
        font_size = _font_size(size_in, font_size)
        chain_a, chain_b, _, _ = self._pair_keys(chain_a, chain_b)
        name_a, name_b = self._chain_name(chain_a), self._chain_name(chain_b)

        gray_matrix = None
        if gray and value == "contact":
            gray_cutoff = _gray_limit(cutoff, gray_cutoff)
            gray_matrix = self._gray_matrix(chain_a, chain_b, gray_cutoff, moieties)
        if matrix is None:
            if gray_matrix is not None:
                matrix = np.where(gray_matrix <= cutoff, gray_matrix, np.nan)
            else:
                matrix = self.contact_matrix(chain_a, chain_b, cutoff, moieties)[0]
        else:
            matrix = self._checked_matrix(matrix, chain_a, chain_b)

        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi)
        else:
            fig = ax.figure
        self._draw_pair(
            ax, matrix, chain_a, chain_b, value=value, style=style, background=background, color=color,
            cmap=cmap, cutoff=cutoff, font_size=font_size, tick_step=tick_step,
            min_marker_pt=min_marker_pt, map_size_in=size_in,
            title=title or f"Contacts {name_a}-{name_b} (cutoff {cutoff} nm)",
            gray_matrix=gray_matrix, gray_cutoff=gray_cutoff, n_ticks=n_ticks, show_ticks=show_ticks,
            tick_size=tick_size, tick_rotation=tick_rotation, xlabel=xlabel, ylabel=ylabel,
        )
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, ax

    def plot_contacts_two(self, chain_a="A", chain_b="B", cutoff=DEFAULT_CUTOFF,
                          moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca", gray_cutoff=None,
                          inter_color=CONTACT_MAP_COLOR, intra_color=TWO_INTRA_COLOR,
                          background=CONTACT_MAP_BACKGROUND, box_color=TWO_BOX_COLOR, box_lw=TWO_BOX_LW,
                          plddt=False, min_separation=COMPLETE_MIN_SEPARATION,
                          size_in=CONTACT_MAP_SIZE_IN, font_size=None,
                          tick_step=CONTACT_MAP_TICK_STEP, n_ticks=CONTACT_MAP_N_TICKS, show_ticks=True,
                          tick_size=None, tick_rotation=0, xlabel=None, ylabel=None,
                          min_marker_pt=CONTACT_MAP_MIN_MARKER_PT, show_color_legend=False,
                          show_distance_legend=False,
                          dpi=CONTACT_MAP_DPI, title=None, ax=None, save_path=None):
        """Plot two chains on one contact map, with solid contact squares, a dashed chain boundary, and a thin square around each chain.

        Each axis lists chain A's residues followed by chain B's. Every residue pair up to gray_cutoff is
        drawn in grayscale by distance. Within cutoff, intra-chain contacts are drawn in intra_color and
        inter-chain contacts (the interface) in inter_color. If plddt is True, AlphaFold pLDDT bands are
        drawn along the top and left edges, one cell per residue.

        Args:
            chain_a: First chain, standard or author ID. Placed first on both axes.
            chain_b: Second chain, standard or author ID. Placed after chain_a.
            cutoff: Distance cutoff in nm for colored contacts.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            gray_cutoff: Upper distance limit of the grayscale background in nm. None uses CONTACT_GRAY_FACTOR times cutoff.
            inter_color: Color of inter-chain contacts within cutoff.
            intra_color: Color of intra-chain contacts within cutoff.
            background: Background color.
            box_color: Color of the outline and dashed boundary around each chain.
            box_lw: Outline and dashed boundary width in points.
            plddt: Whether to draw pLDDT tracks. Requires AlphaFold models, where pLDDT is in the B-factor column.
            min_separation: Smallest sequence separation kept for intra-chain contacts.
            size_in: Figure width and height in inches.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Spacing of residue-number ticks, used when n_ticks is None.
            n_ticks: Number of ticks per chain, spread evenly along each chain. Overrides tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for the quoted chain IDs.
            ylabel: Y-axis title, or None for the quoted chain IDs.
            min_marker_pt: Smallest marker side in points.
            show_color_legend: Whether to draw the legend of intra-chain, inter-chain, and pLDDT colors.
            show_distance_legend: Whether to draw the grayscale distance colorbar.
            dpi: Figure resolution.
            title: Optional title, or None for a default.
            ax: Existing Matplotlib Axes to draw on, or None to create a figure.
            save_path: Optional file path to save the figure to.

        Returns:
            Tuple of (figure, axes).
        """
        font_size = _font_size(size_in, font_size)
        k1, k2, label1, label2 = self._pair_keys(chain_a, chain_b)
        if k1 == k2:
            raise ValueError("plot_contacts_two needs two different chains.")
        name1, name2 = self._chain_name(k1), self._chain_name(k2)
        gray_cutoff = _gray_limit(cutoff, gray_cutoff)
        n1, n2 = self.chains[k1]["sequence_length"], self.chains[k2]["sequence_length"]
        size = n1 + n2

        contacts = self.compute_contacts(cutoff=gray_cutoff, contact_type="all", chain_ids=[k1, k2],
                                         scheme=scheme, moieties=moieties)
        block = np.full((size, size), np.nan)
        if not contacts.empty:
            gi, gj, distance = self._block_entries(contacts, label1, n1)
            block[gi, gj] = distance
            block[gj, gi] = distance
        index = np.arange(size)
        same_block = (index < n1)[:, None] == (index < n1)[None, :]
        far_enough = np.abs(index[:, None] - index[None, :]) >= min_separation
        close = block <= cutoff
        gray_mask = ~np.isnan(block) & (~same_block | far_enough)
        intra_close = close & same_block & far_enough
        inter_close = close & ~same_block

        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi)
        else:
            fig = ax.figure
        margin = TWO_TRACK_OFFSET + 0.5 if plddt else 0.5
        ax.set_facecolor(background)
        ax.set_xlim(-margin, size - 0.5)
        ax.set_ylim(size - 0.5, -margin)
        ax.set_aspect("equal", adjustable="box")

        marker_pt = _marker_points(size_in, size, size, min_marker_pt)
        marker_area = marker_pt ** 2
        rows, cols = np.nonzero(gray_mask)
        background_scatter = ax.scatter(cols, rows, c=block[rows, cols], cmap=CONTACT_GRAY_CMAP,
                                        vmin=0, vmax=gray_cutoff, s=marker_area, marker="s",
                                        linewidths=0, zorder=1)
        for mask, color in ((intra_close, intra_color), (inter_close, inter_color)):
            rows, cols = np.nonzero(mask)
            ax.scatter(cols, rows, s=marker_area, marker="s", color=color, linewidths=0, zorder=2)

        if plddt:
            values = np.concatenate([self._residue_plddt(k1), self._residue_plddt(k2)])
            colors = plddt_colors(values)
            positions = np.arange(size)
            ax.scatter(positions, np.full(size, -TWO_TRACK_OFFSET), c=colors, s=marker_area, marker="s",
                       linewidths=0, clip_on=False, zorder=2)
            ax.scatter(np.full(size, -TWO_TRACK_OFFSET), positions, c=colors, s=marker_area, marker="s",
                       linewidths=0, clip_on=False, zorder=2)

        boundaries = np.array([n1 - 0.5])
        ax.hlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)
        ax.vlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)

        for start, length in ((0, n1), (n1, n2)):
            ax.add_patch(Rectangle((start - 0.5, start - 0.5), length, length, fill=False,
                                   edgecolor=box_color, linewidth=box_lw, zorder=3))

        p1, labels1 = self._tick_labels(k1, n1, tick_step, n_ticks)
        p2, labels2 = self._tick_labels(k2, n2, tick_step, n_ticks)
        ticks = np.concatenate([p1, p2[1:] + n1])
        tick_labels = labels1 + labels2[1:]
        default_axis = f"{name1} | {name2}"
        self._style_axes(ax, ticks, tick_labels, ticks, tick_labels, show_ticks=show_ticks,
                         tick_size=tick_size, tick_rotation=tick_rotation,
                         xlabel=xlabel or default_axis, ylabel=ylabel or default_axis, font_size=font_size)
        ax.set_title(title or f"Contacts {name1} and {name2} (cutoff {cutoff} nm)",
                     fontsize=font_size * CONTACT_TITLE_RATIO)

        if show_color_legend:
            handles = [Patch(facecolor=intra_color, label="intra-chain"),
                       Patch(facecolor=inter_color, label="inter-chain")]
            if plddt:
                handles += [Patch(facecolor=color, label=f"pLDDT {label}") for label, color in PLDDT_BANDS]
            ax.legend(handles=handles, loc="upper center",
                      bbox_to_anchor=(0.5, _legend_offset(size_in, font_size, len(handles), 2)), ncol=2,
                      fontsize=font_size, frameon=False)
        if show_distance_legend:
            colorbar = fig.colorbar(background_scatter, ax=ax, fraction=0.046, pad=0.04)
            colorbar.set_label("distance (nm)", fontsize=font_size)
            colorbar.ax.tick_params(labelsize=font_size * CONTACT_TICK_RATIO)
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, ax

    def plot_contacts_three(self, chain_a="A", chain_b="B", chain_c="C", cutoff=DEFAULT_CUTOFF,
                            moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca", gray_cutoff=None,
                            pair_colors=THREE_PAIR_COLORS, intra_color=TWO_INTRA_COLOR,
                            background=CONTACT_MAP_BACKGROUND, box_color=TWO_BOX_COLOR, box_lw=TWO_BOX_LW,
                            plddt=False, min_separation=COMPLETE_MIN_SEPARATION,
                            size_in=CONTACT_MAP_SIZE_IN, font_size=None,
                            tick_step=CONTACT_MAP_TICK_STEP, n_ticks=CONTACT_MAP_N_TICKS, show_ticks=True,
                            tick_size=None, tick_rotation=0, xlabel=None, ylabel=None,
                            min_marker_pt=CONTACT_MAP_MIN_MARKER_PT, show_color_legend=False,
                            show_distance_legend=False,
                            dpi=CONTACT_MAP_DPI, title=None, ax=None, save_path=None):
        """Plot three chains on one contact map, with solid contact squares, dashed chain boundaries, and a thin square around each chain.

        Each axis lists chain A's residues, then chain B's, then chain C's. Every residue pair up to
        gray_cutoff is drawn in grayscale by distance. Within cutoff, intra-chain contacts are drawn in
        intra_color, and inter-chain contacts are drawn in one color per chain pair (A-B, A-C, B-C), so
        each interface can be told apart. If plddt is True, AlphaFold pLDDT bands are drawn along the top
        and left edges, one cell per residue.

        Args:
            chain_a: First chain, standard or author ID. Placed first on both axes.
            chain_b: Second chain, standard or author ID. Placed after chain_a.
            chain_c: Third chain, standard or author ID. Placed after chain_b.
            cutoff: Distance cutoff in nm for colored contacts.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            gray_cutoff: Upper distance limit of the grayscale background in nm. None uses CONTACT_GRAY_FACTOR times cutoff.
            pair_colors: Colors of the inter-chain contacts for the pairs (A-B, A-C, B-C).
            intra_color: Color of intra-chain contacts within cutoff.
            background: Background color.
            box_color: Color of the outline and dashed boundaries around each chain.
            box_lw: Outline and dashed boundary width in points.
            plddt: Whether to draw pLDDT tracks. Requires AlphaFold models, where pLDDT is in the B-factor column.
            min_separation: Smallest sequence separation kept for intra-chain contacts.
            size_in: Figure width and height in inches.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Spacing of residue-number ticks, used when n_ticks is None.
            n_ticks: Number of ticks per chain, spread evenly along each chain. Overrides tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for the quoted chain IDs.
            ylabel: Y-axis title, or None for the quoted chain IDs.
            min_marker_pt: Smallest marker side in points.
            show_color_legend: Whether to draw the legend of intra-chain, inter-chain, and pLDDT colors.
            show_distance_legend: Whether to draw the grayscale distance colorbar.
            dpi: Figure resolution.
            title: Optional title, or None for a default.
            ax: Existing Matplotlib Axes to draw on, or None to create a figure.
            save_path: Optional file path to save the figure to.

        Returns:
            Tuple of (figure, axes).
        """
        font_size = _font_size(size_in, font_size)
        keys = [self._require_chain(c) for c in (chain_a, chain_b, chain_c)]
        if len(set(keys)) != 3:
            raise ValueError("plot_contacts_three needs three different chains.")
        labels = [self.chains[k]["chain_id"] for k in keys]
        names = [self._chain_name(k) for k in keys]
        sizes = [self.chains[k]["sequence_length"] for k in keys]
        offsets = np.cumsum([0] + sizes[:-1])
        gray_cutoff = _gray_limit(cutoff, gray_cutoff)
        size = sum(sizes)

        contacts = self.compute_contacts(cutoff=gray_cutoff, contact_type="all", chain_ids=keys,
                                         scheme=scheme, moieties=moieties)
        block = np.full((size, size), np.nan)
        if not contacts.empty:
            gi, gj, distance = self._multi_block_entries(contacts, labels, offsets)
            block[gi, gj] = distance
            block[gj, gi] = distance
        block_of = np.repeat(np.arange(3), sizes)
        pair_code = np.full((3, 3), -1)
        for p, (a, b) in enumerate(THREE_PAIRS):
            pair_code[a, b] = pair_code[b, a] = p
        index = np.arange(size)
        same_block = block_of[:, None] == block_of[None, :]
        pair_index = pair_code[block_of[:, None], block_of[None, :]]
        far_enough = np.abs(index[:, None] - index[None, :]) >= min_separation
        close = block <= cutoff
        gray_mask = ~np.isnan(block) & (~same_block | far_enough)
        intra_close = close & same_block & far_enough
        inter_close = close & ~same_block

        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi)
        else:
            fig = ax.figure
        margin = TWO_TRACK_OFFSET + 0.5 if plddt else 0.5
        ax.set_facecolor(background)
        ax.set_xlim(-margin, size - 0.5)
        ax.set_ylim(size - 0.5, -margin)
        ax.set_aspect("equal", adjustable="box")

        marker_pt = _marker_points(size_in, size, size, min_marker_pt)
        marker_area = marker_pt ** 2
        rows, cols = np.nonzero(gray_mask)
        background_scatter = ax.scatter(cols, rows, c=block[rows, cols], cmap=CONTACT_GRAY_CMAP,
                                        vmin=0, vmax=gray_cutoff, s=marker_area, marker="s",
                                        linewidths=0, zorder=1)
        rows, cols = np.nonzero(intra_close)
        ax.scatter(cols, rows, s=marker_area, marker="s", color=intra_color, linewidths=0, zorder=2)
        for p, color in enumerate(pair_colors):
            rows, cols = np.nonzero(inter_close & (pair_index == p))
            ax.scatter(cols, rows, s=marker_area, marker="s", color=color, linewidths=0, zorder=2)

        if plddt:
            values = np.concatenate([self._residue_plddt(k) for k in keys])
            colors = plddt_colors(values)
            positions = np.arange(size)
            ax.scatter(positions, np.full(size, -TWO_TRACK_OFFSET), c=colors, s=marker_area, marker="s",
                       linewidths=0, clip_on=False, zorder=2)
            ax.scatter(np.full(size, -TWO_TRACK_OFFSET), positions, c=colors, s=marker_area, marker="s",
                       linewidths=0, clip_on=False, zorder=2)

        boundaries = offsets[1:] - 0.5
        ax.hlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)
        ax.vlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)

        for start, length in zip(offsets, sizes):
            ax.add_patch(Rectangle((start - 0.5, start - 0.5), length, length, fill=False,
                                   edgecolor=box_color, linewidth=box_lw, zorder=3))

        ticks, tick_labels = [], []
        for index, (k, offset) in enumerate(zip(keys, offsets)):
            positions, chain_labels = self._tick_labels(k, self.chains[k]["sequence_length"], tick_step, n_ticks)
            if index:
                positions, chain_labels = positions[1:], chain_labels[1:]
            ticks.append(positions + offset)
            tick_labels += chain_labels
        ticks = np.concatenate(ticks)
        default_axis = " | ".join(names)
        self._style_axes(ax, ticks, tick_labels, ticks, tick_labels, show_ticks=show_ticks,
                         tick_size=tick_size, tick_rotation=tick_rotation,
                         xlabel=xlabel or default_axis, ylabel=ylabel or default_axis, font_size=font_size)
        ax.set_title(title or f"Contacts {'-'.join(names)} (cutoff {cutoff} nm)",
                     fontsize=font_size * CONTACT_TITLE_RATIO)

        if show_color_legend:
            handles = [Patch(facecolor=intra_color, label="intra-chain")]
            handles += [Patch(facecolor=pair_colors[p], label=f"{names[a]}-{names[b]}")
                        for p, (a, b) in enumerate(THREE_PAIRS)]
            if plddt:
                handles += [Patch(facecolor=color, label=f"pLDDT {label}") for label, color in PLDDT_BANDS]
            ax.legend(handles=handles, loc="upper center",
                      bbox_to_anchor=(0.5, _legend_offset(size_in, font_size, len(handles), 2)), ncol=2,
                      fontsize=font_size, frameon=False)
        if show_distance_legend:
            colorbar = fig.colorbar(background_scatter, ax=ax, fraction=0.046, pad=0.04)
            colorbar.set_label("distance (nm)", fontsize=font_size)
            colorbar.ax.tick_params(labelsize=font_size * CONTACT_TICK_RATIO)
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, ax

    def plot_contacts_four(self, chain_a="A", chain_b="B", chain_c="C", chain_d="D", cutoff=DEFAULT_CUTOFF,
                           moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca", gray_cutoff=None,
                           pair_colors=FOUR_PAIR_COLORS, intra_color=TWO_INTRA_COLOR,
                           background=CONTACT_MAP_BACKGROUND, box_color=TWO_BOX_COLOR, box_lw=TWO_BOX_LW,
                           plddt=False, min_separation=COMPLETE_MIN_SEPARATION,
                           size_in=CONTACT_FOUR_SIZE_IN, font_size=None,
                           tick_step=CONTACT_MAP_TICK_STEP, n_ticks=CONTACT_MAP_N_TICKS, show_ticks=True,
                           tick_size=None, tick_rotation=0, xlabel=None, ylabel=None,
                           min_marker_pt=CONTACT_MAP_MIN_MARKER_PT, show_color_legend=False,
                           show_distance_legend=False,
                           dpi=CONTACT_MAP_DPI, title=None, ax=None, save_path=None):
        """Plot four chains on one contact map, with solid contact squares, dashed chain boundaries, and a thin square around each chain.

        Args:
            chain_a: First chain, standard or author ID. Placed first on both axes.
            chain_b: Second chain, standard or author ID. Placed after chain_a.
            chain_c: Third chain, standard or author ID. Placed after chain_b.
            chain_d: Fourth chain, standard or author ID. Placed after chain_c.
            cutoff: Distance cutoff in nm for colored contacts.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            gray_cutoff: Upper distance limit of the grayscale background in nm. None uses CONTACT_GRAY_FACTOR times cutoff.
            pair_colors: Colors of the inter-chain contacts for the pairs in FOUR_PAIRS order.
            intra_color: Color of intra-chain contacts within cutoff.
            background: Background color.
            box_color: Color of the outline and dashed boundaries around each chain.
            box_lw: Outline and dashed boundary width in points.
            plddt: Whether to draw pLDDT tracks. Requires AlphaFold models, where pLDDT is in the B-factor column.
            min_separation: Smallest sequence separation kept for intra-chain contacts.
            size_in: Figure width and height in inches.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Spacing of residue-number ticks, used when n_ticks is None.
            n_ticks: Number of ticks per chain, spread evenly along each chain. Overrides tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for the quoted chain IDs.
            ylabel: Y-axis title, or None for the quoted chain IDs.
            min_marker_pt: Smallest marker side in points.
            show_color_legend: Whether to draw the legend of intra-chain, inter-chain, and pLDDT colors.
            show_distance_legend: Whether to draw the grayscale distance colorbar.
            dpi: Figure resolution.
            title: Optional title, or None for a default.
            ax: Existing Matplotlib Axes to draw on, or None to create a figure.
            save_path: Optional file path to save the figure to.

        Returns:
            Tuple of (figure, axes).
        """
        font_size = _font_size(size_in, font_size)
        keys = [self._require_chain(c) for c in (chain_a, chain_b, chain_c, chain_d)]
        if len(set(keys)) != 4:
            raise ValueError("plot_contacts_four needs four different chains.")
        labels = [self.chains[k]["chain_id"] for k in keys]
        names = [self._chain_name(k) for k in keys]
        sizes = [self.chains[k]["sequence_length"] for k in keys]
        offsets = np.cumsum([0] + sizes[:-1])
        gray_cutoff = _gray_limit(cutoff, gray_cutoff)
        size = sum(sizes)

        contacts = self.compute_contacts(cutoff=gray_cutoff, contact_type="all", chain_ids=keys,
                                         scheme=scheme, moieties=moieties)
        block = np.full((size, size), np.nan)
        if not contacts.empty:
            gi, gj, distance = self._multi_block_entries(contacts, labels, offsets)
            block[gi, gj] = distance
            block[gj, gi] = distance
        block_of = np.repeat(np.arange(4), sizes)
        pair_code = np.full((4, 4), -1)
        for p, (a, b) in enumerate(FOUR_PAIRS):
            pair_code[a, b] = pair_code[b, a] = p
        index = np.arange(size)
        same_block = block_of[:, None] == block_of[None, :]
        pair_index = pair_code[block_of[:, None], block_of[None, :]]
        far_enough = np.abs(index[:, None] - index[None, :]) >= min_separation
        close = block <= cutoff
        gray_mask = ~np.isnan(block) & (~same_block | far_enough)
        intra_close = close & same_block & far_enough
        inter_close = close & ~same_block

        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi)
        else:
            fig = ax.figure
        margin = TWO_TRACK_OFFSET + 0.5 if plddt else 0.5
        ax.set_facecolor(background)
        ax.set_xlim(-margin, size - 0.5)
        ax.set_ylim(size - 0.5, -margin)
        ax.set_aspect("equal", adjustable="box")

        marker_pt = _marker_points(size_in, size, size, min_marker_pt)
        marker_area = marker_pt ** 2
        rows, cols = np.nonzero(gray_mask)
        background_scatter = ax.scatter(cols, rows, c=block[rows, cols], cmap=CONTACT_GRAY_CMAP,
                                        vmin=0, vmax=gray_cutoff, s=marker_area, marker="s",
                                        linewidths=0, zorder=1)
        rows, cols = np.nonzero(intra_close)
        ax.scatter(cols, rows, s=marker_area, marker="s", color=intra_color, linewidths=0, zorder=2)
        for p, color in enumerate(pair_colors):
            rows, cols = np.nonzero(inter_close & (pair_index == p))
            ax.scatter(cols, rows, s=marker_area, marker="s", color=color, linewidths=0, zorder=2)

        if plddt:
            values = np.concatenate([self._residue_plddt(k) for k in keys])
            colors = plddt_colors(values)
            positions = np.arange(size)
            ax.scatter(positions, np.full(size, -TWO_TRACK_OFFSET), c=colors, s=marker_area, marker="s",
                       linewidths=0, clip_on=False, zorder=2)
            ax.scatter(np.full(size, -TWO_TRACK_OFFSET), positions, c=colors, s=marker_area, marker="s",
                       linewidths=0, clip_on=False, zorder=2)

        boundaries = offsets[1:] - 0.5
        ax.hlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)
        ax.vlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)

        for start, length in zip(offsets, sizes):
            ax.add_patch(Rectangle((start - 0.5, start - 0.5), length, length, fill=False,
                                   edgecolor=box_color, linewidth=box_lw, zorder=3))

        ticks, tick_labels = [], []
        for k_index, (k, offset) in enumerate(zip(keys, offsets)):
            positions, chain_labels = self._tick_labels(k, self.chains[k]["sequence_length"], tick_step, n_ticks)
            if k_index:
                positions, chain_labels = positions[1:], chain_labels[1:]
            ticks.append(positions + offset)
            tick_labels += chain_labels
        ticks = np.concatenate(ticks)
        default_axis = " | ".join(names)
        self._style_axes(ax, ticks, tick_labels, ticks, tick_labels, show_ticks=show_ticks,
                         tick_size=tick_size, tick_rotation=tick_rotation,
                         xlabel=xlabel or default_axis, ylabel=ylabel or default_axis, font_size=font_size)
        ax.set_title(title or f"Contacts {'-'.join(names)} (cutoff {cutoff} nm)",
                     fontsize=font_size * CONTACT_TITLE_RATIO)

        if show_color_legend:
            handles = [Patch(facecolor=intra_color, label="intra-chain")]
            handles += [Patch(facecolor=pair_colors[p], label=f"{names[a]}-{names[b]}")
                        for p, (a, b) in enumerate(FOUR_PAIRS)]
            if plddt:
                handles += [Patch(facecolor=color, label=f"pLDDT {label}") for label, color in PLDDT_BANDS]
            ax.legend(handles=handles, loc="upper center",
                      bbox_to_anchor=(0.5, _legend_offset(size_in, font_size, len(handles), 2)), ncol=2,
                      fontsize=font_size, frameon=False)
        if show_distance_legend:
            colorbar = fig.colorbar(background_scatter, ax=ax, fraction=0.046, pad=0.04)
            colorbar.set_label("distance (nm)", fontsize=font_size)
            colorbar.ax.tick_params(labelsize=font_size * CONTACT_TICK_RATIO)
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, ax

    def plot_contact_panel(self, chain_a="A", chain_b=None, cutoff=DEFAULT_CUTOFF,
                           moieties=DEFAULT_CONTACT_MOIETIES, style="scatter",
                           background=CONTACT_MAP_BACKGROUND, color=CONTACT_MAP_COLOR,
                           size_in=CONTACT_PANEL_SIZE_IN, font_size=None,
                           tick_step=CONTACT_MAP_TICK_STEP, n_ticks=CONTACT_MAP_N_TICKS, show_ticks=True,
                           tick_size=None, min_marker_pt=CONTACT_MAP_MIN_MARKER_PT,
                           dpi=CONTACT_MAP_DPI, save_path=None):
        """Plot a 2x2 summary: contact map, distance histogram, per-residue profile, and sequence separation.

        Args:
            chain_a: First chain, standard or author ID.
            chain_b: Second chain, or None for an intra-chain summary.
            cutoff: Distance cutoff in nm.
            moieties: Biomolecule types to include.
            style: Map style, 'scatter' or 'heatmap'.
            background: Background color for the map.
            color: Contact color for the map. Profile and histogram lines use the palette.
            size_in: Figure width and height in inches.
            font_size: Base font size in points, or None to scale with each panel's size.
            tick_step: Spacing of residue-number ticks on the map, used when n_ticks is None.
            n_ticks: Number of ticks per axis on the map. Overrides tick_step.
            show_ticks: Whether to draw tick marks and labels on the map.
            tick_size: Tick label font size on the map, or None for a fixed share of font_size.
            min_marker_pt: Smallest marker side in points.
            dpi: Figure resolution.
            save_path: Optional file path to save the figure to.

        Returns:
            Tuple of (figure, axes) where axes is a 2x2 array.
        """
        font_size = _font_size(size_in / 2, font_size)
        chain_a, chain_b, label_a, label_b = self._pair_keys(chain_a, chain_b)
        name_a, name_b = self._chain_name(chain_a), self._chain_name(chain_b)
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all",
                                         chain_ids=sorted({chain_a, chain_b}), moieties=moieties)
        matrix = self._matrix_from_contacts(contacts, chain_a, chain_b)
        keys = [chain_a] if chain_a == chain_b else [chain_a, chain_b]
        labels = [self.chains[k]["chain_id"] for k in keys]
        names = [self._chain_name(k) for k in keys]

        fig, axes = plt.subplots(2, 2, figsize=(size_in, size_in), dpi=dpi)
        ax_map, ax_dist, ax_profile, ax_sep = axes.ravel()

        self._draw_pair(
            ax_map, matrix, chain_a, chain_b, value="contact", style=style, background=background,
            color=color, cmap=CONTACT_MAP_CMAP, cutoff=cutoff, font_size=font_size, tick_step=tick_step,
            min_marker_pt=min_marker_pt, map_size_in=size_in * 0.5,
            title=f"Contacts {name_a}-{name_b}", n_ticks=n_ticks, show_ticks=show_ticks, tick_size=tick_size,
        )

        distances = self._pair_rows(contacts, label_a, label_b)["min_distance_nm"].to_numpy()
        ax_dist.hist(distances, bins=CONTACT_HIST_BINS, range=(0, cutoff), color=color)
        ax_dist.set_xlabel("min distance (nm)", fontsize=font_size)
        ax_dist.set_ylabel("contacts", fontsize=font_size)
        ax_dist.set_title("Contact distances", fontsize=font_size * CONTACT_TITLE_RATIO)

        profile = self._residue_profile(contacts, moieties)
        separation = self._separation_frame(contacts)
        for index, (label, name) in enumerate(zip(labels, names)):
            line_color = CONTACT_PANEL_LINE_COLORS[index % len(CONTACT_PANEL_LINE_COLORS)]
            rows = profile[profile["chain_id"] == label]
            ax_profile.step(rows["resnum"], rows["n_contacts"], where="mid", color=line_color, label=f"chain {name}")
            values = separation.loc[separation["chain_id"] == label, "separation"]
            ax_sep.hist(values, bins=CONTACT_HIST_BINS, histtype="step", color=line_color, label=f"chain {name}")
        ax_profile.set_xlabel("residue number", fontsize=font_size)
        ax_profile.set_ylabel("contacts per residue", fontsize=font_size)
        ax_profile.set_title("Per-residue contacts", fontsize=font_size * CONTACT_TITLE_RATIO)
        ax_sep.set_xlabel("sequence separation |i - j|", fontsize=font_size)
        ax_sep.set_ylabel("intra-chain contacts", fontsize=font_size)
        ax_sep.set_title("Sequence separation", fontsize=font_size * CONTACT_TITLE_RATIO)
        if len(labels) > 1:
            ax_profile.legend(fontsize=font_size, frameon=False)
            ax_sep.legend(fontsize=font_size, frameon=False)

        for ax in (ax_dist, ax_profile, ax_sep):
            ax.tick_params(labelsize=font_size * CONTACT_TICK_RATIO)
        fig.suptitle(f"Contact summary {name_a}-{name_b} (cutoff {cutoff} nm)",
                     fontsize=font_size * CONTACT_SUPTITLE_RATIO)
        fig.tight_layout()
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, axes