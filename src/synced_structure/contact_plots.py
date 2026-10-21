"""Contact map plots: pair maps with grayscale background contacts, multi-chain block maps with chain outlines and optional pLDDT tracks, and a summary panel."""

###################################### IMPORTS AND CONSTANTS

import itertools
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle

from .constants import COMPLETE_MIN_SEPARATION, CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES
from .plot_style import (
    CONTACT_AXIS_TITLE_RATIO, CONTACT_FONT_MAX_PT, CONTACT_FONT_MIN_PT, CONTACT_FONT_PT_PER_IN, CONTACT_GRAY_CMAP,
    CONTACT_GRAY_FACTOR, CONTACT_HIST_BINS, CONTACT_LEGEND_ROW_RATIO, CONTACT_MAP_BACKGROUND, CONTACT_MAP_CMAP,
    CONTACT_MAP_COLOR, CONTACT_MAP_DPI, CONTACT_MAP_MIN_MARKER_PT, CONTACT_MAP_N_TICKS, CONTACT_MAP_SIZE_IN,
    CONTACT_MAP_STYLES, CONTACT_MAP_TICK_STEP, CONTACT_MAP_VALUES, CONTACT_MULTI_SIZE_IN, CONTACT_PANEL_LINE_COLORS,
    CONTACT_PANEL_SIZE_IN, CONTACT_SUPTITLE_RATIO, CONTACT_TICK_BLOCK_RATIO, CONTACT_TICK_RATIO, CONTACT_TITLE_RATIO,
    MULTI_BOX_COLOR, MULTI_BOX_LW, MULTI_INTRA_COLOR, MULTI_PAIR_COLORS, MULTI_TRACK_OFFSET, PLDDT_BANDS,
    PLDDT_MISSING_COLOR,
)

###################################### LAYOUT HELPERS


def _font_size(size_in, font_size=None):
    """Return the base font size in points, proportional to the plot size unless given explicitly."""
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
    """Return the marker side in points so that one residue is at least min_marker_pt wide."""
    return max(min_marker_pt, map_size_in * 72 * 0.8 / max(n_cols, n_rows, 1))


def _gray_limit(cutoff, gray_cutoff):
    """Return the upper distance in nm of the grayscale background: gray_cutoff, or CONTACT_GRAY_FACTOR times cutoff."""
    if gray_cutoff is None:
        return cutoff * CONTACT_GRAY_FACTOR
    return max(gray_cutoff, cutoff)


def _tick_positions(n, tick_step, n_ticks):
    """Return tick positions along an axis of n residues: n_ticks evenly spread, or every tick_step residues.

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
    """Return AlphaFold pLDDT band colors for per-residue values in 0 to 100, gray where a value is NaN."""
    values = np.asarray(values, dtype=float)
    colors = np.where(values >= 90, PLDDT_BANDS[0][1],
                      np.where(values >= 70, PLDDT_BANDS[1][1],
                               np.where(values >= 50, PLDDT_BANDS[2][1], PLDDT_BANDS[3][1])))
    return [PLDDT_MISSING_COLOR if np.isnan(v) else c for v, c in zip(values, colors)]

###################################### DRAWING HELPERS


class ContactPlotMixin:
    def _residue_plddt(self, chain_key):
        """Return per-residue pLDDT for a chain from the C-alpha B-factor (atom mean otherwise), NaN where missing.

        AlphaFold models store pLDDT in the B-factor column; for experimental structures these are B-factors.
        """
        values = np.full(self.chains[chain_key]["sequence_length"], np.nan)
        for position, residue in enumerate(self.chains[chain_key]["residues"]):
            if "CA" in residue:
                values[position] = residue["CA"].get_bfactor()
            else:
                atoms = list(residue.get_atoms())
                if atoms:
                    values[position] = np.mean([atom.get_bfactor() for atom in atoms])
        return values

    def _chain_name(self, chain_key):
        """Return the quoted author chain ID used in plot text, such as 'A'."""
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
        return positions, [str(residues[p].id[1]) if p < len(residues) else str(p + 1) for p in positions]

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
            xlabel, ylabel: Axis titles.
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
        """Draw a pair matrix on ax over an optional grayscale background, and return a mappable or None.

        Distance values use a fixed scale from 0 to cutoff, so plots with the same cutoff are comparable.

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
            return ax.imshow(np.ma.masked_invalid(matrix), cmap=colormap, vmin=0, vmax=cutoff, interpolation="nearest")
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
                         tick_rotation=tick_rotation, xlabel=xlabel or self._chain_name(chain_b),
                         ylabel=ylabel or self._chain_name(chain_a), font_size=font_size)
        ax.set_title(title, fontsize=font_size * CONTACT_TITLE_RATIO)

    def _checked_matrix(self, matrix, chain_a, chain_b):
        """Return a precomputed matrix as float after checking its shape matches the two chains."""
        matrix = np.asarray(matrix, dtype=float)
        expected = (self.chains[chain_a]["sequence_length"], self.chains[chain_b]["sequence_length"])
        if matrix.shape != expected:
            raise ValueError(f"matrix shape {matrix.shape} does not match chains {expected}.")
        return matrix

###################################### PUBLIC PLOTS

    def plot_contact_map(self, chain_a="A", chain_b=None, cutoff=CONTACTS_CUTOFF_NM,
                         moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca", value="contact", style="scatter",
                         gray=True, gray_cutoff=None,
                         background=CONTACT_MAP_BACKGROUND, color=CONTACT_MAP_COLOR, cmap=CONTACT_MAP_CMAP,
                         size_in=CONTACT_MAP_SIZE_IN, font_size=None,
                         tick_step=CONTACT_MAP_TICK_STEP, n_ticks=CONTACT_MAP_N_TICKS, show_ticks=True,
                         tick_size=None, tick_rotation=0, xlabel=None, ylabel=None,
                         min_marker_pt=CONTACT_MAP_MIN_MARKER_PT,
                         dpi=CONTACT_MAP_DPI, title=None, ax=None, save_path=None, matrix=None):
        """Plot a fixed-aspect contact map for one chain or a chain pair, with each axis titled by its chain ID.

        With value 'contact' and gray True, every residue pair up to gray_cutoff is drawn in light gray,
        darker when closer, and contacts within cutoff are drawn on top in color.

        Args:
            chain_a: Chain on the y-axis, standard or author ID.
            chain_b: Chain on the x-axis, or None for an intra-chain map.
            cutoff: Distance cutoff in nm for the colored contacts.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            value: 'contact' or 'distance' (minimum distances in nm).
            style: 'scatter' or 'heatmap'.
            gray: Whether to draw the grayscale background of all contacts up to gray_cutoff.
            gray_cutoff: Upper distance limit of the background in nm. None uses CONTACT_GRAY_FACTOR times cutoff.
            background: Background color.
            color: Contact color when value is 'contact'.
            cmap: Colormap name when value is 'distance'.
            size_in: Figure width and height in inches.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Spacing of residue-number ticks, used when n_ticks is None.
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
        gray_matrix = None
        if gray and value == "contact":
            gray_cutoff = _gray_limit(cutoff, gray_cutoff)
            gray_matrix = self.contact_matrix(chain_a, chain_b, gray_cutoff, moieties, scheme)["matrix"]
        if matrix is not None:
            matrix = self._checked_matrix(matrix, chain_a, chain_b)
        elif gray_matrix is not None:
            matrix = np.where(gray_matrix <= cutoff, gray_matrix, np.nan)
        else:
            matrix = self.contact_matrix(chain_a, chain_b, cutoff, moieties, scheme)["matrix"]
        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi)
        else:
            fig = ax.figure
        self._draw_pair(
            ax, matrix, chain_a, chain_b, value=value, style=style, background=background, color=color,
            cmap=cmap, cutoff=cutoff, font_size=font_size, tick_step=tick_step,
            min_marker_pt=min_marker_pt, map_size_in=size_in,
            title=title or f"Contacts {self._chain_name(chain_a)}-{self._chain_name(chain_b)} (cutoff {cutoff} nm)",
            gray_matrix=gray_matrix, gray_cutoff=gray_cutoff, n_ticks=n_ticks, show_ticks=show_ticks,
            tick_size=tick_size, tick_rotation=tick_rotation, xlabel=xlabel, ylabel=ylabel,
        )
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, ax

    def plot_contacts_multi(self, chain_ids=("A", "B"), cutoff=CONTACTS_CUTOFF_NM,
                            moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca", gray_cutoff=None,
                            pair_colors=MULTI_PAIR_COLORS, intra_color=MULTI_INTRA_COLOR,
                            background=CONTACT_MAP_BACKGROUND, box_color=MULTI_BOX_COLOR, box_lw=MULTI_BOX_LW,
                            plddt=False, min_separation=COMPLETE_MIN_SEPARATION,
                            size_in=None, font_size=None,
                            tick_step=CONTACT_MAP_TICK_STEP, n_ticks=CONTACT_MAP_N_TICKS, show_ticks=True,
                            tick_size=None, tick_rotation=0, xlabel=None, ylabel=None,
                            min_marker_pt=CONTACT_MAP_MIN_MARKER_PT, show_color_legend=False,
                            show_distance_legend=False,
                            dpi=CONTACT_MAP_DPI, title=None, ax=None, save_path=None):
        """Plot several chains on one block contact map, with dashed chain boundaries and a thin square around each chain.

        Each axis lists the chains' residues in order. Every residue pair up to gray_cutoff is drawn in
        grayscale by distance. Within cutoff, intra-chain contacts are drawn in intra_color and
        inter-chain contacts in one color per chain pair, in itertools.combinations order. If plddt is
        True, AlphaFold pLDDT bands are drawn along the top and left edges.

        Args:
            chain_ids: Two or more different chains, standard or author IDs, in axis order.
            cutoff: Distance cutoff in nm for colored contacts.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            gray_cutoff: Upper distance limit of the grayscale background in nm. None uses CONTACT_GRAY_FACTOR times cutoff.
            pair_colors: Inter-chain contact colors, one per chain pair.
            intra_color: Color of intra-chain contacts within cutoff.
            background: Background color.
            box_color: Color of the outlines and dashed boundaries.
            box_lw: Outline and dashed boundary width in points.
            plddt: Whether to draw pLDDT tracks, read from the B-factor column.
            min_separation: Smallest sequence separation kept for intra-chain contacts.
            size_in: Figure width and height in inches, or None for a size by number of chains.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Spacing of residue-number ticks, used when n_ticks is None.
            n_ticks: Number of ticks per chain. Overrides tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for the quoted chain IDs.
            ylabel: Y-axis title, or None for the quoted chain IDs.
            min_marker_pt: Smallest marker side in points.
            show_color_legend: Whether to draw the legend of intra-chain, chain-pair, and pLDDT colors.
            show_distance_legend: Whether to draw the grayscale distance colorbar.
            dpi: Figure resolution.
            title: Optional title, or None for a default.
            ax: Existing Matplotlib Axes to draw on, or None to create a figure.
            save_path: Optional file path to save the figure to.

        Returns:
            Tuple of (figure, axes).
        """
        keys = [self._require_chain(c) for c in chain_ids]
        if len(keys) < 2 or len(set(keys)) != len(keys):
            raise ValueError("plot_contacts_multi needs two or more different chains.")
        pairs = list(itertools.combinations(range(len(keys)), 2))
        if len(pair_colors) < len(pairs):
            raise ValueError(f"pair_colors needs {len(pairs)} colors for {len(keys)} chains.")
        size_in = size_in or CONTACT_MULTI_SIZE_IN.get(len(keys), max(CONTACT_MULTI_SIZE_IN.values()))
        font_size = _font_size(size_in, font_size)
        names = [self._chain_name(k) for k in keys]
        gray_cutoff = _gray_limit(cutoff, gray_cutoff)
        block = self._block_matrix(keys, gray_cutoff, scheme, moieties)
        matrix, sizes, offsets = block["matrix"], block["sizes"], block["offsets"]
        size = sum(sizes)

        block_of = np.repeat(np.arange(len(keys)), sizes)
        pair_code = np.full((len(keys), len(keys)), -1)
        for p, (a, b) in enumerate(pairs):
            pair_code[a, b] = pair_code[b, a] = p
        index = np.arange(size)
        same_block = block_of[:, None] == block_of[None, :]
        pair_index = pair_code[block_of[:, None], block_of[None, :]]
        far_enough = np.abs(index[:, None] - index[None, :]) >= min_separation
        close = matrix <= cutoff
        gray_mask = ~np.isnan(matrix) & (~same_block | far_enough)
        intra_close = close & same_block & far_enough
        inter_close = close & ~same_block

        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi)
        else:
            fig = ax.figure
        margin = MULTI_TRACK_OFFSET + 0.5 if plddt else 0.5
        ax.set_facecolor(background)
        ax.set_xlim(-margin, size - 0.5)
        ax.set_ylim(size - 0.5, -margin)
        ax.set_aspect("equal", adjustable="box")

        marker_area = _marker_points(size_in, size, size, min_marker_pt) ** 2
        rows, cols = np.nonzero(gray_mask)
        background_scatter = ax.scatter(cols, rows, c=matrix[rows, cols], cmap=CONTACT_GRAY_CMAP,
                                        vmin=0, vmax=gray_cutoff, s=marker_area, marker="s", linewidths=0, zorder=1)
        masks = [(intra_close, intra_color)] + [(inter_close & (pair_index == p), pair_colors[p]) for p in range(len(pairs))]
        for mask, color in masks:
            rows, cols = np.nonzero(mask)
            ax.scatter(cols, rows, s=marker_area, marker="s", color=color, linewidths=0, zorder=2)

        if plddt:
            colors = plddt_colors(np.concatenate([self._residue_plddt(k) for k in keys]))
            track = np.full(size, -MULTI_TRACK_OFFSET)
            for x, y in ((index, track), (track, index)):
                ax.scatter(x, y, c=colors, s=marker_area, marker="s", linewidths=0, clip_on=False, zorder=2)

        boundaries = offsets[1:] - 0.5
        ax.hlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)
        ax.vlines(boundaries, -0.5, size - 0.5, colors=box_color, linestyles="--", linewidth=box_lw, zorder=3)
        for start, length in zip(offsets, sizes):
            ax.add_patch(Rectangle((start - 0.5, start - 0.5), length, length, fill=False,
                                   edgecolor=box_color, linewidth=box_lw, zorder=3))

        ticks, tick_labels = [], []
        for k_index, (k, offset, n) in enumerate(zip(keys, offsets, sizes)):
            positions, chain_labels = self._tick_labels(k, n, tick_step, n_ticks)
            if k_index:
                positions, chain_labels = positions[1:], chain_labels[1:]
            ticks.append(positions + offset)
            tick_labels += chain_labels
        ticks = np.concatenate(ticks)
        default_axis = " | ".join(names)
        self._style_axes(ax, ticks, tick_labels, ticks, tick_labels, show_ticks=show_ticks,
                         tick_size=tick_size, tick_rotation=tick_rotation,
                         xlabel=xlabel or default_axis, ylabel=ylabel or default_axis, font_size=font_size)
        ax.set_title(title or f"Contacts {'-'.join(names)} (cutoff {cutoff} nm)", fontsize=font_size * CONTACT_TITLE_RATIO)

        if show_color_legend:
            handles = [Patch(facecolor=intra_color, label="intra-chain")]
            handles += [Patch(facecolor=pair_colors[p], label=f"{names[a]}-{names[b]}") for p, (a, b) in enumerate(pairs)]
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

    def plot_contact_panel(self, chain_a="A", chain_b=None, cutoff=CONTACTS_CUTOFF_NM,
                           moieties=DEFAULT_CONTACT_MOIETIES, scheme="ca", style="scatter",
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
            scheme: Contact scheme, 'ca' or 'sidechain'.
            style: Map style, 'scatter' or 'heatmap'.
            background: Background color for the map.
            color: Contact color for the map and histogram.
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
        pair = self.contact_matrix(chain_a, chain_b, cutoff, moieties, scheme)
        contacts = pair["contacts"]
        keys = list(dict.fromkeys((chain_a, chain_b)))
        stats = self._statistics_from_contacts(contacts, moieties, chain_keys=keys)

        fig, axes = plt.subplots(2, 2, figsize=(size_in, size_in), dpi=dpi)
        ax_map, ax_dist, ax_profile, ax_sep = axes.ravel()
        self._draw_pair(
            ax_map, pair["matrix"], chain_a, chain_b, value="contact", style=style, background=background,
            color=color, cmap=CONTACT_MAP_CMAP, cutoff=cutoff, font_size=font_size, tick_step=tick_step,
            min_marker_pt=min_marker_pt, map_size_in=size_in * 0.5,
            title=f"Contacts {name_a}-{name_b}", n_ticks=n_ticks, show_ticks=show_ticks, tick_size=tick_size,
        )
        distances = self._pair_rows(contacts, label_a, label_b)["min_distance_nm"].to_numpy()
        ax_dist.hist(distances, bins=CONTACT_HIST_BINS, range=(0, cutoff), color=color)
        ax_dist.set_xlabel("min distance (nm)", fontsize=font_size)
        ax_dist.set_ylabel("contacts", fontsize=font_size)
        ax_dist.set_title("Contact distances", fontsize=font_size * CONTACT_TITLE_RATIO)

        profile, separation = stats["profile"], stats["separation"]
        for index, key in enumerate(keys):
            label, name = self.chains[key]["chain_id"], self._chain_name(key)
            line_color = CONTACT_PANEL_LINE_COLORS[index % len(CONTACT_PANEL_LINE_COLORS)]
            rows = profile[(profile["chain_id"] == label) & profile["sequence_index"].notna()]
            ax_profile.step(rows["resnum"], rows["n_contacts"], where="mid", color=line_color, label=f"chain {name}")
            values = separation.loc[separation["chain_id"] == label, "separation"]
            ax_sep.hist(values, bins=CONTACT_HIST_BINS, histtype="step", color=line_color, label=f"chain {name}")
        ax_profile.set_xlabel("residue number", fontsize=font_size)
        ax_profile.set_ylabel("contacts per residue", fontsize=font_size)
        ax_profile.set_title("Per-residue contacts", fontsize=font_size * CONTACT_TITLE_RATIO)
        ax_sep.set_xlabel("sequence separation |i - j|", fontsize=font_size)
        ax_sep.set_ylabel("intra-chain contacts", fontsize=font_size)
        ax_sep.set_title("Sequence separation", fontsize=font_size * CONTACT_TITLE_RATIO)
        if len(keys) > 1:
            ax_profile.legend(fontsize=font_size, frameon=False)
            ax_sep.legend(fontsize=font_size, frameon=False)
        for ax in (ax_dist, ax_profile, ax_sep):
            ax.tick_params(labelsize=font_size * CONTACT_TICK_RATIO)
        fig.suptitle(f"Contact summary {name_a}-{name_b} (cutoff {cutoff} nm)", fontsize=font_size * CONTACT_SUPTITLE_RATIO)
        fig.tight_layout()
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        return fig, axes
