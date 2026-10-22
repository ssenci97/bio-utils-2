"""Contact map plots: pair maps with grayscale background contacts, optionally cropped to the interface, and multi-chain block maps with chain outlines and optional pLDDT tracks."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

import itertools
from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import ScaledTranslation

from ._types import Moieties, PathLike
from .base import StructureState
from .constants import COMPLETE_MIN_SEPARATION, CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES
from .plot_style import (
    AXIS_LW, CONTACT_CROP_PAD, CONTACT_FONT_MAX_PT, CONTACT_FONT_MIN_PT, CONTACT_FONT_PT_PER_IN, CONTACT_GRAY_CMAP,
    CONTACT_GRAY_FACTOR, CONTACT_MAP_BACKGROUND, CONTACT_MAP_CMAP, CONTACT_MAP_COLOR, CONTACT_MAP_DPI,
    CONTACT_MAP_MIN_MARKER_PT, CONTACT_MAP_N_TICKS, CONTACT_MAP_SIZE_IN, CONTACT_MAP_STYLES, CONTACT_MAP_TICK_STEP,
    CONTACT_MAP_VALUES, CONTACT_MAX_ASPECT, CONTACT_MULTI_SIZE_IN, CONTACT_TICK_RATIO, CONTACT_TITLE_RATIO,
    COLORBAR_PAD, COLORBAR_WIDTH, LEGEND_GAP_RATIO, MAP_FILL_FRACTION, MULTI_BOX_COLOR, MULTI_BOX_LW, MULTI_INTRA_COLOR,
    MULTI_PAIR_COLORS, MULTI_TICK_MIN_GAP, MULTI_TICK_SPACING_EM, MULTI_TRACK_OFFSET, PLDDT_BAND_EDGES, PLDDT_BANDS, PLDDT_MISSING_COLOR,
    POINTS_PER_INCH, TICK_CHAR_EM, TICK_LENGTH_PT, TICK_STEPS,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.cm import ScalarMappable
    from matplotlib.colorbar import Colorbar
    from matplotlib.figure import Figure, SubFigure
    from matplotlib.transforms import Transform
    from matplotlib.typing import ColorType

Aspect = Literal["auto", "equal"] | float

###################################### LAYOUT HELPERS


def _font_size(size_in: float, font_size: float | None = None) -> float:
    """Return the base font size, proportional to the plot size unless given explicitly.

    Args:
        size_in: Plot size in inches.
        font_size: Explicit font size in points, or None.

    Returns:
        The font size in points.
    """
    if font_size is not None:
        return float(font_size)
    return float(np.clip(CONTACT_FONT_PT_PER_IN * size_in, CONTACT_FONT_MIN_PT, CONTACT_FONT_MAX_PT))


def _below_axis_transform(ax: Axes, gap_pt: float) -> Transform:
    """Return an axes-coordinate transform shifted down past the x tick labels and title.

    The shift is measured on the drawn x axis in points, so it holds when a layout engine resizes the axes.

    Args:
        ax: Matplotlib axes.
        gap_pt: Extra gap below the x axis in points.

    Returns:
        The shifted transform.
    """
    fig = ax.figure
    renderer = fig.canvas.get_renderer()  # type: ignore[attr-defined]  # defined by every drawing canvas (Agg, PDF)
    ax.apply_aspect()
    xaxis_box = ax.xaxis.get_tightbbox(renderer)
    below_px = ax.get_window_extent(renderer).y0 - xaxis_box.y0 if xaxis_box is not None else 0.0
    shift_in = (below_px * POINTS_PER_INCH / fig.dpi + gap_pt) / POINTS_PER_INCH
    return ax.transAxes + ScaledTranslation(0, -shift_in, fig.dpi_scale_trans)


def _attached_colorbar(ax: Axes, mappable: ScalarMappable, label: str, font_size: float) -> Colorbar:
    """Draw a colorbar flush against the right edge of ax, as tall as its fixed-aspect box.

    Args:
        ax: Matplotlib axes.
        mappable: Artist whose colors the bar shows.
        label: Colorbar label.
        font_size: Label font size in points.

    Returns:
        The colorbar.
    """
    cax = ax.inset_axes((1 + COLORBAR_PAD, 0, COLORBAR_WIDTH, 1))
    colorbar = ax.figure.colorbar(mappable, cax=cax)
    colorbar.set_label(label, fontsize=font_size)
    colorbar.ax.tick_params(labelsize=font_size * CONTACT_TICK_RATIO, length=TICK_LENGTH_PT, width=AXIS_LW)
    colorbar.outline.set_linewidth(AXIS_LW)
    return colorbar


def _marker_points(map_size_in: float, n_cols: int, n_rows: int, min_marker_pt: float) -> float:
    """Return the marker side so that one residue fills its cell, but at least min_marker_pt wide.

    Args:
        map_size_in: Map size in inches.
        n_cols: Number of columns.
        n_rows: Number of rows.
        min_marker_pt: Smallest marker side in points.

    Returns:
        The marker side in points.
    """
    return max(min_marker_pt, map_size_in * POINTS_PER_INCH * MAP_FILL_FRACTION / max(n_cols, n_rows, 1))


def _gray_limit(cutoff: float, gray_cutoff: float | None) -> float:
    """Return the upper distance of the grayscale background.

    Args:
        cutoff: Contact cutoff in nm.
        gray_cutoff: Explicit background limit in nm, or None for CONTACT_GRAY_FACTOR times cutoff.

    Returns:
        The limit in nm, never below cutoff.
    """
    if gray_cutoff is None:
        return cutoff * CONTACT_GRAY_FACTOR
    return max(gray_cutoff, cutoff)


def _tick_positions(numbers: Sequence[int], tick_step: int, n_ticks: int | None) -> np.ndarray:
    """Return the axis positions whose residue numbers are round: about n_ticks, or multiples of tick_step.

    Args:
        numbers: Residue number at each position along the axis.
        tick_step: Residue-number spacing, used when n_ticks is None.
        n_ticks: Approximate number of ticks, or None.

    Returns:
        The tick positions.

    Raises:
        ValueError: If n_ticks or tick_step is below 1.
    """
    array = np.asarray(numbers, dtype=int)
    if array.size == 0:
        return np.array([], dtype=int)
    low, high = array.min(), array.max()
    if n_ticks is not None:
        if n_ticks < 1:
            raise ValueError("n_ticks must be at least 1.")
        values = np.asarray(MaxNLocator(nbins=int(n_ticks), steps=TICK_STEPS, integer=True).tick_values(low, high))
    else:
        if tick_step < 1:
            raise ValueError("tick_step must be at least 1.")
        values = np.arange(-(-low // tick_step) * tick_step, high + 1, tick_step)
    first: dict[int, int] = {}
    for position, number in enumerate(array.tolist()):
        first.setdefault(number, position)
    positions = sorted({first[int(v)] for v in values if int(v) == v and int(v) in first})
    return np.array(positions or [0], dtype=int)


def _crop_window(mask: np.ndarray, n: int, pad: int) -> tuple[int, int]:
    """Return the window around the True entries of a 1D mask, padded and clamped to 0..n.

    Args:
        mask: Boolean mask along the axis.
        n: Axis length.
        pad: Positions added on each side.

    Returns:
        (start, stop), or the full axis when the mask is empty.
    """
    hits = np.flatnonzero(mask)
    if hits.size == 0:
        return 0, n
    return max(int(hits[0]) - pad, 0), min(int(hits[-1]) + pad + 1, n)


def _resolve_aspect(aspect: Aspect | None, n_rows: int, n_cols: int) -> Aspect:
    """Return the axes aspect: 'equal' unless the map is more elongated than CONTACT_MAX_ASPECT.

    Args:
        aspect: Explicit aspect, or None to decide from the shape.
        n_rows: Number of visible rows.
        n_cols: Number of visible columns.

    Returns:
        The explicit aspect, 'equal' or 'auto'.
    """
    if aspect is not None:
        return aspect
    ratio = max(n_rows, n_cols) / max(min(n_rows, n_cols), 1)
    return "auto" if ratio > CONTACT_MAX_ASPECT else "equal"


def plddt_colors(values: Sequence[float] | np.ndarray) -> list[str]:
    """Return AlphaFold pLDDT band colors for per-residue values in 0 to 100.

    Args:
        values: pLDDT values, NaN where missing.

    Returns:
        One color per value, gray where a value is NaN.
    """
    array = np.asarray(values, dtype=float)
    high, mid, low = PLDDT_BAND_EDGES
    colors = np.where(array >= high, PLDDT_BANDS[0][1],
                      np.where(array >= mid, PLDDT_BANDS[1][1],
                               np.where(array >= low, PLDDT_BANDS[2][1], PLDDT_BANDS[3][1])))
    return [PLDDT_MISSING_COLOR if np.isnan(v) else str(c) for v, c in zip(array, colors)]

###################################### DRAWING HELPERS


def _multi_masks(matrix: np.ndarray, sizes: Sequence[int], pairs: Sequence[tuple[int, int]], cutoff: float,
                 min_separation: int) -> tuple[np.ndarray, list[np.ndarray]]:
    """Return boolean masks over the block matrix of a multi-chain map.

    Args:
        matrix: Block distance matrix in nm.
        sizes: Residues per chain, in axis order.
        pairs: Chain index pairs, in color order.
        cutoff: Contact cutoff in nm.
        min_separation: Smallest sequence separation kept for intra-chain contacts.

    Returns:
        (gray, colored): gray marks every computed residue pair for the background. colored lists the intra-chain
        mask first, then one inter-chain mask per entry of pairs, each limited to contacts within cutoff.
    """
    n_total = sum(sizes)
    block_of = np.repeat(np.arange(len(sizes)), sizes)
    pair_code = np.full((len(sizes), len(sizes)), -1)
    for pair_idx, (i, j) in enumerate(pairs):
        pair_code[i, j] = pair_code[j, i] = pair_idx
    index = np.arange(n_total)
    same_block = block_of[:, None] == block_of[None, :]
    pair_index = pair_code[block_of[:, None], block_of[None, :]]
    far_enough = np.abs(index[:, None] - index[None, :]) >= min_separation
    close = matrix <= cutoff
    gray = ~np.isnan(matrix) & (~same_block | far_enough)
    intra_close = close & same_block & far_enough
    inter_close = close & ~same_block
    colored = [intra_close] + [inter_close & (pair_index == pair_idx) for pair_idx in range(len(pairs))]
    return gray, colored


def _multi_legend(ax: Axes, names: Sequence[str], pairs: Sequence[tuple[int, int]],
                  pair_colors: Sequence[ColorType], intra_color: ColorType, plddt: bool, font_size: float) -> None:
    """Draw the legend of intra-chain, chain-pair and optional pLDDT colors on a multi-chain map.

    Args:
        ax: Matplotlib axes.
        names: Chain labels in axis order.
        pairs: Chain index pairs, in color order.
        pair_colors: One color per chain pair.
        intra_color: Color of intra-chain contacts.
        plddt: Whether to add the pLDDT bands.
        font_size: Legend font size in points.
    """
    handles = [Patch(facecolor=intra_color, label="intra-chain")]
    handles += [Patch(facecolor=pair_colors[pair_idx], label=f"{names[i]}-{names[j]}")
                for pair_idx, (i, j) in enumerate(pairs)]
    if plddt:
        handles += [Patch(facecolor=color, label=f"pLDDT {label}") for label, color in PLDDT_BANDS]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0),
              bbox_transform=_below_axis_transform(ax, font_size * LEGEND_GAP_RATIO), ncol=min(4, len(handles)),
              fontsize=font_size, frameon=False, borderaxespad=0, borderpad=0,
              handlelength=1.0, handletextpad=0.4, columnspacing=1.0)


def _draw_multi_boundaries(ax: Axes, offsets: np.ndarray, sizes: Sequence[int], n_total: int, box_color: ColorType,
                           box_lw: float) -> None:
    """Draw dashed chain boundaries and a thin square around each chain of a multi-chain map.

    Args:
        ax: Matplotlib axes.
        offsets: Start index of each chain.
        sizes: Residues per chain.
        n_total: Residues along each axis.
        box_color: Line color.
        box_lw: Line width in points.
    """
    boundaries = offsets[1:] - 0.5
    dashes = (0, (3, 2))
    ax.hlines(boundaries, -0.5, n_total - 0.5, colors=box_color, linestyles=dashes, linewidth=box_lw, zorder=3)
    ax.vlines(boundaries, -0.5, n_total - 0.5, colors=box_color, linestyles=dashes, linewidth=box_lw, zorder=3)
    for start, length in zip(offsets, sizes):
        ax.add_patch(Rectangle((start - 0.5, start - 0.5), length, length, fill=False,
                               edgecolor=box_color, linewidth=box_lw, zorder=3))


class ContactPlotMixin(StructureState):
    """Contact maps and multi-chain block maps."""

    def _residue_bfactor_as_plddt(self, chain_key: str) -> np.ndarray:
        """Return per-residue pLDDT of a chain from the C-alpha B-factor, or the atom mean without C-alpha.

        AlphaFold models store pLDDT in the B-factor column; for experimental structures these are B-factors.

        Args:
            chain_key: Chain key.

        Returns:
            One value per polymer residue, NaN where missing.
        """
        values = np.full(self.chains[chain_key].sequence_length, np.nan)
        for position, residue in enumerate(self.chains[chain_key].residues):
            if "CA" in residue:
                values[position] = residue["CA"].get_bfactor()
            else:
                atoms = list(residue.get_atoms())
                if atoms:
                    values[position] = np.mean([atom.get_bfactor() for atom in atoms])
        return values

    def _quoted_auth_id(self, chain_key: str) -> str:
        """Return the quoted author chain ID used in plot text.

        Args:
            chain_key: Chain key.

        Returns:
            The ID in quotes, such as 'A'.
        """
        return f"'{self.chains[chain_key].auth_chain_id}'"

    def _pair_title(self, chain_a: str, chain_b: str, cutoff: float) -> str:
        """Return the default title of a pair map.

        Args:
            chain_a: Chain key on the y-axis.
            chain_b: Chain key on the x-axis.
            cutoff: Contact cutoff in nm.

        Returns:
            "Chain 'A' contacts" or "'A'–'B' contacts", with the cutoff.
        """
        names = (f"Chain {self._quoted_auth_id(chain_a)}" if chain_a == chain_b
                 else f"{self._quoted_auth_id(chain_a)}–{self._quoted_auth_id(chain_b)}")
        return f"{names} contacts (≤ {cutoff:g} nm)"

    def _tick_labels(self, chain_key: str, start: int, stop: int, tick_step: int,
                     n_ticks: int | None = None) -> tuple[np.ndarray, list[str]]:
        """Return tick positions and residue-number labels for positions start..stop-1 of one axis.

        Ticks fall on round author residue numbers, so they read naturally whatever the numbering offset.

        Args:
            chain_key: Chain key whose residue numbers label the ticks.
            start: First visible sequence position.
            stop: End of the visible window, exclusive.
            tick_step: Residue-number spacing of ticks, used when n_ticks is None.
            n_ticks: Approximate number of ticks along the window, or None.

        Returns:
            (tick positions, tick labels).
        """
        residues = self.chains[chain_key].residues
        numbers = [residues[p].id[1] if p < len(residues) else p + 1 for p in range(start, stop)]
        positions = _tick_positions(numbers, tick_step, n_ticks)
        return positions + start, [str(numbers[p]) for p in positions]

    @staticmethod
    def _style_axes(ax: Axes, x_ticks: Sequence[int] | np.ndarray, x_tick_labels: Sequence[str],
                    y_ticks: Sequence[int] | np.ndarray, y_tick_labels: Sequence[str], *, show_ticks: bool,
                    tick_size: float | None, tick_rotation: float, xlabel: str, ylabel: str,
                    font_size: float) -> None:
        """Apply tick positions, tick labels and axis titles to an axes.

        Args:
            ax: Matplotlib axes.
            x_ticks: X tick positions.
            x_tick_labels: X tick labels.
            y_ticks: Y tick positions.
            y_tick_labels: Y tick labels.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title.
            ylabel: Y-axis title.
            font_size: Base font size in points.
        """
        if show_ticks:
            ax.set_xticks(x_ticks, x_tick_labels)
            ax.set_yticks(y_ticks, y_tick_labels)
            ax.tick_params(labelsize=tick_size if tick_size is not None else font_size * CONTACT_TICK_RATIO,
                           direction="out", length=TICK_LENGTH_PT, width=AXIS_LW)
            ax.tick_params(axis="x", labelrotation=tick_rotation)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(AXIS_LW)
        ax.set_xlabel(xlabel, fontsize=font_size)
        ax.set_ylabel(ylabel, fontsize=font_size)

    @staticmethod
    def _paint(ax: Axes, matrix: np.ndarray, gray_matrix: np.ndarray | None, value: str, style: str,
               background: ColorType, color: ColorType, cmap: str, cutoff: float, gray_cutoff: float | None,
               marker_pt: float) -> ScalarMappable | None:
        """Draw a pair matrix on ax over an optional grayscale background.

        Distance values use a fixed scale from 0 to cutoff, so plots with the same cutoff are comparable.

        Args:
            ax: Matplotlib axes.
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

        Returns:
            The colored artist, or None for scatter contacts drawn in one color.
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
            # transparent gaps let the grayscale background show through
            colormap.set_bad(background if gray_matrix is None else (0, 0, 0, 0))
            return ax.imshow(np.ma.masked_invalid(matrix), cmap=colormap, vmin=0, vmax=cutoff, interpolation="nearest")
        marker_area = marker_pt ** 2
        if gray_matrix is not None:
            rows, cols = np.nonzero(~np.isnan(gray_matrix))
            ax.scatter(cols, rows, c=gray_matrix[rows, cols], cmap=CONTACT_GRAY_CMAP, vmin=0, vmax=gray_cutoff,
                       s=marker_area, marker="s", linewidths=0, zorder=1)
        rows, cols = np.nonzero(~np.isnan(matrix))
        if value == "contact":
            ax.scatter(cols, rows, s=marker_area, marker="s", color=color, linewidths=0, zorder=2)
            return None
        return ax.scatter(cols, rows, c=matrix[rows, cols], cmap=cmap, vmin=0, vmax=cutoff,
                          s=marker_area, marker="s", linewidths=0, zorder=2)

    def _draw_pair(self, ax: Axes, matrix: np.ndarray, chain_a: str, chain_b: str, *, value: str, style: str,
                   background: ColorType, color: ColorType, cmap: str, cutoff: float, font_size: float,
                   tick_step: int, min_marker_pt: float, map_size_in: float, title: str,
                   gray_matrix: np.ndarray | None = None, gray_cutoff: float | None = None,
                   n_ticks: int | None = None, show_ticks: bool = True, tick_size: float | None = None,
                   tick_rotation: float = 0, xlabel: str | None = None, ylabel: str | None = None,
                   crop: bool = False, aspect: Aspect | None = None) -> None:
        """Draw a pair contact map on ax with quoted chain-ID axis titles, optionally cropped to the contacts.

        Args:
            ax: Matplotlib axes.
            matrix: Contact matrix, NaN for non-contacts.
            chain_a: Chain key on the y-axis.
            chain_b: Chain key on the x-axis.
            value: 'contact' or 'distance'.
            style: 'scatter' or 'heatmap'.
            background: Background color.
            color: Contact color, for value 'contact'.
            cmap: Colormap name, for value 'distance'.
            cutoff: Distance cutoff in nm.
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
            crop: Whether to limit both axes to the residues in contact, padded by CONTACT_CROP_PAD.
            aspect: 'equal', 'auto', or None to decide from the window shape.
        """
        n_rows, n_cols = matrix.shape
        row_start, row_stop, col_start, col_stop = 0, n_rows, 0, n_cols
        if crop:
            close = ~np.isnan(matrix)
            row_start, row_stop = _crop_window(close.any(axis=1), n_rows, CONTACT_CROP_PAD)
            col_start, col_stop = _crop_window(close.any(axis=0), n_cols, CONTACT_CROP_PAD)
        n_visible_rows, n_visible_cols = row_stop - row_start, col_stop - col_start
        ax.set_facecolor(background)
        marker_pt = _marker_points(map_size_in, n_visible_cols, n_visible_rows, min_marker_pt)
        mappable = self._paint(ax, matrix, gray_matrix, value, style, background, color, cmap,
                               cutoff, gray_cutoff, marker_pt)
        # imshow resets limits and aspect, so both are applied after painting
        ax.set_xlim(col_start - 0.5, col_stop - 0.5)
        ax.set_ylim(row_stop - 0.5, row_start - 0.5)
        ax.set_aspect(_resolve_aspect(aspect, n_visible_rows, n_visible_cols), adjustable="box")
        if mappable is not None and value == "distance":
            _attached_colorbar(ax, mappable, "min distance (nm)", font_size)
        x_ticks, x_tick_labels = self._tick_labels(chain_b, col_start, col_stop, tick_step, n_ticks)
        y_ticks, y_tick_labels = self._tick_labels(chain_a, row_start, row_stop, tick_step, n_ticks)
        self._style_axes(ax, x_ticks, x_tick_labels, y_ticks, y_tick_labels, show_ticks=show_ticks, tick_size=tick_size,
                         tick_rotation=tick_rotation, xlabel=xlabel or self._quoted_auth_id(chain_b),
                         ylabel=ylabel or self._quoted_auth_id(chain_a), font_size=font_size)
        ax.set_title(title, fontsize=font_size * CONTACT_TITLE_RATIO)

    def _checked_matrix(self, matrix: np.ndarray, chain_a: str, chain_b: str) -> np.ndarray:
        """Return a precomputed matrix as float after checking its shape matches the two chains.

        Args:
            matrix: Precomputed pair matrix.
            chain_a: Chain key of the rows.
            chain_b: Chain key of the columns.

        Returns:
            The matrix as float.

        Raises:
            ValueError: If the shape does not match the chains' sequence lengths.
        """
        matrix = np.asarray(matrix, dtype=float)
        expected = (self.chains[chain_a].sequence_length, self.chains[chain_b].sequence_length)
        if matrix.shape != expected:
            raise ValueError(f"matrix shape {matrix.shape} does not match chains {expected}.")
        return matrix

###################################### PUBLIC PLOTS

    def plot_contact_map(self, chain_a: str = "A", chain_b: str | None = None, cutoff: float = CONTACTS_CUTOFF_NM,
                         moieties: Moieties = DEFAULT_CONTACT_MOIETIES, scheme: str = "ca", value: str = "contact",
                         style: str = "scatter", gray: bool = True, gray_cutoff: float | None = None,
                         background: ColorType = CONTACT_MAP_BACKGROUND, color: ColorType = CONTACT_MAP_COLOR,
                         cmap: str = CONTACT_MAP_CMAP, size_in: float = CONTACT_MAP_SIZE_IN,
                         font_size: float | None = None, tick_step: int = CONTACT_MAP_TICK_STEP,
                         n_ticks: int | None = CONTACT_MAP_N_TICKS, show_ticks: bool = True,
                         tick_size: float | None = None, tick_rotation: float = 0, xlabel: str | None = None,
                         ylabel: str | None = None, min_marker_pt: float = CONTACT_MAP_MIN_MARKER_PT,
                         crop: bool = False, aspect: Aspect | None = None, dpi: float = CONTACT_MAP_DPI,
                         title: str | None = None, ax: Axes | None = None, save_path: PathLike | None = None,
                         matrix: np.ndarray | None = None) -> tuple[Figure | SubFigure, Axes]:
        """Plot a contact map of one chain or a chain pair, each axis titled by its chain ID.

        With gray True, residue pairs up to gray_cutoff are drawn in gray, darker when closer, under the contacts
        within cutoff, which are drawn in color or colored by distance. Ticks fall on round author residue numbers.
        crop=True zooms onto the residues in contact, which makes sparse interfaces readable.

        Args:
            chain_a: Chain on the y-axis, standard or author ID.
            chain_b: Chain on the x-axis, or None for an intra-chain map.
            cutoff: Distance cutoff in nm for the colored contacts.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            value: 'contact' or 'distance' (minimum distances in nm).
            style: 'scatter' or 'heatmap'.
            gray: Whether to draw the grayscale background of all contacts up to gray_cutoff.
            gray_cutoff: Upper distance of the background in nm, or None for CONTACT_GRAY_FACTOR times cutoff.
            background: Background color.
            color: Contact color when value is 'contact'.
            cmap: Colormap name when value is 'distance'.
            size_in: Figure width and height in inches.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Residue-number spacing of ticks, used when n_ticks is None.
            n_ticks: Approximate number of ticks per axis; overrides tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for the quoted chain ID of chain_b.
            ylabel: Y-axis title, or None for the quoted chain ID of chain_a.
            min_marker_pt: Smallest marker side in points for the scatter style.
            crop: Whether to zoom onto the residues in contact, padded by CONTACT_CROP_PAD.
            aspect: 'equal', 'auto', or None for 'equal' unless the map is very elongated (e.g. protein-DNA).
            dpi: Figure resolution.
            title: Title, or None for a default.
            ax: Axes to draw on, or None to create a figure.
            save_path: File to save the figure to, or None.
            matrix: Precomputed pair matrix from contact_matrix, NaN for non-contacts.

        Returns:
            (figure, axes).

        Raises:
            ValueError: If value or style is unknown.
        """
        if value not in CONTACT_MAP_VALUES:
            raise ValueError(f"value must be one of {CONTACT_MAP_VALUES}.")
        if style not in CONTACT_MAP_STYLES:
            raise ValueError(f"style must be one of {CONTACT_MAP_STYLES}.")
        font_size = _font_size(size_in, font_size)
        chain_a, chain_b, _, _ = self._resolve_pair(chain_a, chain_b)
        gray_matrix = None
        if gray:
            gray_cutoff = _gray_limit(cutoff, gray_cutoff)
            gray_matrix = self.contact_matrix(chain_a, chain_b, gray_cutoff, moieties, scheme)["matrix"]
        if matrix is not None:
            matrix = self._checked_matrix(matrix, chain_a, chain_b)
        elif gray_matrix is not None:
            matrix = np.where(gray_matrix <= cutoff, gray_matrix, np.nan)
        else:
            matrix = self.contact_matrix(chain_a, chain_b, cutoff, moieties, scheme)["matrix"]
        fig: Figure | SubFigure
        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi, layout="constrained")
        else:
            fig = ax.figure
        self._draw_pair(
            ax, matrix, chain_a, chain_b, value=value, style=style, background=background, color=color,
            cmap=cmap, cutoff=cutoff, font_size=font_size, tick_step=tick_step,
            min_marker_pt=min_marker_pt, map_size_in=size_in,
            title=title or self._pair_title(chain_a, chain_b, cutoff),
            gray_matrix=gray_matrix, gray_cutoff=gray_cutoff, n_ticks=n_ticks, show_ticks=show_ticks,
            tick_size=tick_size, tick_rotation=tick_rotation, xlabel=xlabel, ylabel=ylabel, crop=crop, aspect=aspect,
        )
        if save_path is not None:
            fig.figure.savefig(save_path, dpi=dpi)
        return fig, ax

    def plot_multi_contact_map(self, chain_ids: Sequence[str] = ("A", "B"), cutoff: float = CONTACTS_CUTOFF_NM,
                               moieties: Moieties = DEFAULT_CONTACT_MOIETIES, scheme: str = "ca",
                               gray_cutoff: float | None = None, pair_colors: Sequence[ColorType] = MULTI_PAIR_COLORS,
                               intra_color: ColorType = MULTI_INTRA_COLOR,
                               background: ColorType = CONTACT_MAP_BACKGROUND, box_color: ColorType = MULTI_BOX_COLOR,
                               box_lw: float = MULTI_BOX_LW, plddt: bool = False,
                               min_separation: int = COMPLETE_MIN_SEPARATION, size_in: float | None = None,
                               font_size: float | None = None, tick_step: int = CONTACT_MAP_TICK_STEP,
                               n_ticks: int | None = CONTACT_MAP_N_TICKS, show_ticks: bool = True,
                               tick_size: float | None = None, tick_rotation: float = 0, xlabel: str | None = None,
                               ylabel: str | None = None, min_marker_pt: float = CONTACT_MAP_MIN_MARKER_PT,
                               show_color_legend: bool = False, show_distance_legend: bool = False,
                               dpi: float = CONTACT_MAP_DPI, title: str | None = None, ax: Axes | None = None,
                               save_path: PathLike | None = None) -> tuple[Figure | SubFigure, Axes]:
        """Plot several chains on one block contact map, with dashed chain boundaries and a square around each chain.

        Residue pairs up to gray_cutoff are drawn in gray by distance. Within cutoff, intra-chain contacts use
        intra_color and inter-chain contacts one color per chain pair, in itertools.combinations order. Chain IDs
        label the blocks along the top; plddt=True adds pLDDT bands along the top and left edges.

        Args:
            chain_ids: Two or more different chains, standard or author IDs, in axis order.
            cutoff: Distance cutoff in nm for colored contacts.
            moieties: Biomolecule types to include.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            gray_cutoff: Upper distance of the background in nm, or None for CONTACT_GRAY_FACTOR times cutoff.
            pair_colors: Inter-chain contact colors, one per chain pair.
            intra_color: Color of intra-chain contacts within cutoff.
            background: Background color.
            box_color: Color of the outlines and dashed boundaries.
            box_lw: Outline and dashed boundary width in points.
            plddt: Whether to draw pLDDT tracks, read from the B-factor column.
            min_separation: Smallest sequence separation kept for intra-chain contacts.
            size_in: Figure width and height in inches, or None for a size by number of chains.
            font_size: Base font size in points, or None to scale with size_in.
            tick_step: Residue-number spacing of ticks, used when n_ticks is None.
            n_ticks: Approximate number of ticks per chain; overrides tick_step.
            show_ticks: Whether to draw tick marks and labels.
            tick_size: Tick label font size, or None for a fixed share of font_size.
            tick_rotation: X tick label rotation in degrees.
            xlabel: X-axis title, or None for 'residue'.
            ylabel: Y-axis title, or None for 'residue'.
            min_marker_pt: Smallest marker side in points.
            show_color_legend: Whether to draw the legend of intra-chain, chain-pair and pLDDT colors.
            show_distance_legend: Whether to draw the grayscale distance colorbar.
            dpi: Figure resolution.
            title: Title, or None for a default.
            ax: Axes to draw on, or None to create a figure.
            save_path: File to save the figure to, or None.

        Returns:
            (figure, axes).

        Raises:
            ValueError: If fewer than two different chains or too few pair colors are given.
        """
        keys = [self._require_chain(c) for c in chain_ids]
        if len(keys) < 2 or len(set(keys)) != len(keys):
            raise ValueError("plot_multi_contact_map needs two or more different chains.")
        pairs = list(itertools.combinations(range(len(keys)), 2))
        if len(pair_colors) < len(pairs):
            raise ValueError(f"pair_colors needs {len(pairs)} colors for {len(keys)} chains.")
        size_in = size_in or CONTACT_MULTI_SIZE_IN.get(len(keys), max(CONTACT_MULTI_SIZE_IN.values()))
        font_size = _font_size(size_in, font_size)
        names = [self._quoted_auth_id(k) for k in keys]
        gray_cutoff = _gray_limit(cutoff, gray_cutoff)
        block = self._block_matrix(keys, gray_cutoff, scheme, moieties)
        matrix, sizes, offsets = block["matrix"], block["sizes"], block["offsets"]
        n_total = sum(sizes)
        gray_mask, colored_masks = _multi_masks(matrix, sizes, pairs, cutoff, min_separation)
        colors = [intra_color] + list(pair_colors[:len(pairs)])

        fig: Figure | SubFigure
        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi, layout="constrained")
        else:
            fig = ax.figure
        margin = MULTI_TRACK_OFFSET + 0.5 if plddt else 0.5
        ax.set_facecolor(background)
        ax.set_xlim(-margin, n_total - 0.5)
        ax.set_ylim(n_total - 0.5, -margin)
        ax.set_aspect("equal", adjustable="box")

        marker_area = _marker_points(size_in, n_total, n_total, min_marker_pt) ** 2
        rows, cols = np.nonzero(gray_mask)
        background_scatter = ax.scatter(cols, rows, c=matrix[rows, cols], cmap=CONTACT_GRAY_CMAP,
                                        vmin=0, vmax=gray_cutoff, s=marker_area, marker="s", linewidths=0, zorder=1)
        for mask, color in zip(colored_masks, colors):
            rows, cols = np.nonzero(mask)
            ax.scatter(cols, rows, s=marker_area, marker="s", color=color, linewidths=0, zorder=2)
        if plddt:
            self._draw_plddt_tracks(ax, keys, n_total, marker_area)
        _draw_multi_boundaries(ax, offsets, sizes, n_total, box_color, box_lw)

        tick_pt = tick_size if tick_size is not None else font_size * CONTACT_TICK_RATIO
        positions_per_pt = n_total / (size_in * MAP_FILL_FRACTION * POINTS_PER_INCH)
        ticks, tick_labels = self._multi_ticks(keys, offsets, sizes, tick_step, n_ticks, tick_pt * positions_per_pt)
        self._style_axes(ax, ticks, tick_labels, ticks, tick_labels, show_ticks=show_ticks,
                         tick_size=tick_size, tick_rotation=tick_rotation,
                         xlabel=xlabel or "residue", ylabel=ylabel or "residue", font_size=font_size)
        # chains too short to fit their ID stay unlabeled rather than overlapping a neighbor
        shown = [i for i, n in enumerate(sizes) if n >= MULTI_TICK_MIN_GAP * n_total]
        chain_axis = ax.secondary_xaxis("top")
        chain_axis.set_xticks([offsets[i] + (sizes[i] - 1) / 2 for i in shown], [names[i] for i in shown])
        chain_axis.tick_params(length=0, labelsize=font_size)
        chain_axis.spines["top"].set_linewidth(AXIS_LW)
        ax.set_title(title or f"{'–'.join(names)} contacts (≤ {cutoff:g} nm)", fontsize=font_size * CONTACT_TITLE_RATIO)

        if show_color_legend:
            _multi_legend(ax, names, pairs, pair_colors, intra_color, plddt, font_size)
        if show_distance_legend:
            _attached_colorbar(ax, background_scatter, "distance (nm)", font_size)
        if save_path is not None:
            fig.figure.savefig(save_path, dpi=dpi)
        return fig, ax

    def _draw_plddt_tracks(self, ax: Axes, keys: Sequence[str], n_total: int, marker_area: float) -> None:
        """Draw pLDDT bands along the top and left edges of a multi-chain map, from each residue's B-factor.

        Args:
            ax: Matplotlib axes.
            keys: Chain keys in axis order.
            n_total: Number of residues along each axis.
            marker_area: Marker area in square points.
        """
        colors = plddt_colors(np.concatenate([self._residue_bfactor_as_plddt(k) for k in keys]))
        index = np.arange(n_total)
        track = np.full(n_total, -MULTI_TRACK_OFFSET)
        for x, y in ((index, track), (track, index)):
            ax.scatter(x, y, c=colors, s=marker_area, marker="s", linewidths=0, clip_on=False, zorder=2)

    def _multi_ticks(self, keys: Sequence[str], offsets: np.ndarray, sizes: Sequence[int], tick_step: int,
                     n_ticks: int | None, em_positions: float) -> tuple[np.ndarray, list[str]]:
        """Return the tick positions and residue-number labels of a multi-chain map, shared by both axes.

        A tick whose label would come closer than MULTI_TICK_SPACING_EM to the previous label is dropped, so labels
        never collide at chain boundaries or inside short chains, whatever the figure size.

        Args:
            keys: Chain keys in axis order.
            offsets: Start index of each chain in the block matrix.
            sizes: Number of residues of each chain.
            tick_step: Residue-number spacing of ticks, used when n_ticks is None.
            n_ticks: Approximate number of ticks per chain, or None.
            em_positions: Width of one tick-label em in matrix positions.

        Returns:
            (tick positions, tick labels).
        """
        def half_width(label: str) -> float:
            """Return half the drawn width of a tick label.

            Args:
                label: Tick label.

            Returns:
                The half width in matrix positions.
            """
            return len(label) * TICK_CHAR_EM * em_positions / 2

        ticks: list[int] = []
        tick_labels: list[str] = []
        for key, offset, n in zip(keys, offsets, sizes):
            positions, chain_labels = self._tick_labels(key, 0, n, tick_step, n_ticks)
            for position, label in zip(positions + offset, chain_labels):
                min_gap = (half_width(tick_labels[-1]) + half_width(label) + MULTI_TICK_SPACING_EM * em_positions
                           if ticks else 0)
                if not ticks or position - ticks[-1] >= min_gap:
                    ticks.append(position)
                    tick_labels.append(label)
        return np.array(ticks, dtype=int), tick_labels
