"""Panel grids of plot functions with lettered panels and separator lines, exported as multi-page PDFs with optional rasterization."""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

import math
import string
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, cast

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase
from matplotlib.collections import Collection
from matplotlib.figure import Figure, SubFigure
from matplotlib.layout_engine import ConstrainedLayoutEngine
from matplotlib.image import AxesImage
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox, ScaledTranslation, Transform

from .plot_style import (
    GRID_GUTTER_MM, GRID_HSPACE, GRID_LABEL_SIZE, GRID_LINE_COLOR, GRID_LINE_WIDTH, GRID_TITLE_SIZE, GRID_WSPACE,
    MM_PER_INCH, PANEL_PDF_NAME, PANEL_SIZE_IN, PDF_RC, RASTER_DPI,
)
from ._types import PathLike
from .utils import log_info

PlotFunc = Callable[..., Any]

###################################### RASTERIZATION


def rasterize_figure(fig: Figure) -> Figure:
    """Rasterize collections and images of every axes in a figure, keeping text and lines vector.

    Args:
        fig: Matplotlib figure, changed in place.

    Returns:
        The same figure.
    """
    for ax in fig.axes:
        for artist in ax.get_children():
            if isinstance(artist, (Collection, AxesImage)):
                artist.set_rasterized(True)
    return fig

###################################### PANEL LAYOUT


def _renderer(fig: Figure | SubFigure) -> RendererBase:
    """Return the renderer of a figure drawn on a canvas that provides one, such as Agg.

    Args:
        fig: Matplotlib figure.

    Returns:
        The canvas renderer.
    """
    renderer: RendererBase = getattr(fig.canvas, "get_renderer")()
    return renderer


def _label_transform(ax: Axes) -> Transform:
    """Return a transform for anchors in axes coordinates, shifted left past the y tick labels and y title of ax.

    The shift does not depend on the axes width, so it holds when constrained layout later resizes the axes.

    Args:
        ax: Panel axes.

    Returns:
        The shifted axes transform.
    """
    fig = ax.figure
    renderer = _renderer(fig)
    box = ax.yaxis.get_tightbbox(renderer)
    left_px = max(ax.get_window_extent(renderer).x0 - box.x0, 0.0) if box is not None else 0.0
    return ax.transAxes + ScaledTranslation(-left_px / fig.dpi, 0, fig.dpi_scale_trans)


def _left_align_title(ax: Axes) -> None:
    """Move a centered axes title to the left edge, so it starts right of the panel letter and never runs under it.

    Args:
        ax: Panel axes.
    """
    title = ax.get_title()
    if title:
        fontsize = ax.title.get_fontsize()
        ax.set_title("")
        ax.set_title(title, loc="left", fontsize=fontsize)


def _draw_separators(fig: Figure, extents: Sequence[Bbox], n_rows: int, ncols: int, n_panels: int) -> None:
    """Draw one straight gray line through the middle of each gap between grid columns and rows.

    Args:
        fig: Matplotlib figure.
        extents: Full extent of each panel in figure coordinates.
        n_rows: Number of grid rows.
        ncols: Number of grid columns.
        n_panels: Number of drawn panels.
    """
    fig_w, fig_h = fig.get_size_inches()
    dx = GRID_GUTTER_MM / MM_PER_INCH / fig_w
    dy = GRID_GUTTER_MM / MM_PER_INCH / fig_h
    bottom, top = min(e.y0 for e in extents) + dy, max(e.y1 for e in extents) - dy
    left, right = min(e.x0 for e in extents) + dx, max(e.x1 for e in extents) - dx
    cells = {divmod(i, ncols): extent for i, extent in enumerate(extents[:n_panels])}
    for c in range(ncols - 1):
        lefts = [e.x1 for (r, col), e in cells.items() if col == c]
        rights = [e.x0 for (r, col), e in cells.items() if col == c + 1]
        if lefts and rights:
            x = (max(lefts) + min(rights)) / 2.0
            fig.add_artist(Line2D([x, x], [bottom, top], transform=fig.transFigure,
                                  color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))
    for r in range(n_rows - 1):
        uppers = [e.y0 for (row, c), e in cells.items() if row == r]
        lowers = [e.y1 for (row, c), e in cells.items() if row == r + 1]
        if uppers and lowers:
            y = (min(uppers) + max(lowers)) / 2.0
            fig.add_artist(Line2D([left, right], [y, y], transform=fig.transFigure,
                                  color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))


def build_panel_grid(plot_funcs: Sequence[PlotFunc], sizes: float | Sequence[float] = PANEL_SIZE_IN, ncols: int = 2,
                     panel_labels: bool = True, main_title: str | None = None, wspace: float = GRID_WSPACE,
                     hspace: float = GRID_HSPACE) -> Figure:
    """Build a figure with one panel per plot function, lettered A, B, ..., separated by thin gray lines.

    Panels are placed by constrained layout, so colorbars, legends, and labels never overlap a neighbor.
    Lettered panels get left-aligned titles next to their letter. Each separator is one straight line
    through the middle of the gap between neighboring columns or rows.

    Args:
        plot_funcs: Callables drawing on a given Axes, called as func(ax=ax).
        sizes: Panel size in inches, or one size per panel.
        ncols: Number of grid columns.
        panel_labels: Whether to letter the panels.
        main_title: Optional figure title.
        wspace: Horizontal gap between panels, as a fraction of the figure width.
        hspace: Vertical gap between panels, as a fraction of the figure height.

    Returns:
        The panel grid figure.
    """
    n_panels = len(plot_funcs)
    panel_sizes = [float(sizes)] * n_panels if isinstance(sizes, (int, float)) else [float(s) for s in sizes]
    n_rows = math.ceil(n_panels / ncols)
    full = panel_sizes + [panel_sizes[-1]] * (n_rows * ncols - n_panels)
    grid = [full[r * ncols:(r + 1) * ncols] for r in range(n_rows)]
    row_heights = [max(row) for row in grid]
    col_widths = [max(grid[r][c] for r in range(n_rows)) for c in range(ncols)]
    fig = plt.figure(figsize=(sum(col_widths), sum(row_heights)), layout="constrained")
    engine = fig.get_layout_engine()
    if isinstance(engine, ConstrainedLayoutEngine):
        engine.set(wspace=wspace, hspace=hspace)
    gs = fig.add_gridspec(n_rows, ncols, width_ratios=col_widths, height_ratios=row_heights)
    if main_title:
        fig.suptitle(main_title, fontsize=GRID_TITLE_SIZE, fontweight="bold")
    axes = []
    for i, func in enumerate(plot_funcs):
        row, col = divmod(i, ncols)
        ax = fig.add_subplot(gs[row, col])
        func(ax=ax)
        if panel_labels and i >= len(string.ascii_uppercase):
            log_info(f"Only {len(string.ascii_uppercase)} panels can be lettered; panel {i + 1} is unlabeled.")
        if panel_labels and i < len(string.ascii_uppercase):
            _left_align_title(ax)
            ax.text(0, 1, string.ascii_uppercase[i], transform=_label_transform(ax), fontsize=GRID_LABEL_SIZE,
                    fontweight="bold", va="bottom", ha="right")
        axes.append(ax)
    fig.canvas.draw()
    renderer = _renderer(fig)
    to_fig = fig.transFigure.inverted()
    extents = [to_fig.transform_bbox(ax.get_tightbbox(renderer) or ax.get_window_extent(renderer)) for ax in axes]
    _draw_separators(fig, extents, n_rows, ncols, len(axes))
    return fig


def export_panels_pdf(grids: Iterable[tuple[Sequence[PlotFunc], str | None]], outdir: PathLike,
                      filename: str = PANEL_PDF_NAME, raster: bool = True, dpi: int = RASTER_DPI,
                      **grid_kwargs: Any) -> Path:
    """Write one PDF page per panel grid, each at the grid's exact figure size.

    Args:
        grids: Iterable of (plot_funcs, main_title) pairs, one per page.
        outdir: Output directory, created if missing.
        filename: PDF file name.
        raster: Whether to rasterize collections and images.
        dpi: Resolution of rasterized content.
        **grid_kwargs: Passed to build_panel_grid.

    Returns:
        Path of the PDF.
    """
    path = Path(outdir).expanduser() / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    # PDF_RC holds valid rc keys; matplotlib types them as a Literal union that a plain dict can't match
    with mpl.rc_context(cast("dict[Any, Any]", PDF_RC)), PdfPages(path) as pdf:
        for plot_funcs, main_title in grids:
            fig = build_panel_grid(plot_funcs, main_title=main_title, **grid_kwargs)
            if raster:
                rasterize_figure(fig)
            pdf.savefig(fig, dpi=dpi)
            plt.close(fig)
    return path
