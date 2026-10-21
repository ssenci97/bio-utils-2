"""Panel grids of plot functions with lettered panels and separator lines, exported as multi-page PDFs with optional rasterization."""

###################################### IMPORTS AND CONSTANTS

import math
import string
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import Collection
from matplotlib.image import AxesImage
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

from .plot_style import (
    GRID_GUTTER_MM, GRID_HSPACE, GRID_LABEL_SIZE, GRID_LINE_COLOR, GRID_LINE_WIDTH, GRID_TITLE_SIZE, GRID_WSPACE,
    MM_PER_INCH, PANEL_PDF_NAME, PANEL_SIZE_IN, RASTER_DPI,
)
from .utils import log_info

###################################### RASTERIZATION


def rasterize_figure(fig):
    """Rasterize collections and images of every axes in a figure, keeping text and lines vector, and return it."""
    for ax in fig.axes:
        for artist in ax.get_children():
            if isinstance(artist, (Collection, AxesImage)):
                artist.set_rasterized(True)
    return fig

###################################### PANEL LAYOUT


def _label_transform(ax):
    """Return a transform for anchors in axes coordinates, shifted left past the y tick labels, y title, and title of ax."""
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    ax.apply_aspect()
    boxes = [ax.yaxis.get_tightbbox(renderer), ax.title.get_window_extent(renderer) if ax.get_title() else None]
    left_px = max([ax.get_window_extent(renderer).x0 - box.x0 for box in boxes if box is not None] + [0.0])
    return ax.transAxes + ScaledTranslation(-left_px / fig.dpi, 0, fig.dpi_scale_trans)


def build_panel_grid(plot_funcs, sizes=PANEL_SIZE_IN, ncols=2, panel_labels=True, main_title=None,
                     wspace=GRID_WSPACE, hspace=GRID_HSPACE):
    """Return a figure with one panel per plot function, lettered A, B, ..., separated by thin gray lines.

    Panels are placed by constrained layout, so colorbars, legends, and labels never overlap a neighbor,
    and each separator runs through the middle of the gap between the panels' full extents.

    Args:
        plot_funcs: Callables drawing on a given Axes, called as func(ax=ax).
        sizes: Panel size in inches, or one size per panel.
        ncols: Number of grid columns.
        panel_labels: Whether to letter the panels.
        main_title: Optional figure title.
        wspace: Horizontal gap between panels, as a fraction of the figure width.
        hspace: Vertical gap between panels, as a fraction of the figure height.
    """
    n_panels = len(plot_funcs)
    sizes = [sizes] * n_panels if isinstance(sizes, (int, float)) else list(sizes)
    n_rows = math.ceil(n_panels / ncols)
    full = sizes + [sizes[-1]] * (n_rows * ncols - n_panels)
    grid = [full[r * ncols:(r + 1) * ncols] for r in range(n_rows)]
    row_heights = [max(row) for row in grid]
    col_widths = [max(grid[r][c] for r in range(n_rows)) for c in range(ncols)]
    fig = plt.figure(figsize=(sum(col_widths), sum(row_heights)), layout="constrained")
    fig.get_layout_engine().set(wspace=wspace, hspace=hspace)
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
            ax.text(0, 1, string.ascii_uppercase[i], transform=_label_transform(ax), fontsize=GRID_LABEL_SIZE,
                    fontweight="bold", va="bottom", ha="right")
        axes.append(ax)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()
    extents = [to_fig.transform_bbox(ax.get_tightbbox(renderer)) for ax in axes]
    fig_w, fig_h = fig.get_size_inches()
    dx = GRID_GUTTER_MM / MM_PER_INCH / fig_w
    dy = GRID_GUTTER_MM / MM_PER_INCH / fig_h
    for c in range(ncols - 1):
        for r in range(n_rows):
            left, right = r * ncols + c, r * ncols + c + 1
            if right < len(axes):
                el, er = extents[left], extents[right]
                x = (el.x1 + er.x0) / 2.0
                fig.add_artist(Line2D([x, x], [min(el.y0, er.y0) + dy, max(el.y1, er.y1) - dy],
                                      transform=fig.transFigure, color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))
    for r in range(n_rows - 1):
        for c in range(ncols):
            top, bottom = r * ncols + c, (r + 1) * ncols + c
            if bottom < len(axes):
                et, eb = extents[top], extents[bottom]
                y = (et.y0 + eb.y1) / 2.0
                fig.add_artist(Line2D([min(et.x0, eb.x0) + dx, max(et.x1, eb.x1) - dx], [y, y],
                                      transform=fig.transFigure, color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))
    return fig


def export_panels_pdf(grids, outdir, filename=PANEL_PDF_NAME, raster=True, dpi=RASTER_DPI, **grid_kwargs):
    """Write one PDF page per panel grid and return the PDF path.

    Args:
        grids: Iterable of (plot_funcs, main_title) pairs, one per page.
        outdir: Output directory, created if missing.
        filename: PDF file name.
        raster: Whether to rasterize collections and images.
        dpi: Resolution of rasterized content.
        **grid_kwargs: Passed to build_panel_grid.
    """
    path = Path(outdir).expanduser() / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        for plot_funcs, main_title in grids:
            fig = build_panel_grid(plot_funcs, main_title=main_title, **grid_kwargs)
            if raster:
                rasterize_figure(fig)
            pdf.savefig(fig, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
    return path
