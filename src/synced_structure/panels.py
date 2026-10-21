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

from .plot_style import (
    GRID_GUTTER_MM, GRID_LABEL_SIZE, GRID_LINE_COLOR, GRID_LINE_WIDTH, GRID_SHARED_GAP_OFFSET, GRID_TITLE_SIZE,
    MM_PER_INCH, PANEL_LABEL_X, PANEL_LABEL_Y, PANEL_PDF_NAME, RASTER_DPI,
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


def build_panel_grid(plot_funcs, sizes=3.5, ncols=2, panel_labels=True, main_title=None, wspace=0.5, hspace=0.5):
    """Return a figure with one panel per plot function, lettered A, B, ..., separated by thin gray lines.

    Args:
        plot_funcs: Callables drawing on a given Axes, called as func(ax=ax).
        sizes: Panel size in inches, or one size per panel.
        ncols: Number of grid columns.
        panel_labels: Whether to letter the panels.
        main_title: Optional figure title.
        wspace: Horizontal spacing between panels.
        hspace: Vertical spacing between panels.
    """
    n_panels = len(plot_funcs)
    sizes = [sizes] * n_panels if isinstance(sizes, (int, float)) else list(sizes)
    n_rows = math.ceil(n_panels / ncols)
    full = sizes + [sizes[-1]] * (n_rows * ncols - n_panels)
    grid = [full[r * ncols:(r + 1) * ncols] for r in range(n_rows)]
    row_heights = [max(row) for row in grid]
    col_widths = [max(grid[r][c] for r in range(n_rows)) for c in range(ncols)]
    fig = plt.figure(figsize=(sum(col_widths), sum(row_heights)))
    gs = fig.add_gridspec(n_rows, ncols, width_ratios=col_widths, height_ratios=row_heights,
                          wspace=wspace, hspace=hspace)
    if main_title:
        fig.suptitle(main_title, fontsize=GRID_TITLE_SIZE, fontweight="bold", y=0.95)
    axes = []
    for i, func in enumerate(plot_funcs):
        row, col = divmod(i, ncols)
        ax = fig.add_subplot(gs[row, col])
        func(ax=ax)
        if panel_labels and i >= len(string.ascii_uppercase):
            log_info(f"Only {len(string.ascii_uppercase)} panels can be lettered; panel {i + 1} is unlabeled.")
        if panel_labels and i < len(string.ascii_uppercase):
            ax.text(PANEL_LABEL_X, PANEL_LABEL_Y, string.ascii_uppercase[i], transform=ax.transAxes, fontsize=GRID_LABEL_SIZE,
                    fontweight="bold", va="bottom", ha="right")
        axes.append(ax)
    fig.canvas.draw()
    fig_w, fig_h = fig.get_size_inches()
    dx = GRID_GUTTER_MM / MM_PER_INCH / fig_w
    dy = GRID_GUTTER_MM / MM_PER_INCH / fig_h
    for c in range(ncols - 1):
        for r in range(n_rows):
            left, right = r * ncols + c, r * ncols + c + 1
            if right < len(axes):
                pl, pr = axes[left].get_position(), axes[right].get_position()
                x = (pl.x1 + pr.x0) / 2.0 - GRID_SHARED_GAP_OFFSET
                fig.add_artist(Line2D([x, x], [pl.y0 + dy, pl.y1 - dy], transform=fig.transFigure,
                                      color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))
    for r in range(n_rows - 1):
        for c in range(ncols):
            top, bottom = r * ncols + c, (r + 1) * ncols + c
            if bottom < len(axes):
                pt, pb = axes[top].get_position(), axes[bottom].get_position()
                y = (pt.y0 + pb.y1) / 2.0
                fig.add_artist(Line2D([pt.x0 + dx, pt.x1 - dx], [y, y], transform=fig.transFigure,
                                      color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))
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
