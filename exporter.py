#!/usr/bin/env python
"""Curated exporters for SyncedStructure tables and figures: CSV tables for chains, residues, moieties, alignments, and interfaces, single figures, and multi-page panel PDFs with lettered panels and optional high-resolution rasterization."""

###################################### IMPORTS AND CONSTANTS

import math
import string
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import Collection
from matplotlib.image import AxesImage
from matplotlib.lines import Line2D
from tqdm import tqdm

from constants import CHAIN_COLUMNS, EXPORT_SUFFIX, MAPPING_COLUMNS

EXPORT_DPI = 300
RASTER_DPI = 600
FIGURE_FORMATS = ("png", "svg")
PANEL_FORMATS = ("pdf", "png")
PANEL_PDF_NAME = "panels.pdf"
GRID_LINE_COLOR = "gray"
GRID_LINE_WIDTH = 0.5
GRID_GUTTER_MM = 0.2
GRID_SHARED_GAP_OFFSET = 0.004
GRID_LABEL_SIZE = 14
GRID_TITLE_SIZE = 16

###################################### RASTERIZATION


def _rasterize(fig):
    for ax in fig.axes:
        for artist in ax.get_children():
            if isinstance(artist, (Collection, AxesImage)):
                artist.set_rasterized(True)
    return fig

###################################### PANEL LAYOUT


def build_panel_grid(plot_funcs, sizes=3.5, ncols=2, panel_labels=True, main_title=None, wspace=0.5, hspace=0.5):
    n = len(plot_funcs)
    sizes = [sizes] * n if isinstance(sizes, (int, float)) else list(sizes)
    nrows = math.ceil(n / ncols)
    full = sizes + [sizes[-1]] * (nrows * ncols - n)
    grid = [full[r * ncols:(r + 1) * ncols] for r in range(nrows)]
    row_heights = [max(row) for row in grid]
    col_widths = [max(grid[r][c] for r in range(nrows)) for c in range(ncols)]
    fig = plt.figure(figsize=(sum(col_widths), sum(row_heights)))
    gs = fig.add_gridspec(nrows, ncols, width_ratios=col_widths, height_ratios=row_heights,
                          wspace=wspace, hspace=hspace)
    if main_title:
        fig.suptitle(main_title, fontsize=GRID_TITLE_SIZE, fontweight="bold", y=0.95)
    axes = []
    for i, func in enumerate(plot_funcs):
        r, c = divmod(i, ncols)
        ax = fig.add_subplot(gs[r, c])
        func(ax=ax)
        if panel_labels and i < len(string.ascii_uppercase):
            ax.text(-0.12, 1.05, string.ascii_uppercase[i], transform=ax.transAxes, fontsize=GRID_LABEL_SIZE,
                    fontweight="bold", va="bottom", ha="right")
        axes.append(ax)
    fig.canvas.draw()
    fig_w, fig_h = fig.get_size_inches()
    dx = GRID_GUTTER_MM / 25.4 / fig_w
    dy = GRID_GUTTER_MM / 25.4 / fig_h
    for c in range(ncols - 1):
        for r in range(nrows):
            left, right = r * ncols + c, r * ncols + c + 1
            if right < len(axes):
                pl, pr = axes[left].get_position(), axes[right].get_position()
                x = (pl.x1 + pr.x0) / 2.0 - GRID_SHARED_GAP_OFFSET
                fig.add_artist(Line2D([x, x], [pl.y0 + dy, pl.y1 - dy], transform=fig.transFigure,
                                      color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))
    for r in range(nrows - 1):
        for c in range(ncols):
            top, bottom = r * ncols + c, (r + 1) * ncols + c
            if bottom < len(axes):
                pt, pb = axes[top].get_position(), axes[bottom].get_position()
                y = (pt.y0 + pb.y1) / 2.0
                fig.add_artist(Line2D([pt.x0 + dx, pt.x1 - dx], [y, y], transform=fig.transFigure,
                                      color=GRID_LINE_COLOR, linewidth=GRID_LINE_WIDTH))
    return fig


def export_panels_pdf(grids, outdir, filename=PANEL_PDF_NAME, raster=True, dpi=RASTER_DPI, **grid_kwargs):
    path = Path(outdir) / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        for plot_funcs, main_title in grids:
            fig = build_panel_grid(plot_funcs, main_title=main_title, **grid_kwargs)
            if raster:
                _rasterize(fig)
            pdf.savefig(fig, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
    return path

###################################### EXPORT HELPERS


class ExportMixin:
    def _outdir(self, outdir=None):
        path = Path(outdir or self.config["save_outdir"])
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _export(self, frame, name, outdir=None):
        path = self._outdir(outdir) / f"{self.structure_id}_{name}{EXPORT_SUFFIX}"
        frame.to_csv(path, index=False)
        return path

    def _save_figure(self, fig, stem, outdir, formats, dpi=EXPORT_DPI, raster=True):
        base = self._outdir(outdir) / stem
        paths = []
        for fmt in formats:
            path = base.with_name(f"{base.name}.{fmt}")
            if fmt == "pdf" and raster:
                _rasterize(fig)
            fig.savefig(path, dpi=RASTER_DPI if fmt == "pdf" and raster else dpi, bbox_inches="tight")
            paths.append(path)
        plt.close(fig)
        return paths

###################################### CHAIN AND RESIDUE TABLES

    def export_chains(self, outdir=None):
        descriptions = self._chain_description_map()
        rows = [
            {
                "structure_id": self.structure_id,
                "chain_id": data["chain_id"],
                "auth_chain_id": data["auth_chain_id"],
                "entity_id": data.get("entity_id"),
                "description": descriptions.get(key),
                "residues": len(data["residues"]),
            }
            for key, data in self.chains.items()
        ]
        return self._export(pd.DataFrame(rows, columns=CHAIN_COLUMNS), "chains", outdir)

    def export_residue_mapping(self, chain_ids=None, outdir=None):
        selected = self.chain_ids if chain_ids is None else [self._require_chain(c) for c in chain_ids]
        frames = [self.get_residue_mapping(c) for c in tqdm(selected, desc="mapping", unit="chain", leave=False)]
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=MAPPING_COLUMNS)
        return self._export(frame, "residue_mapping", outdir)

    def export_residues(self, outdir=None):
        return self._export(self.residue_map, "residues", outdir)

    def export_moieties(self, outdir=None):
        return self._export(self.moiety_summary, "moieties", outdir)

###################################### ALIGNMENT TABLES

    def export_alignment_stats(self, outdir=None):
        return self._export(self.alignment_stats, "alignment_stats", outdir)

    def export_alignment_positions(self, outdir=None):
        return self._export(self.alignment_positions, "alignment_positions", outdir)

###################################### INTERFACE TABLES

    def export_interface_complex(self, chain_pairs=None, outdir=None):
        return self._export(self.interface_all(chain_pairs)["complex"], "interface_complex", outdir)

    def export_interface_residues(self, chain_pairs=None, outdir=None):
        return self._export(self.interface_all(chain_pairs)["residues"], "interface_residues", outdir)

###################################### FIGURES

    def export_plot(self, plot_func, name, outdir=None, formats=FIGURE_FORMATS, dpi=EXPORT_DPI, raster=True):
        fig, _ = plot_func()
        return self._save_figure(fig, f"{self.structure_id}_{name}", outdir, formats, dpi, raster)

    def export_plot_grid(self, plot_funcs, name, outdir=None, formats=PANEL_FORMATS, dpi=EXPORT_DPI,
                         raster=True, **grid_kwargs):
        fig = build_panel_grid(plot_funcs, **grid_kwargs)
        return self._save_figure(fig, name, outdir, formats, dpi, raster)