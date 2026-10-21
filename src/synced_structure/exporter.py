"""Curated exporters for SyncedStructure tables and figures: TSV tables for chains, residues, moieties, alignments, and interfaces, single or grid figures, and standalone 3D view pages."""

###################################### IMPORTS AND CONSTANTS

from html import escape
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from .constants import CHAIN_COLUMNS, EXPORT_SEP, EXPORT_SUFFIX, MAPPING_COLUMNS
from .panels import rasterize_figure, build_panel_grid
from .plot_style import EXPORT_DPI, FIGURE_FORMATS, PANEL_FORMATS, PDF_RC
from .viz3d import VIZ_PAGE_TEMPLATE

###################################### EXPORT HELPERS


class ExportMixin:
    """Write tables and figures to the output directory."""

    def _outdir(self, outdir=None):
        """Return the output directory, the configured save_outdir by default, creating it if missing."""
        path = Path(outdir or self.config["save_outdir"]).expanduser()
        path.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def _prepare_path(path):
        """Return a file path as a Path, creating its parent directory if missing."""
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _write_table(self, frame, name, outdir=None):
        """Write a DataFrame to '{structure_id}_{name}.tsv' in the output directory and return its path."""
        path = self._outdir(outdir) / f"{self.structure_id}_{name}{EXPORT_SUFFIX}"
        frame.to_csv(path, sep=EXPORT_SEP, index=False)
        return path

    def _save_figure(self, fig, stem, outdir, formats, dpi=EXPORT_DPI, rasterize_pdf=True):
        """Save a figure in each format at its exact figure size, rasterizing PDF collections if raster is True, close it,
        and return the paths."""
        base = self._outdir(outdir) / stem
        paths = []
        for fmt in formats:
            path = base.with_name(f"{base.name}.{fmt}")
            if fmt == "pdf" and rasterize_pdf:
                rasterize_figure(fig)
            with mpl.rc_context(PDF_RC):
                fig.savefig(path, dpi=dpi)
            paths.append(path)
        plt.close(fig)
        return paths

###################################### CHAIN AND RESIDUE TABLES

    def export_chains(self, outdir=None):
        """Export one row per chain with entity, description, and polymer residue count."""
        descriptions = self._entity_descriptions()
        rows = [
            {
                "structure_id": self.structure_id,
                "chain_id": data.chain_id,
                "auth_chain_id": data.auth_chain_id,
                "entity_id": data.entity_id,
                "description": descriptions.get(data.entity_id),
                "residues": len(data.residues),
            }
            for data in self.chains.values()
        ]
        return self._write_table(pd.DataFrame(rows, columns=CHAIN_COLUMNS), "chains", outdir)

    def export_residue_mapping(self, chain_ids=None, outdir=None):
        """Export the residue-to-reference mapping of the selected chains, or all chains."""
        selected = self.chain_ids if chain_ids is None else [self._require_chain(c) for c in chain_ids]
        frames = [self.get_residue_mapping(c) for c in tqdm(selected, desc="mapping", unit="chain", leave=False)]
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=MAPPING_COLUMNS)
        return self._write_table(frame, "residue_mapping", outdir)

    def export_residues(self, outdir=None):
        """Export residue_map."""
        return self._write_table(self.residue_map, "residues", outdir)

    def export_moieties(self, outdir=None):
        """Export moiety_summary."""
        return self._write_table(self.moiety_summary, "moieties", outdir)

###################################### ALIGNMENT TABLES

    def export_alignment_stats(self, outdir=None):
        """Export alignment_stats."""
        return self._write_table(self.alignment_stats, "alignment_stats", outdir)

    def export_alignment_positions(self, outdir=None):
        """Export alignment_positions."""
        return self._write_table(self.alignment_positions, "alignment_positions", outdir)

###################################### INTERFACE TABLES

    def export_interface_complex(self, chain_pairs=None, outdir=None):
        """Export complex-level interface metrics for the given chain pairs, or all pairs."""
        return self._write_table(self.interface_all(chain_pairs)["complex"], "interface_complex", outdir)

    def export_interface_residues(self, chain_pairs=None, outdir=None):
        """Export residue-level interface metrics for the given chain pairs, or all pairs."""
        return self._write_table(self.interface_all(chain_pairs)["residues"], "interface_residues", outdir)

###################################### FIGURES

    def export_plot(self, plot_func, name, outdir=None, formats=FIGURE_FORMATS, dpi=EXPORT_DPI, raster=True):
        """Save the figure returned by plot_func() as '{structure_id}_{name}' in each format and return the paths."""
        fig, _ = plot_func()
        return self._save_figure(fig, f"{self.structure_id}_{name}", outdir, formats, dpi, raster)

    def export_plot_grid(self, plot_funcs, name, outdir=None, formats=PANEL_FORMATS, dpi=EXPORT_DPI,
                         raster=True, **grid_kwargs):
        """Save a panel grid of plot functions as 'name' in each format and return the paths."""
        fig = build_panel_grid(plot_funcs, **grid_kwargs)
        return self._save_figure(fig, name, outdir, formats, dpi, raster)

###################################### 3D VIEWS

    def export_structure_3d(self, name="structure_3d", outdir=None, title=None, **view_kwargs):
        """Save show_structure_3d as a standalone HTML page '{structure_id}_{name}.html' and return its path.

        Args:
            name: File name suffix after the structure ID.
            outdir: Output directory, the configured save_outdir by default.
            title: Optional title shown above the view and used as the page title.
            **view_kwargs: Passed to show_structure_3d.
        """
        view = self.show_structure_3d(title=title, **view_kwargs)
        page = VIZ_PAGE_TEMPLATE.replace("__TITLE__", escape(str(title or self.structure_id)))
        path = self._outdir(outdir) / f"{self.structure_id}_{name}.html"
        path.write_text(page.replace("__BODY__", view.data), encoding="utf-8")
        return path
