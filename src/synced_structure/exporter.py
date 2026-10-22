"""Curated exporters for SyncedStructure tables and figures: TSV tables for chains, residues, moieties, alignments, and interfaces, single or grid figures, and standalone 3D view pages."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

from collections.abc import Callable, Sequence
from html import escape
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from ._types import ChainPair, ChainRefs, PathLike
from .base import StructureState
from .constants import CHAIN_COLUMNS, EXPORT_SEP, EXPORT_SUFFIX, MAPPING_COLUMNS
from .panels import build_panel_grid, rasterize_figure
from .plot_style import EXPORT_DPI, FIGURE_FORMATS, PANEL_FORMATS, PDF_RC
from .viz3d import VIZ_PAGE_TEMPLATE

if TYPE_CHECKING:
    from matplotlib.figure import Figure

###################################### EXPORT HELPERS


class ExportMixin(StructureState):
    """Write tables and figures to the output directory."""

    def _resolve_outdir(self, outdir: PathLike | None = None) -> Path:
        """Resolve the output directory, creating it if missing.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The directory.
        """
        path = Path(outdir or self.config["save_outdir"]).expanduser()
        path.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def _prepare_path(path: PathLike) -> Path:
        """Resolve a file path, creating its parent directory if missing.

        Args:
            path: File path.

        Returns:
            The path, with the user directory expanded.
        """
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _write_table(self, frame: pd.DataFrame, name: str, outdir: PathLike | None = None) -> Path:
        """Write a table to '{structure_id}_{name}.tsv' in the output directory.

        Args:
            frame: Table to write.
            name: File name suffix after the structure ID.
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path.
        """
        path = self._resolve_outdir(outdir) / f"{self.structure_id}_{name}{EXPORT_SUFFIX}"
        frame.to_csv(path, sep=EXPORT_SEP, index=False)
        return path

    def _save_figure(self, fig: Figure, stem: str, outdir: PathLike | None, formats: Sequence[str],
                     dpi: int = EXPORT_DPI, rasterize_pdf: bool = True) -> list[Path]:
        """Save a figure at its exact size in each format, then close it.

        Args:
            fig: Figure to save.
            stem: File name without extension.
            outdir: Output directory, or None for the configured save_outdir.
            formats: File formats, such as 'pdf' or 'png'.
            dpi: Resolution of raster output.
            rasterize_pdf: Whether to rasterize collections and images in PDFs.

        Returns:
            The written paths.
        """
        base = self._resolve_outdir(outdir) / stem
        paths = []
        for fmt in formats:
            path = base.with_name(f"{base.name}.{fmt}")
            if fmt == "pdf" and rasterize_pdf:
                rasterize_figure(fig)
            # PDF_RC holds valid rc keys; matplotlib types them as a Literal union that a plain dict can't match
            with mpl.rc_context(cast("dict[Any, Any]", PDF_RC)):
                fig.savefig(path, dpi=dpi)
            paths.append(path)
        plt.close(fig)
        return paths

###################################### CHAIN AND RESIDUE TABLES

    def export_chain_table(self, outdir: PathLike | None = None) -> Path:
        """Export one row per chain with entity, description, and polymer residue count.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_chains.tsv'.
        """
        descriptions = self._entity_descriptions()
        rows = [
            {
                "structure_id": self.structure_id,
                "chain_id": data.chain_id,
                "auth_chain_id": data.auth_chain_id,
                "entity_id": data.entity_id,
                "description": descriptions.get(data.entity_id or ""),
                "residues": len(data.residues),
            }
            for data in self.chains.values()
        ]
        return self._write_table(pd.DataFrame(rows, columns=CHAIN_COLUMNS), "chains", outdir)

    def export_chain_residue_table(self, chain_ids: ChainRefs | None = None, outdir: PathLike | None = None) -> Path:
        """Export the residue-to-reference mapping of the selected chains.

        Args:
            chain_ids: Chains to export, or None for all.
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_residue_mapping.tsv'.
        """
        selected = self.chain_ids if chain_ids is None else [self._chain_key(c) for c in chain_ids]
        frames = [self.chain_residue_table(c) for c in tqdm(selected, desc="mapping", unit="chain", leave=False)]
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=MAPPING_COLUMNS)
        return self._write_table(frame, "residue_mapping", outdir)

    def export_residue_table(self, outdir: PathLike | None = None) -> Path:
        """Export residue_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_residues.tsv'.
        """
        return self._write_table(self.residue_table, "residues", outdir)

    def export_moiety_table(self, outdir: PathLike | None = None) -> Path:
        """Export moiety_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_moieties.tsv'.
        """
        return self._write_table(self.moiety_table, "moieties", outdir)

###################################### ALIGNMENT TABLES

    def export_alignment_table(self, outdir: PathLike | None = None) -> Path:
        """Export alignment_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_alignment_stats.tsv'.
        """
        return self._write_table(self.alignment_table, "alignment_stats", outdir)

    def export_alignment_column_table(self, outdir: PathLike | None = None) -> Path:
        """Export alignment_column_table.

        Args:
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written path, '{structure_id}_alignment_positions.tsv'.
        """
        return self._write_table(self.alignment_column_table, "alignment_positions", outdir)

###################################### INTERFACE TABLES

    def export_interface_tables(self, chain_pairs: Sequence[ChainPair] | None = None,
                                outdir: PathLike | None = None) -> list[Path]:
        """Export complex-level and residue-level interface metrics.

        Args:
            chain_pairs: Chain pairs to export, or None for all pairs.
            outdir: Output directory, or None for the configured save_outdir.

        Returns:
            The written paths, '{structure_id}_interface_complex.tsv' and '{structure_id}_interface_residues.tsv'.
        """
        tables = self.interface_tables(chain_pairs)
        return [self._write_table(tables["complex"], "interface_complex", outdir),
                self._write_table(tables["residues"], "interface_residues", outdir)]

###################################### FIGURES

    def export_plot(self, plot_func: Callable[[], tuple[Figure, Any]], name: str, outdir: PathLike | None = None,
                    formats: Sequence[str] = FIGURE_FORMATS, dpi: int = EXPORT_DPI, raster: bool = True) -> list[Path]:
        """Save the figure returned by plot_func() as '{structure_id}_{name}' in each format.

        Args:
            plot_func: Callable returning (figure, axes).
            name: File name suffix after the structure ID.
            outdir: Output directory, or None for the configured save_outdir.
            formats: File formats.
            dpi: Resolution of raster output.
            raster: Whether to rasterize collections and images in PDFs.

        Returns:
            The written paths.
        """
        fig, _ = plot_func()
        return self._save_figure(fig, f"{self.structure_id}_{name}", outdir, formats, dpi, raster)

    def export_plot_grid(self, plot_funcs: Sequence[Callable[..., Any]], name: str, outdir: PathLike | None = None,
                         formats: Sequence[str] = PANEL_FORMATS, dpi: int = EXPORT_DPI, raster: bool = True,
                         **grid_kwargs: Any) -> list[Path]:
        """Save a panel grid of plot functions as 'name' in each format.

        Args:
            plot_funcs: Plot functions, one per panel, as taken by build_panel_grid.
            name: File name without extension.
            outdir: Output directory, or None for the configured save_outdir.
            formats: File formats.
            dpi: Resolution of raster output.
            raster: Whether to rasterize collections and images in PDFs.
            **grid_kwargs: Passed to build_panel_grid.

        Returns:
            The written paths.
        """
        fig = build_panel_grid(plot_funcs, **grid_kwargs)
        return self._save_figure(fig, name, outdir, formats, dpi, raster)

###################################### 3D VIEWS

    def export_structure_3d(self, name: str = "structure_3d", outdir: PathLike | None = None, title: str | None = None,
                            **view_kwargs: Any) -> Path:
        """Save show_structure_3d as a standalone HTML page '{structure_id}_{name}.html'.

        Args:
            name: File name suffix after the structure ID.
            outdir: Output directory, or None for the configured save_outdir.
            title: Optional title shown above the view and used as the page title.
            **view_kwargs: Passed to show_structure_3d.

        Returns:
            The written path.
        """
        view = self.show_structure_3d(title=title, **view_kwargs)
        page = VIZ_PAGE_TEMPLATE.replace("__TITLE__", escape(str(title or self.structure_id)))
        path = self._resolve_outdir(outdir) / f"{self.structure_id}_{name}.html"
        path.write_text(page.replace("__BODY__", view.data), encoding="utf-8")
        return path
