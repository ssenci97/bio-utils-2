"""Synchronized Biopython structures: chain mapping, contacts, contact maps, interfaces, and 3D views."""

from .panels import build_panel_grid, export_panels_pdf
from .structure import SyncedStructure
from .viz3d import side_by_side

__version__ = "0.1.0"
__all__ = ["SyncedStructure", "build_panel_grid", "export_panels_pdf", "side_by_side"]
