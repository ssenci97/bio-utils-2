#!/usr/bin/env python3
"""Render the PNG figures of the guides into docs/assets/ from the bundled 1AAY structure.

Writes contact_map_1aay.png, partners_two_plus_one.png, highlight_dim_alpha.png, chain_graph_1brs.png (from 1BRS) and, when the docking example has been
run, docking_caprieval_1aay.png. The 3D views are interactive HTML, so highlight_dim_alpha.png is a static 2D picture
of what highlight= and dim_alpha= do.

Run by hand with:  python docs/scripts/guide_images.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

from synced_structure import SyncedStructure  # noqa: E402

###################################### Configuration
ROOT = Path(__file__).resolve().parents[2]
CIF = ROOT / "tests" / "data" / "1AAY.cif"
BRS = ROOT / "tests" / "data" / "1BRS.cif"
CAPRIEVAL = ROOT / "examples" / "output" / "docking_1AAY_oracle" / "tables" / "1AAY_oracle_caprieval_models.tsv"
ASSETS = ROOT / "docs" / "assets"
DPI = 130
DIM_ALPHA = 0.25
HIGHLIGHT = {"A": range(120, 131), "B": range(3, 7), "C": range(51, 55)}
CHAIN_COLORS = {"A": "#1f77b4", "B": "#d62728", "C": "#ff7f0e"}
ANCHORS = ("CA", "C4'", "P")


def _trace(structure: SyncedStructure, key: str) -> pd.DataFrame:
    """Return residue numbers and x, y, z of one anchor atom per residue of a chain.

    Args:
        structure: The structure.
        key: Chain key.

    Returns:
        A table with columns resnum, x, y and z.
    """
    rows = []
    for residue in structure.chains[key].residues:
        atom = next((residue[name] for name in ANCHORS if name in residue), None)
        if atom is not None:
            rows.append({"resnum": residue.id[1], "x": atom.coord[0], "y": atom.coord[1], "z": atom.coord[2]})
    return pd.DataFrame(rows)


def contact_map(structure: SyncedStructure) -> Path:
    """Draw the contact map of all chains of 1AAY.

    Args:
        structure: The structure.

    Returns:
        The PNG path.
    """
    fig, _ = structure.plot_multi_contact_map(structure.chain_ids)
    path = ASSETS / "contact_map_1aay.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def highlight_picture(structure: SyncedStructure) -> Path:
    """Draw the chain traces with highlighted residues opaque and all others at DIM_ALPHA.

    Args:
        structure: The structure.

    Returns:
        The PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    for ax, dim in zip(axes, (1.0, DIM_ALPHA)):
        for key in structure.chain_ids:
            trace = _trace(structure, key)
            chosen = trace["resnum"].isin(HIGHLIGHT.get(key, ()))
            ax.plot(trace["x"], trace["y"], color=CHAIN_COLORS[key], alpha=dim, lw=1.2)
            ax.scatter(trace.loc[~chosen, "x"], trace.loc[~chosen, "y"], s=10, color=CHAIN_COLORS[key], alpha=dim)
            if dim < 1:
                ax.scatter(trace.loc[chosen, "x"], trace.loc[chosen, "y"], s=34, color=CHAIN_COLORS[key],
                           edgecolors="black", linewidths=0.6, zorder=3)
        ax.set_aspect("equal")
        ax.axis("off")
    axes[0].set_title("show_structure_3d()", family="monospace", fontsize=10)
    axes[1].set_title(f"highlight={{...}}, dim_alpha={DIM_ALPHA}", family="monospace", fontsize=10)
    fig.suptitle("1AAY chain traces (A protein, B and C DNA); static 2D picture of the 3D option", fontsize=9)
    path = ASSETS / "highlight_dim_alpha.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def partners_picture() -> Path:
    """Draw how three chains become two HADDOCK3 molecules.

    Returns:
        The PNG path.
    """
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.4, 4.9)
    ax.axis("off")
    chains = (("chain B", 3.0, "#d62728"), ("chain C", 1.9, "#ff7f0e"), ("chain A", 0.2, "#1f77b4"))
    for label, y, color in chains:
        ax.add_patch(FancyBboxPatch((0.4, y), 2.0, 0.8, boxstyle="round,pad=0.05", fc=color, alpha=0.9, ec="black"))
        ax.text(1.4, y + 0.4, label, ha="center", va="center", color="white", weight="bold")
        ax.annotate("", xy=(5.4, y + 0.4), xytext=(2.6, y + 0.4), arrowprops={"arrowstyle": "->"})
    molecules = ((1.6, 4.0, "segid A (molecule 1)\nrigid body of 2 chains", "#d62728"),
                 (-0.1, 1.2, "segid B (molecule 2)\nrigid body of 1 chain", "#1f77b4"))
    for low, high, label, color in molecules:
        ax.add_patch(FancyBboxPatch((5.5, low), 4.0, high - low, boxstyle="round,pad=0.05", fc="none", ec=color, lw=2.2))
        ax.text(7.5, (low + high) / 2, label, ha="center", va="center", fontsize=10)
    ax.text(5.0, 4.75, 'prepare_haddock_docking(["B", "C"], "A")', ha="center", family="monospace", fontsize=10)
    path = ASSETS / "partners_two_plus_one.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def chain_graph_picture() -> Path:
    """Draw the chain contact graph of 1BRS, three barnase-barstar complexes.

    Returns:
        The PNG path.
    """
    fig, _ = SyncedStructure.from_path(BRS).plot_chain_graph()
    path = ASSETS / "chain_graph_1brs.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def caprieval_picture() -> Path | None:
    """Draw HADDOCK score against i-RMSD from the docking example's caprieval table, when it exists.

    Returns:
        The PNG path, or None when the example has not been run.
    """
    if not CAPRIEVAL.is_file():
        return None
    models = pd.read_csv(CAPRIEVAL, sep="\t")
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    points = ax.scatter(models["irmsd"], models["score"], c=models["dockq"], cmap="viridis", s=60, edgecolors="black")
    ax.set_xlabel("i-RMSD (Å)")
    ax.set_ylabel("HADDOCK score")
    ax.set_title("1AAY oracle re-docking, downsized run", fontsize=10)
    fig.colorbar(points, ax=ax, label="DockQ")
    path = ASSETS / "docking_caprieval_1aay.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def render() -> list[Path]:
    """Render every figure.

    Returns:
        The PNG paths written.
    """
    structure = SyncedStructure.from_path(CIF)
    paths = [contact_map(structure), highlight_picture(structure), partners_picture(), chain_graph_picture(),
             caprieval_picture()]
    return [path for path in paths if path is not None]


if __name__ == "__main__":
    for written in render():
        print(written.relative_to(ROOT))
