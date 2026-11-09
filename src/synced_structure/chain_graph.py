"""Chain contact graphs: chains as nodes, joined by edges as thick and dark as the number of residue contacts between them."""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

from collections.abc import Sequence
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from ._types import Moieties
from .base import StructureState
from .constants import ANGSTROM_PER_NM, CONTACTS_CUTOFF_NM, DEFAULT_CONTACT_MOIETIES
from .plot_style import (
    CHAIN_GRAPH_COUNT_BOX, CHAIN_GRAPH_DESCRIPTION_CHARS, CHAIN_GRAPH_EDGE_ALPHA, CHAIN_GRAPH_EDGE_COLOR,
    CHAIN_GRAPH_EDGE_LW, CHAIN_GRAPH_LABELS, CHAIN_GRAPH_LAYOUTS, CHAIN_GRAPH_NODE_EDGE_COLOR, CHAIN_GRAPH_NODE_EDGE_LW,
    CHAIN_GRAPH_NODE_PT2, CHAIN_GRAPH_PAD, CHAIN_GRAPH_SIZE_IN, CHAIN_GRAPH_SPRING_ITERATIONS, CHAIN_GRAPH_TEXT_COLOR,
    CONTACT_MAP_DPI,
)
from .utils import blocked_distances
from .viz3d import CHAIN_PALETTE

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure, SubFigure

CHAIN_GRAPH_COLUMNS = ["chain_key_a", "chain_a", "auth_chain_a", "chain_key_b", "chain_b", "auth_chain_b",
                       "contact_pairs"]
# distances held at once when counting contacts between two chains, bounding memory to a few MB
CHAIN_GRAPH_BLOCK_VALUES = 200_000

###################################### LAYOUTS


def _circular_layout(n: int) -> np.ndarray:
    """Return n points evenly spaced on the unit circle, the first at the top.

    Args:
        n: Number of points.

    Returns:
        (n, 2) positions.
    """
    if n == 1:
        return np.zeros((1, 2))
    angles = np.pi / 2 - 2 * np.pi * np.arange(n) / n
    return np.column_stack([np.cos(angles), np.sin(angles)])


def _spring_layout(weights: np.ndarray, iterations: int = CHAIN_GRAPH_SPRING_ITERATIONS) -> np.ndarray:
    """Return a deterministic force-directed (Fruchterman-Reingold) layout; every edge pulls alike, so edge lengths stay even.

    Args:
        weights: (n, n) symmetric edge weights, 0 for no edge.
        iterations: Number of cooling steps.

    Returns:
        (n, 2) positions within the unit disc.
    """
    n = len(weights)
    pos = _circular_layout(n)
    if n < 3:
        return pos
    attraction = (weights > 0).astype(float)
    k = 1 / np.sqrt(n)
    for step in range(iterations):
        delta = pos[:, None] - pos[None]
        dist = np.maximum(np.linalg.norm(delta, axis=-1), 1e-3)
        magnitude = k * k / dist - attraction * dist * dist / k
        np.fill_diagonal(magnitude, 0.0)
        force = (magnitude / dist)[..., None] * delta
        force = force.sum(axis=1) - 0.1 * k * pos
        length = np.maximum(np.linalg.norm(force, axis=1), 1e-9)
        temperature = 0.1 * (1 - step / iterations)
        pos = pos + force / length[:, None] * np.minimum(length, temperature)[:, None]
    pos -= pos.mean(axis=0)
    radius = np.linalg.norm(pos, axis=1).max()
    return pos / radius if radius > 0 else pos


def _scale(values: np.ndarray, bounds: tuple[float, float]) -> np.ndarray:
    """Map values linearly onto bounds, the smallest to the lower bound.

    Args:
        values: Values to map.
        bounds: (lower, upper) output range.

    Returns:
        The mapped values; all equal values map to the upper bound.
    """
    low, high = bounds
    if not len(values):
        return values.astype(float)
    span = values.max() - values.min()
    if span == 0:
        return np.full(len(values), high, float)
    return low + (values - values.min()) / span * (high - low)

###################################### CHAIN GRAPH MIXIN


class ChainGraphMixin(StructureState):
    """Graphs of chains joined by their inter-chain contacts."""

    def chain_graph_table(self, cutoff: float = CONTACTS_CUTOFF_NM, scheme: str = "ca",
                          moieties: Moieties = DEFAULT_CONTACT_MOIETIES,
                          contacts: pd.DataFrame | None = None) -> pd.DataFrame:
        """Return the number of residue contacts between each pair of chains in contact.

        Args:
            cutoff: Contact cutoff in nm.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            moieties: Biomolecule types to include.
            contacts: Precomputed contact_table(), or None to compute inter-chain contacts.

        Returns:
            One row per chain pair with CHAIN_GRAPH_COLUMNS, most contacts first.
        """
        order = {key: index for index, key in enumerate(self.chains)}
        counts: dict[tuple[str, str], int] = {}
        if contacts is None and scheme == "ca":
            counts = self._point_contact_counts(cutoff, moieties)
        else:
            if contacts is None:
                contacts = self.contact_table(cutoff=cutoff, contact_type="inter", scheme=scheme, moieties=moieties)
            inter = contacts[contacts["chain_key_a"] != contacts["chain_key_b"]]
            for key_a, key_b in zip(inter["chain_key_a"], inter["chain_key_b"]):
                pair = (key_a, key_b) if order[key_a] <= order[key_b] else (key_b, key_a)
                counts[pair] = counts.get(pair, 0) + 1
        if not counts:
            return pd.DataFrame(columns=CHAIN_GRAPH_COLUMNS)
        rows = [{"chain_key_a": a, "chain_a": self.chains[a].chain_id, "auth_chain_a": self.chains[a].auth_chain_id,
                 "chain_key_b": b, "chain_b": self.chains[b].chain_id, "auth_chain_b": self.chains[b].auth_chain_id,
                 "contact_pairs": n} for (a, b), n in counts.items()]
        table = pd.DataFrame(rows, columns=CHAIN_GRAPH_COLUMNS)
        return table.sort_values("contact_pairs", ascending=False, kind="stable").reset_index(drop=True)

    def _point_contact_counts(self, cutoff: float, moieties: Moieties) -> dict[tuple[str, str], int]:
        """Count contact-point pairs between chains, as contact_table(scheme='ca') would, in bounded memory.

        Chain pairs whose bounding boxes lie farther apart than the cutoff are skipped, and distances are computed a
        block of rows at a time, so large assemblies need neither all residue pairs nor much memory.

        Args:
            cutoff: Contact cutoff in nm.
            moieties: Biomolecule types to include.

        Returns:
            {(chain key, chain key) in chain order: number of residue pairs within the cutoff}, pairs in contact only.
        """
        cutoff_a = cutoff * ANGSTROM_PER_NM
        points: dict[str, list[np.ndarray]] = {}
        for record in self._contact_records(moieties=tuple(moieties)):
            points.setdefault(record["key"], []).append(record["point"])
        coords = {key: np.asarray(values, np.float32) for key, values in points.items()}
        boxes = {key: (xyz.min(axis=0), xyz.max(axis=0)) for key, xyz in coords.items()}
        keys = list(coords)
        counts: dict[tuple[str, str], int] = {}
        for i, key_a in enumerate(keys):
            for key_b in keys[i + 1:]:
                (low_a, high_a), (low_b, high_b) = boxes[key_a], boxes[key_b]
                gap = np.maximum(0, np.maximum(low_a - high_b, low_b - high_a))
                if float(np.linalg.norm(gap)) > cutoff_a:
                    continue
                xa, xb = coords[key_a], coords[key_b]
                block = max(1, CHAIN_GRAPH_BLOCK_VALUES // len(xb))
                n = sum(int((dist <= cutoff_a).sum()) for _, dist in blocked_distances(xa, xb, block))
                if n:
                    counts[(key_a, key_b)] = n
        return counts

    def _node_label(self, key: str, labels: str, describe: bool, descriptions: dict[str, str]) -> tuple[str, str]:
        """Return the label of a chain node and the entity description shown under it.

        Args:
            key: Chain key.
            labels: 'key', 'auth', or 'label'.
            describe: Whether to add the entity description.
            descriptions: Entity descriptions by entity ID.

        Returns:
            (label, description or '').
        """
        record = self.chains[key]
        label = {"key": key, "auth": record.auth_chain_id, "label": record.chain_id}[labels]
        text = descriptions.get(record.entity_id or "", "") if describe else ""
        if len(text) > CHAIN_GRAPH_DESCRIPTION_CHARS:
            text = text[:CHAIN_GRAPH_DESCRIPTION_CHARS - 1].rstrip() + "…"
        return label, text

    def plot_chain_graph(self, cutoff: float = CONTACTS_CUTOFF_NM, scheme: str = "ca",
                         moieties: Moieties = DEFAULT_CONTACT_MOIETIES, layout: str = "spring",
                         min_contacts: int = 1, labels: str = "key", describe: bool = True,
                         show_counts: bool = True, colors: dict[str, str] | None = None,
                         chain_palette: Sequence[str] = CHAIN_PALETTE, size_in: float = CHAIN_GRAPH_SIZE_IN,
                         dpi: float = CONTACT_MAP_DPI, title: str | None = None, ax: Axes | None = None,
                         table: pd.DataFrame | None = None) -> tuple[Figure | SubFigure, Axes]:
        """Plot the chains as nodes joined by edges as thick and dark as the number of contacts between them.

        Node colors match show_structure_3d; node area grows with the number of residues of the chain. Chains
        without contacts are drawn unconnected. For a PDB entry: SyncedStructure.from_rcsb('6R3Q').plot_chain_graph().

        Args:
            cutoff: Contact cutoff in nm.
            scheme: Contact scheme, 'ca' or 'sidechain'.
            moieties: Biomolecule types to include.
            layout: 'spring' (chains in contact drawn close) or 'circular'.
            min_contacts: Fewest contacts drawn as an edge.
            labels: Chain ID shown in each node: 'key', 'auth', or 'label'.
            describe: Whether to show the entity description under each node.
            show_counts: Whether to write the number of contacts on each edge.
            colors: Chain ID to color; other chains take colors from chain_palette.
            chain_palette: Colors for chains without an explicit color, in order.
            size_in: Figure width and height in inches.
            dpi: Figure resolution on screen; export_plot sets the saved resolution.
            title: Title, or None for a default.
            ax: Axes to draw on, or None to create a figure.
            table: Precomputed chain_graph_table(), or None to compute it.

        Returns:
            (figure, axes).

        Raises:
            ValueError: If layout or labels is unknown.
        """
        if layout not in CHAIN_GRAPH_LAYOUTS:
            raise ValueError(f"layout must be one of {CHAIN_GRAPH_LAYOUTS}.")
        if labels not in CHAIN_GRAPH_LABELS:
            raise ValueError(f"labels must be one of {CHAIN_GRAPH_LABELS}.")
        if table is None:
            table = self.chain_graph_table(cutoff=cutoff, scheme=scheme, moieties=moieties)
        table = table[table["contact_pairs"] >= min_contacts]
        keys = list(self.chains)
        index = {key: i for i, key in enumerate(keys)}
        weights = np.zeros((len(keys), len(keys)))
        for key_a, key_b, count in zip(table["chain_key_a"], table["chain_key_b"], table["contact_pairs"]):
            weights[index[key_a], index[key_b]] = weights[index[key_b], index[key_a]] = count
        pos = _spring_layout(weights) if layout == "spring" else _circular_layout(len(keys))

        fig: Figure | SubFigure
        if ax is None:
            fig, ax = plt.subplots(figsize=(size_in, size_in), dpi=dpi, layout="constrained")
        else:
            fig = ax.figure
        edges = table.sort_values("contact_pairs", kind="stable")
        counts = edges["contact_pairs"].to_numpy(float)
        widths, alphas = _scale(counts, CHAIN_GRAPH_EDGE_LW), _scale(counts, CHAIN_GRAPH_EDGE_ALPHA)
        for (key_a, key_b, count), width, alpha in zip(
                edges[["chain_key_a", "chain_key_b", "contact_pairs"]].itertuples(index=False), widths, alphas):
            (xa, ya), (xb, yb) = pos[index[key_a]], pos[index[key_b]]
            ax.plot([xa, xb], [ya, yb], color=CHAIN_GRAPH_EDGE_COLOR, lw=width, alpha=alpha,
                    solid_capstyle="round", zorder=1)
            if show_counts:
                ax.text((xa + xb) / 2, (ya + yb) / 2, str(int(count)), ha="center", va="center", fontsize=7,
                        color=CHAIN_GRAPH_TEXT_COLOR, bbox=CHAIN_GRAPH_COUNT_BOX, zorder=2)

        palette = self._chain_colors(colors, chain_palette)
        sizes = _scale(np.sqrt([len(self.chains[key].residues) for key in keys]), CHAIN_GRAPH_NODE_PT2)
        ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=[palette[key] for key in keys], zorder=3,
                   edgecolors=CHAIN_GRAPH_NODE_EDGE_COLOR, linewidths=CHAIN_GRAPH_NODE_EDGE_LW)
        descriptions = self._entity_descriptions() if describe else {}
        for key, (x, y), size in zip(keys, pos, sizes):
            label, text = self._node_label(key, labels, describe, descriptions)
            ax.text(x, y, label, ha="center", va="center", fontsize=9, fontweight="bold",
                    color=CHAIN_GRAPH_TEXT_COLOR, zorder=4)
            if text:
                ax.annotate(text, (x, y), xytext=(0, -np.sqrt(size) / 2 - 2), textcoords="offset points",
                            ha="center", va="top", fontsize=6.5, color=CHAIN_GRAPH_TEXT_COLOR,
                            bbox=CHAIN_GRAPH_COUNT_BOX, zorder=4)

        if len(counts):
            shown = sorted({int(counts.min()), int(counts.max())})
            styles = list(zip(CHAIN_GRAPH_EDGE_LW, CHAIN_GRAPH_EDGE_ALPHA))[-len(shown):]
            handles = [Line2D([], [], color=CHAIN_GRAPH_EDGE_COLOR, lw=lw, alpha=alpha, solid_capstyle="round")
                       for lw, alpha in styles]
            ax.legend(handles, [str(c) for c in shown], title="contacts", loc="upper center",
                      bbox_to_anchor=(0.5, 0.0), ncol=len(shown), frameon=False, fontsize=7, title_fontsize=7)
        ax.set_xlim(-1 - CHAIN_GRAPH_PAD, 1 + CHAIN_GRAPH_PAD)
        ax.set_ylim(-1 - CHAIN_GRAPH_PAD, 1 + CHAIN_GRAPH_PAD)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title or f"{self.structure_id} chain contacts ({cutoff:g} nm, {scheme})", fontsize=10)
        return fig, ax
