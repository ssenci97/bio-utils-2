#!/usr/bin/env python
"""Jupyter visualization for synchronized structures: clickable py3Dmol views using mmCIF input, with a title inside the view, per-chain colors, representation styles, backgrounds, sizes, rotation angles, and side-by-side layouts."""

###################################### IMPORTS AND CONSTANTS

import copy
import io
import time
from html import escape

import py3Dmol
from Bio.PDB import MMCIFIO
from IPython.display import HTML

from helpers import finish_log

SEED = 42
VIZ_SIZE_PX = 400
VIZ_BACKGROUND = "white"
VIZ_STYLES = ("cartoon", "stick", "sphere", "line")
VIZ_FORMAT = "cif"
TITLE_FONT_PX = 16
SIDE_BY_SIDE_GAP_PX = 16
CHAIN_PALETTE = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                 "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
CHAIN_LABEL_JS = (
    "function(atom,viewer,event,container){"
    "if(!atom.label){"
    "atom.label=viewer.addLabel('chain '+atom.chain+' '+atom.resn+' '+atom.resi,"
    "{position:atom,backgroundColor:'black',fontColor:'white',fontSize:12});"
    "}else{viewer.removeLabel(atom.label);delete atom.label;}"
    "viewer.render();}"
)

###################################### COLOR HELPERS


class Viz3DMixin:
    def _chain_colors(self, colors=None):
        """Return {chain_id: color} for loaded chains, ignoring unloaded chains and filling gaps from the palette.

        Args:
            colors: Dict mapping chain IDs to colors, or None.
        """
        colors = {str(key): value for key, value in (colors or {}).items()}
        ignored = sorted(set(colors) - set(self.chain_ids))
        if ignored:
            print(f"[PYTHON-INFO] {self.structure_id}: colors ignored for unloaded chains {ignored}")
        return {cid: colors.get(cid, CHAIN_PALETTE[index % len(CHAIN_PALETTE)])
                for index, cid in enumerate(self.chain_ids)}

###################################### STRUCTURE TEXT

    def _mmcif_string(self, chain_id=None):
        """Return the structure, or one chain, as mmCIF text, which has no atom serial limit.

        Args:
            chain_id: Chain to export, or None for all chains.
        """
        selected = {self._require_chain(chain_id)} if chain_id is not None else set(self.chains)
        structure = copy.deepcopy(self.structure)
        model = next(structure.get_models())
        for chain in list(model.get_chains()):
            if str(chain.id) not in selected:
                model.detach_child(chain.id)
        writer = MMCIFIO()
        writer.set_structure(structure)
        output = io.StringIO()
        writer.save(output)
        return output.getvalue()

###################################### RENDERER

    def _interactive_view(self, angles, chain_id, palette, style, background, size):
        """Return a clickable py3Dmol view with per-chain colors, rotated by the given angles.

        Args:
            angles: Rotation (x, y, z) in degrees.
            chain_id: Chain to show, or None for all chains.
            palette: Dict mapping chain IDs to colors.
            style: Representation, one of VIZ_STYLES.
            background: Background color.
            size: Width and height of the view in pixels.
        """
        view = py3Dmol.view(width=size, height=size)
        view.addModel(self._mmcif_string(chain_id=chain_id), VIZ_FORMAT)
        for cid, color in palette.items():
            view.setStyle({"chain": cid}, {style: {"color": color}})
        view.setBackgroundColor(background)
        view.setClickable({}, True, CHAIN_LABEL_JS)
        for angle, axis in zip(angles, "xyz"):
            if angle:
                view.rotate(angle, axis=axis)
        view.zoomTo()
        return view

###################################### PUBLIC ENTRY POINT

    def jupy_viz(self, angles=(0, 0, 0), chain_id=None, colors=None, title=None,
                 style="cartoon", background=VIZ_BACKGROUND, size=VIZ_SIZE_PX):
        """Show the structure in a Jupyter notebook as a clickable py3Dmol view with the title inside the view.

        Args:
            angles: Rotation (x, y, z) in degrees.
            chain_id: Optional single chain to show, or None for all chains.
            colors: Optional dict mapping chain IDs to colors, e.g. {'A': 'red', 'B': '#00aa00'}. Unloaded chains are ignored.
            title: Optional title text rendered above the view in the same output.
            style: Representation, one of 'cartoon', 'stick', 'sphere', or 'line'.
            background: Background color, e.g. 'white', 'black', or '#f0f0f0'.
            size: Width and height of the view in pixels.

        Returns:
            An IPython HTML object containing the title and the py3Dmol view.
        """
        if len(angles) != 3:
            raise ValueError("angles must contain three values (x, y, z) in degrees.")
        if style not in VIZ_STYLES:
            raise ValueError(f"style must be one of {VIZ_STYLES}.")
        start = time.perf_counter()
        view = self._interactive_view(angles, chain_id, self._chain_colors(colors), style, background, size)
        html = view._make_html()
        if title:
            heading = f'<div style="font-family:sans-serif;font-size:{TITLE_FONT_PX}px;font-weight:bold;margin-bottom:4px">{escape(str(title))}</div>'
            html = heading + html
        finish_log(start, f"Rendered interactive view of {self.structure_id} at angles {tuple(angles)}")
        return HTML(html)

###################################### LAYOUT


def side_by_side(*outputs, gap_px=SIDE_BY_SIDE_GAP_PX):
    """Return one HTML object that places several jupy_viz outputs in a row, wrapping if the notebook is narrow.

    Args:
        *outputs: HTML objects returned by jupy_viz, in left-to-right order.
        gap_px: Horizontal and vertical spacing between views in pixels.
    """
    if not outputs:
        raise ValueError("side_by_side needs at least one view.")
    cells = "".join(f'<div style="flex:0 0 auto">{output.data}</div>' for output in outputs)
    return HTML(f'<div style="display:flex;flex-wrap:wrap;align-items:flex-start;gap:{gap_px}px">{cells}</div>')