"""Jupyter 3D views of synchronized structures using 3Dmol.js.

Views are clickable, with a centered title above them, per-chain colors, residue highlighting with the other residues drawn at a reduced alpha, representation styles, backgrounds, sizes, rotation angles with reset and current readout buttons, and side-by-side layouts.
"""

from __future__ import annotations

###################################### IMPORTS AND CONSTANTS

import hashlib
import json
import math
import time
import uuid
from collections.abc import Sequence
from html import escape

from IPython.display import HTML

from .base import StructureState
from .constants import PALETTE
from .utils import log_completion, log_info

VIZ_SIZE_PX = 500
VIZ_PANEL_SIZE_PX = 350
VIZ_BACKGROUND = "white"
VIZ_STYLES = ("cartoon", "stick", "sphere", "line")
VIZ_FORMAT = "cif"
VIZ_DIM_ALPHA = 0.25
VIZ_JS_URL = "https://3dmol.org/build/3Dmol-min.js"
DEFAULT_SHOW_ANGLES = True
TITLE_FONT_PX = 16
ANGLES_FONT_PX = 12
SIDE_BY_SIDE_GAP_PX = 16
CHAIN_PALETTE = tuple(PALETTE)
CHAIN_LABEL_JS = (
    "function(atom,viewer,event,container){"
    "if(!atom.label){"
    "atom.label=viewer.addLabel('chain '+atom.chain+' '+atom.resn+' '+atom.resi,"
    "{position:atom,backgroundColor:'black',fontColor:'white',fontSize:12});"
    "}else{viewer.removeLabel(atom.label);delete atom.label;}"
    "viewer.render();}"
)
ANGLES_JS = (
    "(function(q){"
    "var n=Math.hypot(q[0],q[1],q[2],q[3])||1;"
    "var x=q[0]/n,y=q[1]/n,z=q[2]/n,w=q[3]/n;"
    "var r=Math.atan2(2*(w*x+y*z),1-2*(x*x+y*y));"
    "var p=Math.asin(Math.max(-1,Math.min(1,2*(w*y-z*x))));"
    "var h=Math.atan2(2*(w*z+x*y),1-2*(y*y+z*z));"
    "return [r,p,h].map(function(t){return Math.round(t*180/Math.PI)});})"
)
VIZ_PAGE_TEMPLATE = """<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title></head>
<body style="margin:16px;background:white">__BODY__</body></html>
"""
VIZ_TEMPLATE = """<div id="__ELEMENT__" class="viz3d-view" style="position:relative;width:__SIZE__px;height:__SIZE__px"><div id="__ELEMENT__-view" style="width:100%;height:100%"></div>__OVERLAY__</div>
<script>
(function(){
var id="__ELEMENT__",text=__TEXT__,pal=__PALETTE__,pose=__POSE__,hl=__HIGHLIGHT__,dim=__DIM__,sticks=__STICKS__;
if(!window.$3Dmol&&!document.getElementById("3dmol-js")){var s=document.createElement("script");s.id="3dmol-js";s.src="__JS_URL__";document.head.appendChild(s);}
function box(){return document.getElementById(id+"-view");}
function out(){return document.getElementById(id+"-angles");}
function build(){
if(!window.$3Dmol){return setTimeout(build,100);}
var b=box();
b.innerHTML="";
var v=$3Dmol.createViewer(b,{});
b.__viewer=v;
v.addModel(__CIF__,"__FORMAT__");
for(var c in pal){
if(!hl){v.setStyle({chain:c},{__STYLE__:{color:pal[c]}});continue;}
v.setStyle({chain:c},{__STYLE__:{color:pal[c],opacity:dim}});
var r=hl[c];
if(!r||!r.length){continue;}
v.setStyle({chain:c,resi:r},{__STYLE__:{color:pal[c]}});
if(sticks){
var sel=v.selectedAtoms({chain:c,resi:r}),prot={},ser=[];
sel.forEach(function(a){if(a.atom==="CA"){prot[a.resi]=1;}});
sel.forEach(function(a){if(!(prot[a.resi]&&["N","C","O","OXT"].indexOf(a.atom)>=0)){ser.push(a.serial);}});
if(ser.length){v.addStyle({serial:ser},{stick:{color:pal[c]}});}
}
}
v.setBackgroundColor(__BACKGROUND__);
v.setClickable({},true,__LABEL__);
v.zoomTo();
var view=Array.prototype.slice.call(v.getView());
for(var i=0;i<4;i++){view[4+i]=pose[i];}
v.setView(view);
v.render();
}
function readCurrent(){
var v=box().__viewer;
if(!v){return "loading, click again";}
var raw=Array.prototype.slice.call(v.getView());
console.log("[JS-INFO] getView raw",raw);
var a=__ANGLES__(raw.slice(4,8).map(Number));
if(a.some(function(t){return !isFinite(t)})){return "unreadable, see console";}
return "x "+a[0]+"\u00b0 y "+a[1]+"\u00b0 z "+a[2]+"\u00b0";
}
build();
if(__SHOW__){
document.getElementById(id+"-reset").onclick=function(){
try{build();out().textContent=text;}
catch(e){console.error("[JS-INFO] reset failed",e);}
};
document.getElementById(id+"-current").onclick=function(){
try{out().textContent=readCurrent();}
catch(e){console.error("[JS-INFO] angle read failed",e);out().textContent="unreadable, see console";}
};
}
})();
</script>"""

###################################### ROTATION


Quaternion = tuple[float, float, float, float]


def _quaternion(angles: Sequence[float]) -> Quaternion:
    """Return the (x, y, z, w) quaternion of x, y, z rotations applied in x, y, z order.

    Args:
        angles: Rotations about x, y and z in degrees.

    Returns:
        The quaternion as (x, y, z, w).
    """
    def axis(index: int, degrees: float) -> Quaternion:
        """Return the quaternion of one rotation about a coordinate axis.

        Args:
            index: Axis index, 0 for x, 1 for y, 2 for z.
            degrees: Rotation in degrees.

        Returns:
            The quaternion as (x, y, z, w).
        """
        half = math.radians(degrees) / 2
        x, y, z = (math.sin(half) if j == index else 0.0 for j in range(3))
        return x, y, z, math.cos(half)

    def mul(p: Quaternion, q: Quaternion) -> Quaternion:
        """Return the Hamilton product p * q.

        Args:
            p: Left quaternion.
            q: Right quaternion.

        Returns:
            The product quaternion.
        """
        x1, y1, z1, w1 = p
        x2, y2, z2, w2 = q
        return (w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2, w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2, w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2)

    return mul(mul(axis(2, angles[2]), axis(1, angles[1])), axis(0, angles[0]))

###################################### HIGHLIGHT


def _residue_numbers(spec: Sequence[int] | str) -> list[int]:
    """Return sorted residue numbers from a list or a text such as "10-25,40".

    Args:
        spec: Residue numbers, or text of numbers and inclusive ranges separated by commas.

    Returns:
        The sorted unique residue numbers.

    Raises:
        ValueError: If the text holds something other than non-negative numbers and ranges.
    """
    if not isinstance(spec, str):
        return sorted({int(number) for number in spec})
    numbers: set[int] = set()
    for part in spec.replace(" ", "").split(","):
        first, _, last = part.partition("-")
        if not first.isdigit() or (last and not last.isdigit()):
            raise ValueError(f"Invalid residue selection {part!r}; use numbers and ranges such as '10-25,40'.")
        numbers.update(range(int(first), int(last or first) + 1))
    return sorted(numbers)

###################################### COLOR HELPERS


class Viz3DMixin(StructureState):
    """Jupyter 3Dmol.js views of the structure."""

    def _chain_colors(self, colors: dict[str, str] | None = None,
                      chain_palette: Sequence[str] = CHAIN_PALETTE) -> dict[str, str]:
        """Return a color per loaded chain, ignoring unloaded chains and filling gaps from the palette.

        Args:
            colors: Chain ID to color, or None.
            chain_palette: Colors for chains without an explicit color, in order.

        Returns:
            Chain key to color, for every loaded chain.
        """
        colors = {str(key): value for key, value in (colors or {}).items()}
        ignored = sorted(set(colors) - set(self.chain_ids))
        if ignored:
            log_info(f"{self.structure_id}: colors ignored for unloaded chains {ignored}")
        return {cid: colors.get(cid, chain_palette[index % len(chain_palette)])
                for index, cid in enumerate(self.chain_ids)}

    def _highlight_numbers(self, highlight: dict[str, Sequence[int] | str] | None) -> dict[str, list[int]] | None:
        """Return the residues to highlight per loaded chain key, ignoring unloaded chains.

        Args:
            highlight: Chain ID (key, label or author) to residue numbers or text such as "10-25,40", or None.

        Returns:
            Chain key to sorted PDB residue numbers, or None when nothing is to be highlighted.
        """
        if not highlight:
            return None
        selected: dict[str, list[int]] = {}
        ignored = []
        for chain_id, spec in highlight.items():
            try:
                key = self._chain_key(chain_id)
            except KeyError:
                ignored.append(str(chain_id))
                continue
            selected[key] = sorted(set(selected.get(key, [])) | set(_residue_numbers(spec)))
        if ignored:
            log_info(f"{self.structure_id}: highlight ignored for unloaded chains {sorted(ignored)}")
        return selected

###################################### OVERLAY

    def _angles_overlay_html(self, element: str, text: str, angles_font_px: int = ANGLES_FONT_PX) -> str:
        """Return the bottom right overlay with the angle text and the reset and current buttons.

        Args:
            element: DOM id prefix of the view.
            text: Angle text to show initially.
            angles_font_px: Font size of the overlay in pixels.

        Returns:
            The overlay HTML.
        """
        button_style = f'cursor:pointer;background:none;border:none;color:white;font-size:{angles_font_px}px;padding:0'
        return (f'<div style="position:absolute;right:6px;bottom:6px;z-index:1;display:flex;align-items:center;gap:6px;'
                f'background:rgba(128,128,128,0.75);color:white;font-family:monospace;font-size:{angles_font_px}px;'
                f'padding:2px 6px;border-radius:3px"><span id="{element}-angles">{escape(text)}</span>'
                f'<button id="{element}-reset" title="Reset to input angles" style="{button_style}">&#8635;</button>'
                f'<button id="{element}-current" title="Show current angles" style="{button_style}">&#9673;</button></div>')

###################################### RENDERER

    def _viewer_html(self, angles: Sequence[float], chain_id: str | None, palette: dict[str, str], style: str,
                     background: str, size: int, show_angles: bool, viz_format: str = VIZ_FORMAT,
                     js_url: str = VIZ_JS_URL, label_js: str = CHAIN_LABEL_JS, angles_js: str = ANGLES_JS,
                     angles_font_px: int = ANGLES_FONT_PX, highlight: dict[str, list[int]] | None = None,
                     dim_alpha: float = VIZ_DIM_ALPHA, highlight_sticks: bool = False) -> str:
        """Return HTML for a 3Dmol viewer posed at the input angles, with reset and current angle buttons.

        Args:
            angles: Rotation (x, y, z) in degrees.
            chain_id: Chain to show, or None for all chains.
            palette: Chain ID to color.
            style: Representation, one of VIZ_STYLES.
            background: Background color.
            size: Width and height of the view in pixels.
            show_angles: Whether to show the angles and buttons in the bottom right corner.
            viz_format: Structure format passed to 3Dmol.js.
            js_url: URL of the 3Dmol.js script.
            label_js: JavaScript click handler that toggles atom labels.
            angles_js: JavaScript that converts a quaternion into x, y, z degrees.
            angles_font_px: Font size of the angle overlay in pixels.
            highlight: Chain key to residue numbers drawn at full opacity, or None to draw everything normally.
            dim_alpha: Opacity of the residues that are not highlighted.
            highlight_sticks: Whether highlighted residues also show their side chains as sticks.

        Returns:
            The viewer HTML with its script.
        """
        key = (self.structure_id, tuple(angles), chain_id, sorted(palette.items()), style, background, size, show_angles,
               viz_format, js_url, label_js, angles_js, angles_font_px, sorted((highlight or {}).items()), dim_alpha,
               highlight_sticks)
        element = "viz_" + hashlib.sha1(repr(key).encode()).hexdigest()[:12]
        text = f"x {angles[0]}\u00b0 y {angles[1]}\u00b0 z {angles[2]}\u00b0"
        overlay = self._angles_overlay_html(element, text, angles_font_px) if show_angles else ""
        replacements = {
            "__ELEMENT__": element, "__SIZE__": str(size), "__OVERLAY__": overlay, "__JS_URL__": js_url,
            "__FORMAT__": viz_format, "__PALETTE__": json.dumps(palette), "__STYLE__": style,
            "__BACKGROUND__": json.dumps(background), "__LABEL__": label_js,
            "__POSE__": json.dumps(list(_quaternion(angles))), "__SHOW__": "true" if show_angles else "false",
            "__ANGLES__": angles_js, "__TEXT__": json.dumps(text),
            "__HIGHLIGHT__": json.dumps(highlight), "__DIM__": json.dumps(dim_alpha),
            "__STICKS__": "true" if highlight_sticks else "false",
            "__CIF__": json.dumps(self._mmcif_text(None if chain_id is None else [chain_id], key_as_author=True)),
        }
        html = VIZ_TEMPLATE
        for key_name, value in replacements.items():
            html = html.replace(key_name, value)
        return html

###################################### PUBLIC ENTRY POINT

    def show_structure_3d(self, angles: Sequence[float] = (0, 0, 0), chain_id: str | None = None,
                          colors: dict[str, str] | None = None, title: str | None = None, style: str = "cartoon",
                          background: str = VIZ_BACKGROUND, size: int = VIZ_SIZE_PX,
                          show_angles: bool = DEFAULT_SHOW_ANGLES, chain_palette: Sequence[str] = CHAIN_PALETTE,
                          viz_format: str = VIZ_FORMAT, js_url: str = VIZ_JS_URL, label_js: str = CHAIN_LABEL_JS,
                          angles_js: str = ANGLES_JS, title_font_px: int = TITLE_FONT_PX,
                          angles_font_px: int = ANGLES_FONT_PX, highlight: dict[str, Sequence[int] | str] | None = None,
                          dim_alpha: float = VIZ_DIM_ALPHA, highlight_sticks: bool = False) -> HTML:
        """Show the structure in Jupyter as a clickable 3Dmol.js view, with an optional centered title.

        Args:
            angles: Rotation (x, y, z) in degrees, applied as the pose and shown when show_angles is True.
            chain_id: Single chain to show, or None for all chains.
            colors: Chain ID to color, e.g. {'A': 'red'}; unloaded chains are ignored.
            title: Title shown centered above the view.
            style: Representation: 'cartoon', 'stick', 'sphere' or 'line'.
            background: Background color, e.g. 'white' or '#f0f0f0'.
            size: Width and height of the view in pixels.
            show_angles: Whether to show the angles with reset and read-current buttons.
            chain_palette: Colors for chains without an explicit color, in order.
            viz_format: Structure format passed to 3Dmol.js.
            js_url: URL of the 3Dmol.js script.
            label_js: JavaScript click handler that toggles atom labels.
            angles_js: JavaScript that converts a quaternion into x, y, z degrees.
            title_font_px: Font size of the title in pixels.
            angles_font_px: Font size of the angle overlay in pixels.
            highlight: Chain ID to the residues to highlight, as PDB residue numbers or text such as "10-25,40",
                e.g. {'A': '10-25', 'B': [4, 7]}. Every other residue, also of chains not listed, is drawn at dim_alpha.
                Unloaded chains are ignored. None draws everything normally.
            dim_alpha: Opacity from 0 to 1 of the residues that are not highlighted; only used with highlight.
            highlight_sticks: Whether highlighted residues also show their side chains as sticks.

        Returns:
            An IPython HTML object with the title and the view.

        Raises:
            ValueError: If angles does not hold three values, style is unknown, dim_alpha is outside 0 to 1, or a
                highlight text is malformed.
        """
        if len(angles) != 3:
            raise ValueError("angles must contain three values (x, y, z) in degrees.")
        if style not in VIZ_STYLES:
            raise ValueError(f"style must be one of {VIZ_STYLES}.")
        if not 0 <= dim_alpha <= 1:
            raise ValueError("dim_alpha must be between 0 and 1.")
        start = time.perf_counter()
        palette = self._chain_colors(colors, chain_palette)
        html = self._viewer_html(angles, chain_id, palette, style, background, size, show_angles,
                                 viz_format, js_url, label_js, angles_js, angles_font_px,
                                 self._highlight_numbers(highlight), dim_alpha, highlight_sticks)
        if title:
            heading = (f'<div class="viz3d-title" style="width:{size}px;margin:0 auto 4px auto;text-align:center;font-family:sans-serif;'
                       f'font-size:{title_font_px}px;font-weight:bold">{escape(str(title))}</div>')
            html = heading + html
        log_completion(start, f"Rendered view of {self.structure_id} at angles {tuple(angles)}")
        return HTML(html)


###################################### LAYOUT


def show_side_by_side_3d(*outputs: HTML, gap_px: int = SIDE_BY_SIDE_GAP_PX, size_px: int = VIZ_PANEL_SIZE_PX) -> HTML:
    """Place several show_structure_3d outputs in a row, wrapping when the notebook is narrow.

    Every view becomes a size_px square panel, whatever size it was rendered at.

    Args:
        *outputs: Views returned by show_structure_3d, left to right.
        gap_px: Spacing between views in pixels.
        size_px: Width and height of each panel in pixels.

    Returns:
        One HTML object holding all the views.

    Raises:
        ValueError: If no view is given.
    """
    if not outputs:
        raise ValueError("show_side_by_side_3d needs at least one view.")
    grid_id = f"viz3d-panels-{uuid.uuid4().hex[:8]}"
    style = (f"<style>#{grid_id} .viz3d-view{{width:{size_px}px !important;height:{size_px}px !important}}"
             f"#{grid_id} .viz3d-title{{width:{size_px}px !important}}</style>")
    cells = "".join(f'<div style="flex:0 0 auto">{output.data}</div>' for output in outputs)
    return HTML(f'{style}<div id="{grid_id}" style="display:flex;flex-wrap:wrap;align-items:flex-start;'
                f'gap:{gap_px}px">{cells}</div>')
