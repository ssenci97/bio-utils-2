#!/usr/bin/env python
"""Jupyter visualization for synchronized structures: clickable py3Dmol views using mmCIF input, with a centered title inside the view, per-chain colors, representation styles, backgrounds, sizes, rotation angles with reset and current readout buttons, and side-by-side layouts."""

###################################### IMPORTS AND CONSTANTS

import copy
import hashlib
import io
import json
import math
import time
from html import escape

from Bio.PDB import MMCIFIO
from IPython.display import HTML

from helpers import finish_log

VIZ_SIZE_PX = 400
VIZ_BACKGROUND = "white"
VIZ_STYLES = ("cartoon", "stick", "sphere", "line")
VIZ_FORMAT = "cif"
VIZ_JS_URL = "https://3dmol.org/build/3Dmol-min.js"
SHOW_ANGLES = True
TITLE_FONT_PX = 16
ANGLES_FONT_PX = 12
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
ANGLES_JS = (
    "(function(q){"
    "var n=Math.hypot(q[0],q[1],q[2],q[3])||1;"
    "var x=q[0]/n,y=q[1]/n,z=q[2]/n,w=q[3]/n;"
    "var r=Math.atan2(2*(w*x+y*z),1-2*(x*x+y*y));"
    "var p=Math.asin(Math.max(-1,Math.min(1,2*(w*y-z*x))));"
    "var h=Math.atan2(2*(w*z+x*y),1-2*(y*y+z*z));"
    "return [r,p,h].map(function(t){return Math.round(t*180/Math.PI)});})"
)
VIZ_TEMPLATE = """<div id="__ELEMENT__" style="position:relative;width:__SIZE__px;height:__SIZE__px"><div id="__ELEMENT__-view" style="width:100%;height:100%"></div>__OVERLAY__</div>
<script>
(function(){
var id="__ELEMENT__",text=__TEXT__,pal=__PALETTE__,pose=__POSE__;
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
for(var c in pal){v.setStyle({chain:c},{__STYLE__:{color:pal[c]}});}
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
console.log("[PYTHON-INFO] getView raw",raw);
var a=__ANGLES__(raw.slice(4,8).map(Number));
if(a.some(function(t){return !isFinite(t)})){return "unreadable, see console";}
return "x "+a[0]+"\u00b0 y "+a[1]+"\u00b0 z "+a[2]+"\u00b0";
}
build();
if(__SHOW__){
document.getElementById(id+"-reset").onclick=function(){
try{build();out().textContent=text;}
catch(e){console.error("[PYTHON-INFO] reset failed",e);}
};
document.getElementById(id+"-current").onclick=function(){
try{out().textContent=readCurrent();}
catch(e){console.error("[PYTHON-INFO] angle read failed",e);out().textContent="unreadable, see console";}
};
}
})();
</script>"""

###################################### ROTATION


def _quaternion(angles):
    def axis(index, degrees):
        half = math.radians(degrees) / 2
        return tuple(math.sin(half) if j == index else 0.0 for j in range(3)) + (math.cos(half),)

    def mul(p, q):
        x1, y1, z1, w1 = p
        x2, y2, z2, w2 = q
        return (w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2, w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2, w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2)

    return mul(mul(axis(2, angles[2]), axis(1, angles[1])), axis(0, angles[0]))

###################################### COLOR HELPERS


class Viz3DMixin:
    def _chain_colors(self, colors=None, chain_palette=CHAIN_PALETTE):
        """Return {chain_id: color} for loaded chains, ignoring unloaded chains and filling gaps from the palette.

        Args:
            colors: Dict mapping chain IDs to colors, or None.
            chain_palette: Colors used for chains without an explicit color, in order.
        """
        colors = {str(key): value for key, value in (colors or {}).items()}
        ignored = sorted(set(colors) - set(self.chain_ids))
        if ignored:
            print(f"[PYTHON-INFO] {self.structure_id}: colors ignored for unloaded chains {ignored}")
        return {cid: colors.get(cid, chain_palette[index % len(chain_palette)])
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

###################################### OVERLAY

    def _angles_overlay(self, element, text, angles_font_px=ANGLES_FONT_PX):
        """Return the bottom right overlay with the angle text and the reset and current buttons.

        Args:
            element: DOM id prefix of the view.
            text: Angle text to show initially.
            angles_font_px: Font size of the overlay in pixels.
        """
        button_style = f'cursor:pointer;background:none;border:none;color:white;font-size:{angles_font_px}px;padding:0'
        return (f'<div style="position:absolute;right:6px;bottom:6px;z-index:1;display:flex;align-items:center;gap:6px;'
                f'background:rgba(128,128,128,0.75);color:white;font-family:monospace;font-size:{angles_font_px}px;'
                f'padding:2px 6px;border-radius:3px"><span id="{element}-angles">{escape(text)}</span>'
                f'<button id="{element}-reset" title="Reset to input angles" style="{button_style}">&#8635;</button>'
                f'<button id="{element}-current" title="Show current angles" style="{button_style}">&#9673;</button></div>')

###################################### RENDERER

    def _viewer_html(self, angles, chain_id, palette, style, background, size, show_angles,
                     viz_format=VIZ_FORMAT, js_url=VIZ_JS_URL, label_js=CHAIN_LABEL_JS,
                     angles_js=ANGLES_JS, angles_font_px=ANGLES_FONT_PX):
        """Return HTML for a 3Dmol viewer posed at the input angles, with reset and current angle buttons.

        Args:
            angles: Rotation (x, y, z) in degrees.
            chain_id: Chain to show, or None for all chains.
            palette: Dict mapping chain IDs to colors.
            style: Representation, one of VIZ_STYLES.
            background: Background color.
            size: Width and height of the view in pixels.
            show_angles: If True, display the angles and the buttons in the bottom right corner.
            viz_format: Structure format passed to py3Dmol.
            js_url: URL of the py3Dmol script.
            label_js: JavaScript click handler that toggles atom labels.
            angles_js: JavaScript that converts a quaternion into x, y, z degrees.
            angles_font_px: Font size of the angle overlay in pixels.
        """
        key = (self.structure_id, tuple(angles), chain_id, sorted(palette.items()), style, background, size, show_angles)
        element = "viz_" + hashlib.sha1(repr(key).encode()).hexdigest()[:12]
        text = f"x {angles[0]}\u00b0 y {angles[1]}\u00b0 z {angles[2]}\u00b0"
        overlay = self._angles_overlay(element, text, angles_font_px) if show_angles else ""
        replacements = {
            "__ELEMENT__": element, "__SIZE__": str(size), "__OVERLAY__": overlay, "__JS_URL__": js_url,
            "__FORMAT__": viz_format, "__PALETTE__": json.dumps(palette), "__STYLE__": style,
            "__BACKGROUND__": json.dumps(background), "__LABEL__": label_js,
            "__POSE__": json.dumps(list(_quaternion(angles))), "__SHOW__": "true" if show_angles else "false",
            "__ANGLES__": angles_js, "__TEXT__": json.dumps(text),
            "__CIF__": json.dumps(self._mmcif_string(chain_id=chain_id)),
        }
        html = VIZ_TEMPLATE
        for key_name, value in replacements.items():
            html = html.replace(key_name, value)
        return html

###################################### PUBLIC ENTRY POINT

    def jupy_viz(self, angles=(0, 0, 0), chain_id=None, colors=None, title=None,
                 style="cartoon", background=VIZ_BACKGROUND, size=VIZ_SIZE_PX, show_angles=SHOW_ANGLES,
                 chain_palette=CHAIN_PALETTE, viz_format=VIZ_FORMAT, js_url=VIZ_JS_URL, label_js=CHAIN_LABEL_JS,
                 angles_js=ANGLES_JS, title_font_px=TITLE_FONT_PX, angles_font_px=ANGLES_FONT_PX):
        """Show the structure in a Jupyter notebook as a clickable py3Dmol view with the centered title inside the view.

        Args:
            angles: Rotation (x, y, z) in degrees, applied as the pose and displayed when show_angles is True.
            chain_id: Optional single chain to show, or None for all chains.
            colors: Optional dict mapping chain IDs to colors, e.g. {'A': 'red', 'B': '#00aa00'}. Unloaded chains are ignored.
            title: Optional title text rendered centered above the view in the same output.
            style: Representation, one of 'cartoon', 'stick', 'sphere', or 'line'.
            background: Background color, e.g. 'white', 'black', or '#f0f0f0'.
            size: Width and height of the view in pixels.
            show_angles: If True, display the angles with a reset button (rebuilds the input pose) and a current button (reads the live orientation).
            chain_palette: Colors used for chains without an explicit color, in order.
            viz_format: Structure format passed to py3Dmol.
            js_url: URL of the py3Dmol script.
            label_js: JavaScript click handler that toggles atom labels.
            angles_js: JavaScript that converts a quaternion into x, y, z degrees.
            title_font_px: Font size of the title in pixels.
            angles_font_px: Font size of the angle overlay in pixels.

        Returns:
            An IPython HTML object containing the title and the py3Dmol view.
        """
        if len(angles) != 3:
            raise ValueError("angles must contain three values (x, y, z) in degrees.")
        if style not in VIZ_STYLES:
            raise ValueError(f"style must be one of {VIZ_STYLES}.")
        start = time.perf_counter()
        palette = self._chain_colors(colors, chain_palette)
        html = self._viewer_html(angles, chain_id, palette, style, background, size, show_angles,
                                 viz_format, js_url, label_js, angles_js, angles_font_px)
        if title:
            heading = (f'<div style="width:{size}px;margin:0 auto 4px auto;text-align:center;font-family:sans-serif;'
                       f'font-size:{title_font_px}px;font-weight:bold">{escape(str(title))}</div>')
            html = heading + html
        finish_log(start, f"Rendered view of {self.structure_id} at angles {tuple(angles)}")
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