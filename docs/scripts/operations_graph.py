#!/usr/bin/env python3
"""Render docs/assets/operations.png: what a SyncedStructure is made of, and which methods change or read each part.

The SyncedStructure cluster shows its parts, with the per-chain ChainRecord nested inside. Method groups on the left
point into the parts they change; groups on the right read parts and point to what they return. HaddockDocking is
one table-style box listing its own groups. Each part's subtitle says what it is built from. The PNG fits a 1920 x 1080 screen. The mkdocs hook in docs/hooks.py runs
this on every build, and tests/test_docs.py checks that every public name appears exactly once.

Run by hand with:  python docs/scripts/operations_graph.py
"""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from graphviz import Digraph

###################################### Configuration
FONT = "Arial"
MONO = "Courier"
BLACK = "#000000"
EDGE = "#89939D"
# Graphviz scales the drawing down to MAX_SIZE_IN at DPI: at most 1920 x 1080 px
DPI = 100
MAX_SIZE_IN = (19.2, 10.79)
OUTPUT = Path(__file__).resolve().parents[1] / "assets" / "operations"
WRAP_CHARS = 44                 # longest line of method names, in characters
TITLE_PT, PART_PT, METHOD_PT, CLUSTER_PT = 14, 12, 11, 18

# domain: (border, fill of method boxes, fill of data boxes)
DOMAINS = {
    "input": ("#397C91", "#E4EFF4", "#C9DFE8"),
    "write": ("#C46A43", "#F8EEE8", "#F2DED3"),
    "read": ("#8064A2", "#EFE9F4", "#DED2E9"),
    "docking": ("#648568", "#E9F2E7", "#D3E5D0"),
}
CLUSTER_FILL, CLUSTER_BORDER = "#EEF1F8", "#65749A"
RECORD_FILL = "#E1E6F2"
PART_FILL = "#FFFFFF"

###################################### Composition

# parts of a SyncedStructure: id → (attribute names, description)
PARTS = {
    "identity": ("structure_id · file_format · key_by", "identity"),
    "mmcif": ("mmcif_dict", "mmCIF metadata"),
    "ids": ("label ↔ author IDs", "from mmcif_dict atom_site"),
    "entities": ("entity sequences · accessions", "from mmcif_dict entity_poly, struct_ref"),
    "config": ("config · aligners", "YAML settings"),
    "structure": ("structure → model", "Biopython atoms"),
}
# fields of each ChainRecord in chains[key]
# fields of each ChainRecord in chains[key]; subtitles say what each is built from
RECORD = {
    "sequence": ("sequence", "observed, from the chain's residues"),
    "reference": ("reference_sequence", "from entity sequences"),
    "mapping": ("sequence_index_to_ref_index", "sequence aligned to reference_sequence"),
    "notes": ("entity_id · external_ref · refseq_id · metadata", "from entities, or set by you"),
    "chain_ids": ("chain_id · auth_chain_id", "from label ↔ author IDs"),
    "chain": ("chain · residues · hetero · moieties", "the Biopython chain in model"),
}

# inputs and outputs: id → (label, domain)
DATA = {
    "pdb_id": ("PDB ID", "input"),
    "file": ("PDB / mmCIF file", "input"),
    "biopython": ("Biopython Structure", "input"),
    "values": ("Sequences, IDs & maps", "read"),
    "tables": ("Tables (DataFrame)", "read"),
    "arrays": ("Strings & matrices", "read"),
    "new": ("New SyncedStructure", "read"),
    "figures": ("Figures & 3D views", "read"),
    "files": ("Files: TSV PDF PNG HTML CIF PDB", "read"),
}

# method groups: id → (title, domain, inputs, parts read, parts changed, outputs, methods)
OPERATIONS = {
    "rcsb": ("RCSB", "input", [], [], [], ["pdb_id"], ["search_rcsb", "fetch_rcsb_citations"]),
    "load": ("Load", "input", ["pdb_id", "file", "biopython"], [], ["identity", "structure", "mmcif", "config"], [],
             ["from_rcsb", "from_path", "from_biopython"]),
    "remap": ("Reference & mapping", "write", [], [], ["reference", "mapping"], [], [
        "set_reference_sequence", "set_reference_sequences", "remap_chain", "remap_all_chains"]),
    "annotate": ("Annotate", "write", [], [], ["notes"], [], [
        "set_entity_id", "set_accession", "set_reference_id", "set_chain_metadata"]),
    "chainset": ("Chain set", "write", [], [], ["structure", "chain_ids"], [], [
        "add_chain", "remove_chain", "replace_chain", "rename_chains", "rename_chains_by_sequence"]),
    "renumber": ("Numbering", "write", [], ["mapping"], ["chain"], [],
                 ["renumber_chain", "apply_reference_numbering"]),
    "lookup": ("Look up", "read", [], ["sequence", "reference", "notes"], [], ["values"], [
        "chain_ids", "get_sequence", "get_reference_sequence", "get_accession", "get_reference_id",
        "find_chains_by_sequence", "chain_ids_by_moiety"]),
    "tabulate": ("Chains, residues & alignment", "read", [], ["mapping", "notes"], [],
                 ["tables"], [
        "chain_table", "chain_residue_table", "residue_table", "moiety_table", "alignment_table",
        "alignment_column_table"]),
    "contacts": ("Contacts (nm)", "read", [], ["chain", "mapping"], [], ["tables", "arrays", "values"], [
        "contact_table", "contact_stat_tables", "contact_string", "reference_contact_string", "masked_sequence",
        "reference_masked_sequence", "contact_matrix", "block_contact_matrix", "contact_point_trajectory",
        "chain_partner_map"]),
    "interfaces": ("Interfaces (Å)", "read", [], ["chain", "mapping"], [], ["tables", "values", "new"], [
        "interface_metrics", "interface_tables", "interface_residue_table", "interface_partner_map",
        "chain_partner_metrics", "chain_partner_table", "interface_complex", "chain_partner_structures"]),
    "subset": ("Subset", "read", [], ["chain"], [], ["new"], ["subset_chains"]),
    "show": ("Plot & 3D", "read", [], ["chain"], [], ["figures"], [
        "plot_contact_map", "plot_multi_contact_map", "build_panel_grid", "show_structure_3d",
        "show_side_by_side_3d"]),
    "export": ("Export", "read", ["figures"], ["structure", "notes"], [], ["files"], [
        "export_structure", "mmcif_string", "pdb_string", "export_chain_table", "export_chain_residue_table",
        "export_residue_table", "export_moiety_table", "export_alignment_table", "export_alignment_column_table",
        "export_contact_strings", "export_contact_stat_tables", "export_interface_tables", "export_plot",
        "export_plot_grid", "export_structure_3d", "export_panels_pdf"]),
    "prepare": ("Prepare docking", "docking", [], ["chain"], [], ["docking"],
                ["prepare_haddock_docking", "HaddockDocking", "core_count"]),
}

# HaddockDocking box: (row title, methods), top to bottom
DOCKING = [
    ("Restrain & configure (in place)", [
        "add_ambiguous_restraints", "add_unambiguous_restraint", "add_patch_pair_restraints",
        "add_contact_restraints", "add_scan_residues", "clear_restraints", "set_reference"]),
    ("Inspect set up", [
        "HaddockDocking.residue_table", "restraint_table", "contact_patch_table", "partner_contact_map",
        "restraint_text", "workflow_modules", "config_toml", "input_paths", "run_dir"]),
    ("Write & run haddock3", ["export_inputs", "validate_config", "run_haddock"]),
    ("Results → tables, figures, SyncedStructures", [
        "step_paths", "score_table", "caprieval_tables", "alascan_tables", "docked_model_table",
        "docked_model_paths", "run_caprieval", "plot_caprieval", "load_docked_model", "load_docked_models"]),
]

# columns left to right (rankdir LR), each ordered top to bottom; parts are placed by the cluster
ROWS = [
    ["rcsb", "pdb_id", "file", "biopython"],
    ["load", "remap", "annotate", "chainset", "renumber"],
    ["lookup", "tabulate", "contacts", "interfaces", "subset", "show", "export", "prepare"],
    ["values", "tables", "arrays", "new", "figures", "files", "docking"],
]


def documented_names() -> list[str]:
    """Return every method name in the graph, as listed.

    Returns:
        Names in graph order; 'HaddockDocking.x' marks a HaddockDocking method that shares a SyncedStructure name.
    """
    return ([name for *_, methods in OPERATIONS.values() for name in methods]
            + [name for _, methods in DOCKING for name in methods])

###################################### Graph construction


def _methods(names: list[str], wrap: int = WRAP_CHARS) -> str:
    """Return method names as small monospace lines, wrapped at a width.

    Args:
        names: Method names; 'Class.method' prints as 'method' and a class name as 'Class(…)'.
        wrap: Longest line in characters; a longer single name gets its own line.

    Returns:
        HTML-like label text.
    """
    shown = [f"{name}(…)" if name[0].isupper() and "." not in name else name.split(".")[-1] for name in names]
    lines: list[str] = []
    for name in shown:
        if lines and len(lines[-1]) + 3 + len(name) <= wrap:
            lines[-1] += f" · {name}"
        else:
            lines.append(name)
    return f'<FONT FACE="{MONO}" POINT-SIZE="{METHOD_PT}">{"<BR/>".join(escape(line) for line in lines)}</FONT>'


def _part(graph: Digraph, node: str, attributes: str, note: str, fill: str) -> None:
    """Add a part of the structure: attribute names in monospace over a short description.

    Args:
        graph: Graph or cluster to add to.
        node: Node ID.
        attributes: Attribute names.
        note: Description, or empty.
        fill: Fill color.
    """
    text = f'<FONT FACE="{MONO}" POINT-SIZE="{PART_PT}"><B>{escape(attributes)}</B></FONT>'
    if note:
        text += f'<BR/><FONT POINT-SIZE="{METHOD_PT}" COLOR="#555555">{escape(note)}</FONT>'
    graph.node(node, label=f"<{text}>", fillcolor=fill, color=CLUSTER_BORDER, penwidth="1.2", shape="box",
               style="filled,rounded")


def _docking_label() -> str:
    """Return the HaddockDocking box as an HTML table: its name, then one row per method group.

    Returns:
        The Graphviz label.
    """
    border, fill, head = DOMAINS["docking"]
    rows = [f'<TR><TD BGCOLOR="{head}"><FONT POINT-SIZE="{TITLE_PT + 2}"><B>HaddockDocking</B></FONT></TD></TR>']
    rows += [f'<TR><TD BGCOLOR="{fill}"><B>{escape(title)}</B><BR/>{_methods(methods, int(WRAP_CHARS * 1.4))}</TD></TR>'
             for title, methods in DOCKING]
    return (f'<<TABLE BORDER="2" CELLBORDER="0" CELLSPACING="0" CELLPADDING="5" COLOR="{border}" STYLE="rounded">'
            f'{"".join(rows)}</TABLE>>')


def build_graph() -> Digraph:
    """Return the composition and operations graph.

    Returns:
        A graphviz.Digraph.
    """
    from graphviz import Digraph

    graph = Digraph("synced_structure", engine="dot", format="png")
    graph.attr(bgcolor="#ffffff", rankdir="LR", ranksep="0.5", nodesep="0.14", splines="spline", pad="0.15",
               dpi=str(DPI), size=f"{MAX_SIZE_IN[0]},{MAX_SIZE_IN[1]}", outputorder="edgesfirst", fontname=FONT)
    graph.attr("node", fontname=FONT, style="filled,rounded", shape="box", fontcolor=BLACK, margin="0.1,0.05")
    graph.attr("edge", color=EDGE, penwidth="1.1", arrowsize="0.6", arrowhead="vee")
    with graph.subgraph(name="cluster_structure") as cluster:
        cluster.attr(label="<<B>SyncedStructure</B>>", fontsize=str(CLUSTER_PT), style="filled,rounded",
                     fillcolor=CLUSTER_FILL, color=CLUSTER_BORDER, penwidth="2.5", margin="12")
        with cluster.subgraph(name="parts") as column:          # one column of structure-level parts
            column.attr(rank="same")
            for node, (attributes, note) in PARTS.items():
                _part(column, node, attributes, note, PART_FILL)
        with cluster.subgraph(name="cluster_record") as record:
            record.attr(label="<<B>chains[key] → ChainRecord</B>>", fontsize=str(TITLE_PT), style="filled,rounded",
                        fillcolor=RECORD_FILL, color=CLUSTER_BORDER, penwidth="1.5", margin="10")
            with record.subgraph(name="fields") as column:      # one column of ChainRecord fields
                column.attr(rank="same")
                for node, (attributes, note) in RECORD.items():
                    _part(column, node, attributes, note, PART_FILL)
    for node, (label, domain) in DATA.items():
        border, _, fill = DOMAINS[domain]
        graph.node(node, label=f"<<B>{escape(label)}</B>>", fontsize=str(TITLE_PT), fillcolor=fill, color=border,
                   penwidth="2")
    graph.node("docking", label=_docking_label(), shape="plain", style="")
    for node, (title, domain, inputs, reads, writes, outputs, methods) in OPERATIONS.items():
        border, fill, _ = DOMAINS[domain]
        graph.node(node, label=f"<<B>{escape(title)}</B><BR/>{_methods(methods)}>", fontsize=str(TITLE_PT),
                   fillcolor=fill, color=border, penwidth="1.4")
        for source in inputs:
            graph.edge(source, node, constraint="false" if source == "figures" else "true")
        for part in writes:
            graph.edge(node, part, color=border, penwidth="1.4")
        for part in reads:
            graph.edge(part, node)
        for target in outputs:
            graph.edge(node, target, color=border)
    for index, row in enumerate(ROWS):
        with graph.subgraph(name=f"row_{index}") as rank:
            rank.attr(rank="same")
            for node in row:
                rank.node(node)
        for left, right in zip(row[:-1], row[1:]):
            graph.edge(left, right, style="invis", weight="10")
    graph.attr(label=f'<<FONT POINT-SIZE="{METHOD_PT + 1}" COLOR="#555555">Left: methods that change the parts they '
                     f'point to. Right: methods that read the parts pointing at them.</FONT>>', labelloc="b")
    return graph


def render(output: Path = OUTPUT) -> Path:
    """Render the graph to a PNG, replacing the file only when its bytes change.

    Args:
        output: Path without suffix.

    Returns:
        The PNG path.
    """
    png = output.with_suffix(".png")
    data = build_graph().pipe(format="png")
    if not png.is_file() or png.read_bytes() != data:
        png.parent.mkdir(parents=True, exist_ok=True)
        png.write_bytes(data)
    return png


if __name__ == "__main__":
    print(render())
