"""Build the single-page HTML API reference of synced_structure from its signatures and docstrings.

Usage:
    python docs/build_docs.py                  # writes docs/index.html
    python docs/build_docs.py -o site/api.html

Every function, method, property and class of the package is listed with its signature and its Google-style
docstring parsed into arguments, returns and raises. Internal (underscore) names are included and hidden by a
toggle. The page is self-contained apart from web fonts, and the output is deterministic, so it can be
committed and served by GitHub Pages from the docs folder.
"""

from __future__ import annotations

import argparse
import ast
import html
import importlib
import inspect
import re
import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))   # document the code in src/, never a stale installed copy

import synced_structure  # noqa: E402
from synced_structure.base import StructureState  # noqa: E402
from synced_structure.constants import DEFAULT_CONFIG  # noqa: E402

PACKAGE = "synced_structure"
REPO = "https://github.com/ssenci97/bio-utils-2"
SOURCE_URL = f"{REPO}/blob/main/src/{PACKAGE}"
OUTPUT = Path(__file__).resolve().parent / "index.html"

# Sidebar groups: (group title, [(module, section title)])
LAYOUT = [
    ("Structure", [("structure", "SyncedStructure"), ("loaders", "Loading"), ("chains", "Chains"),
                   ("metadata", "mmCIF metadata"), ("residues", "Residues and numbering"),
                   ("alignment", "Reference alignment")]),
    ("Analyses", [("contacts", "Contacts"), ("contact_matrix", "Contact matrices"),
                  ("contact_stats", "Contact statistics"), ("interface_metrics", "Interfaces")]),
    ("Outputs", [("contact_plots", "Contact maps"), ("panels", "Panel grids"), ("viz3d", "3D views"),
                 ("exporter", "Export")]),
    ("Docking", [("docking", "HADDOCK3 docking")]),
    ("Utilities", [("rcsb", "RCSB search"), ("sequence", "Sequences"), ("utils", "Configuration and logging"),
                   ("_types", "Type aliases")]),
]
GUIDE_SECTIONS = [("overview", "Overview"), ("conventions", "Conventions"), ("concepts", "Concepts"),
                  ("configuration", "Configuration")]
SECTION_HEADERS = {"Args", "Arguments", "Returns", "Yields", "Raises", "Attributes", "Example", "Examples",
                   "Note", "Notes"}
CONFIG_HELP = {
    "save_outdir": "Default output folder of every export_ method.",
    "aln_match_score": "Local alignment match score (reference alignment tables, sequence matching).",
    "aln_mismatch_score": "Local alignment mismatch score.",
    "aln_open_gap_score": "Local alignment gap opening score.",
    "aln_extend_gap_score": "Local alignment gap extension score.",
    "interface_contact_cutoff_a": "Heavy-atom contact cutoff of interface metrics, Å.",
    "interface_score_cutoff_a": "Distance range of the interface contact score, Å.",
    "interface_block_size": "Rows per block when computing interface atom distances (memory).",
}
KEY_ENTRY_POINTS = ["SyncedStructure.from_rcsb", "SyncedStructure.from_path", "SyncedStructure.chain_table",
                    "SyncedStructure.contact_table", "SyncedStructure.interface_metrics",
                    "SyncedStructure.plot_contact_map", "SyncedStructure.show_structure_3d",
                    "SyncedStructure.prepare_haddock_docking"]


###################################### MODEL

@dataclass
class Doc:
    """A docstring split into summary, description paragraphs and sections."""

    summary: str = ""
    paragraphs: list[str] = field(default_factory=list)
    args: list[tuple[str, str]] = field(default_factory=list)
    returns: str = ""
    raises: list[tuple[str, str]] = field(default_factory=list)
    attributes: list[tuple[str, str]] = field(default_factory=list)
    extra: list[tuple[str, str]] = field(default_factory=list)


@dataclass
class Entry:
    """One documented object: a function, method, property or class."""

    anchor: str
    title: str
    name: str
    kind: str
    obj: Any
    signature: inspect.Signature | None
    doc: Doc
    internal: bool
    source: str
    note: str = ""


@dataclass
class Section:
    """One module of the package, rendered as one section of the page."""

    anchor: str
    title: str
    module: str
    lede: str
    import_line: str
    entries: list[Entry]


###################################### DOCSTRINGS

def parse_docstring(text: str | None) -> Doc:
    """Split a Google-style docstring into its parts.

    Args:
        text: Cleaned docstring, as returned by inspect.getdoc, or None.

    Returns:
        The parsed docstring.
    """
    doc = Doc()
    if not text:
        return doc
    blocks: dict[str, list[str]] = {}
    body: list[str] = []
    current: list[str] = body
    for line in text.splitlines():
        header = line.strip().rstrip(":")
        if line.strip().endswith(":") and not line.startswith(" ") and header in SECTION_HEADERS:
            current = blocks.setdefault(header, [])
            continue
        current.append(line)
    paragraphs = [" ".join(part.split()) for part in "\n".join(body).split("\n\n") if part.strip()]
    if paragraphs:
        doc.summary, doc.paragraphs = paragraphs[0], paragraphs[1:]
    doc.args = _items(blocks.pop("Args", []) + blocks.pop("Arguments", []))
    doc.raises = _items(blocks.pop("Raises", []))
    doc.attributes = _items(blocks.pop("Attributes", []))
    doc.returns = " ".join(" ".join(blocks.pop("Returns", []) + blocks.pop("Yields", [])).split())
    doc.extra = [(name, "\n".join(lines).strip()) for name, lines in blocks.items()]
    return doc


def _items(lines: list[str]) -> list[tuple[str, str]]:
    """Return the (name, description) items of an Args-like section, joining continuation lines.

    Args:
        lines: Section lines, items indented once and continuations deeper.

    Returns:
        One pair per item.
    """
    items: list[tuple[str, str]] = []
    indent = min((len(line) - len(line.lstrip()) for line in lines if line.strip()), default=0)
    for line in lines:
        if not line.strip():
            continue
        if len(line) - len(line.lstrip()) == indent and ":" in line:
            name, _, text = line.strip().partition(":")
            items.append((name.strip(), text.strip()))
        elif items:
            items[-1] = (items[-1][0], f"{items[-1][1]} {line.strip()}")
    return items


###################################### COLLECTION

def _source(obj: Any) -> str:
    """Return 'module.py:line' of an object's definition, or '' when it cannot be found.

    Args:
        obj: Function, class or property getter.

    Returns:
        The location relative to the package folder.
    """
    try:
        target = inspect.unwrap(obj)
        path = Path(inspect.getsourcefile(target) or "")
        line = inspect.getsourcelines(target)[1]
    except (TypeError, OSError):
        return ""
    return f"{path.name}:{line}"


def _signature(obj: Any) -> inspect.Signature | None:
    """Return an object's signature, or None when it has none.

    Args:
        obj: Callable.

    Returns:
        The signature.
    """
    try:
        return inspect.signature(obj)
    except (TypeError, ValueError):
        return None


def _class_members(cls: type, owner: str) -> Iterator[Entry]:
    """Yield the methods and properties a class defines itself, public first.

    Args:
        cls: Class to read.
        owner: Name shown before each member, such as 'SyncedStructure'.

    Yields:
        One entry per member.
    """
    found = []
    for name, raw in cls.__dict__.items():
        if isinstance(raw, property):
            kind, func = "property", raw.fget
        elif isinstance(raw, classmethod):
            kind, func = "classmethod", raw.__func__
        elif isinstance(raw, staticmethod):
            kind, func = "staticmethod", raw.__func__
        elif inspect.isfunction(raw):
            kind, func = "method", raw
        else:
            continue
        if name == "__init__" or (name.startswith("__") and not inspect.getdoc(func)):
            continue
        internal = name.startswith("_") and not name.startswith("__")
        signature = _signature(func)
        if signature is not None and kind != "staticmethod":
            params = list(signature.parameters.values())[1:]   # drop self / cls
            signature = signature.replace(parameters=params)
        title = f"{owner}.{name}"
        found.append(Entry(anchor=title, title=title, name=name, kind=kind, obj=func, signature=signature,
                           doc=parse_docstring(inspect.getdoc(func)), internal=internal, source=_source(func)))
    yield from sorted(found, key=lambda entry: entry.internal)


def _class_entry(cls: type) -> Entry:
    """Return the entry of a class itself: its constructor signature and docstring.

    Args:
        cls: Class to read.

    Returns:
        The class entry.
    """
    signature = _signature(cls) if "__init__" in cls.__dict__ or hasattr(cls, "__dataclass_fields__") else None
    kind = "typeddict" if hasattr(cls, "__total__") else "class"
    if kind == "typeddict":
        fields = [inspect.Parameter(name, inspect.Parameter.KEYWORD_ONLY, annotation=annotation)
                  for name, annotation in cls.__annotations__.items()]
        signature = inspect.Signature(fields)
    return Entry(anchor=cls.__name__, title=cls.__name__, name=cls.__name__, kind=kind, obj=cls,
                 signature=signature, doc=parse_docstring(inspect.getdoc(cls)),
                 internal=cls.__name__.startswith("_"), source=_source(cls))


def collect_section(module_name: str, title: str) -> Section:
    """Collect every documented object defined in one module of the package.

    Mixin members are shown as members of SyncedStructure, since that is where users call them.

    Args:
        module_name: Module name inside the package, such as 'chains'.
        title: Section title.

    Returns:
        The section.
    """
    module = importlib.import_module(f"{PACKAGE}.{module_name}")
    lede = (inspect.getdoc(module) or "").split("\n\n")[0]
    public_api = set(synced_structure.__all__)
    defined = [(name, obj) for name, obj in vars(module).items()
               if getattr(obj, "__module__", None) == module.__name__ and (inspect.isclass(obj) or
                                                                         inspect.isfunction(obj))]
    entries: list[Entry] = []
    for name, obj in defined:
        if not inspect.isclass(obj):
            continue
        is_mixin = issubclass(obj, StructureState) and obj is not synced_structure.SyncedStructure
        if obj is StructureState:
            continue
        if is_mixin:
            members = list(_class_members(obj, "SyncedStructure"))
            for entry in members:
                entry.note = f"from {name}"
            entries.extend(members)
        else:
            entries.append(_class_entry(obj))
            entries.extend(_class_members(obj, name))
    functions = []
    for name, obj in defined:
        if inspect.isfunction(obj):
            title_name = name if name in public_api else f"{module_name}.{name}"
            functions.append(Entry(anchor=f"{module_name}.{name}", title=title_name, name=name, kind="function",
                                   obj=obj, signature=_signature(obj), doc=parse_docstring(inspect.getdoc(obj)),
                                   internal=name.startswith("_"), source=_source(obj)))
    entries.extend(sorted(functions, key=lambda entry: entry.internal))
    if module_name == "_types":
        entries.extend(_type_aliases(module))
    exported = sorted(name for name, _ in defined if name in public_api)
    if exported:
        import_line = f"from {PACKAGE} import {', '.join(exported)}"
    elif any(entry.title.startswith("SyncedStructure.") for entry in entries):
        import_line = f"from {PACKAGE} import SyncedStructure"
    else:
        import_line = f"import {PACKAGE}.{module_name}"
    return Section(anchor=f"mod-{module_name.strip('_')}", title=title, module=module_name, lede=lede,
                   import_line=import_line, entries=entries)


def _type_aliases(module: Any) -> list[Entry]:
    """Return one entry per type alias assigned at the top of a module, with its source expression.

    Args:
        module: Imported module.

    Returns:
        The alias entries.
    """
    path = Path(module.__file__)
    tree = ast.parse(path.read_text())
    entries = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name, value = node.targets[0].id, ast.unparse(node.value)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.value is not None:
            name, value = node.target.id, ast.unparse(node.value)
        else:
            continue
        entries.append(Entry(anchor=f"_types.{name}", title=name, name=name, kind="alias", obj=value,
                             signature=None, doc=Doc(), internal=False, source=f"{path.name}:{node.lineno}"))
    return entries


###################################### RENDERING

def esc(text: Any) -> str:
    """Return text escaped for HTML.

    Args:
        text: Any value.

    Returns:
        The escaped string.
    """
    return html.escape(str(text), quote=True)


class Linker:
    """Turn docstring text into HTML, linking known API names and marking code."""

    IDENT = re.compile(r"(?<![\w.#])([A-Za-z_][A-Za-z0-9_]*)(\(\))?(?![\w])")
    QUOTED = re.compile(r"(?<![\w])'([^'\n]{1,48}?)'(?![\w])")
    TICKED = re.compile(r"`([^`\n]+)`")

    def __init__(self, sections: list[Section]) -> None:
        """Index every entry by its bare name, keeping the first one found.

        Args:
            sections: All collected sections.
        """
        self.targets: dict[str, str] = {}
        for section in sections:
            for entry in section.entries:
                linkable = "_" in entry.name.strip("_") or entry.kind in ("class", "typeddict")
                if linkable and entry.kind != "alias":
                    self.targets.setdefault(entry.name, entry.anchor)

    def __call__(self, text: str) -> str:
        """Return docstring text as HTML.

        Args:
            text: Plain docstring text.

        Returns:
            Escaped text with code spans and links.
        """
        out = html.escape(text, quote=False)
        out = self.TICKED.sub(lambda m: f"<code>{m.group(1)}</code>", out)
        out = self.QUOTED.sub(lambda m: f"<code>'{m.group(1)}'</code>", out)
        return self.IDENT.sub(self._link, out)

    def _link(self, match: re.Match[str]) -> str:
        """Return a link for a known API name, or the text unchanged.

        Args:
            match: Identifier match.

        Returns:
            HTML.
        """
        name = match.group(1)
        anchor = self.targets.get(name)
        if anchor is None:
            return match.group(0)
        return f'<a class="xref" href="#{esc(anchor)}"><code>{esc(name)}{match.group(2) or ""}</code></a>'


def _short_repr(value: Any, limit: int = 48) -> str:
    """Return the repr of a default value, shortened past limit characters.

    Args:
        value: Default value.
        limit: Maximum length.

    Returns:
        The repr.
    """
    text = repr(value)
    return text if len(text) <= limit else text[:limit - 1] + "…"


def _annotation(annotation: Any) -> str:
    """Return an annotation as source text.

    Args:
        annotation: Annotation object or string.

    Returns:
        The text, or '' when there is none.
    """
    if annotation is inspect.Parameter.empty or annotation is inspect.Signature.empty:
        return ""
    return annotation if isinstance(annotation, str) else inspect.formatannotation(annotation)


def render_signature(entry: Entry) -> str:
    """Return the highlighted signature block of an entry, one parameter per line when long.

    Args:
        entry: Entry with a signature.

    Returns:
        HTML.
    """
    signature = entry.signature
    if entry.kind == "alias":
        return f'<pre class="sig"><span class="n">{esc(entry.name)}</span> = <span class="t">{esc(entry.obj)}</span></pre>'
    if signature is None:
        return ""
    parts_plain, parts_html = [], []
    seen_kw_marker = False
    for param in signature.parameters.values():
        prefix = ""
        if param.kind is param.VAR_POSITIONAL:
            prefix, seen_kw_marker = "*", True
        elif param.kind is param.VAR_KEYWORD:
            prefix = "**"
        elif param.kind is param.KEYWORD_ONLY and not seen_kw_marker and entry.kind != "typeddict":
            parts_plain.append("*")
            parts_html.append('<span class="o">*</span>')
            seen_kw_marker = True
        plain = f"{prefix}{param.name}"
        rich = f'<span class="o">{prefix}</span><span class="p">{esc(param.name)}</span>'
        annotation = _annotation(param.annotation)
        if annotation:
            plain += f": {annotation}"
            rich += f': <span class="t">{esc(annotation)}</span>'
        if param.default is not param.empty:
            default = _short_repr(param.default)
            plain += f" = {default}"
            rich += f' = <span class="d">{esc(default)}</span>'
        parts_plain.append(plain)
        parts_html.append(rich)
    returns = _annotation(signature.return_annotation)
    if entry.kind in ("class", "typeddict"):
        returns = ""
    name = entry.name if entry.kind not in ("class", "typeddict") else entry.title
    tail_plain = f" -> {returns}" if returns else ""
    tail_html = f' <span class="o">→</span> <span class="t">{esc(returns)}</span>' if returns else ""
    if entry.kind == "property":
        return f'<pre class="sig"><span class="n">{esc(name)}</span>{tail_html}</pre>'
    one_line = f"{name}({', '.join(parts_plain)}){tail_plain}"
    if len(one_line) <= 88 or not parts_html:
        inner = ", ".join(parts_html)
    else:
        inner = "\n" + ",\n".join(f"    {part}" for part in parts_html) + ",\n"
    return f'<pre class="sig"><span class="n">{esc(name)}</span>({inner}){tail_html}</pre>'


def render_entry(entry: Entry, link: Linker) -> str:
    """Return the card of one entry.

    Args:
        entry: Entry to render.
        link: Docstring text renderer.

    Returns:
        HTML.
    """
    doc = entry.doc
    params = entry.signature.parameters if entry.signature is not None else {}
    returns = _annotation(entry.signature.return_annotation) if entry.signature is not None else ""
    if entry.kind in ("class", "typeddict", "alias") or returns == "None":
        returns = ""
    blocks = []
    if doc.summary:
        blocks.append(f'<p class="summary">{link(doc.summary)}</p>')
    blocks.extend(f"<p>{link(paragraph)}</p>" for paragraph in doc.paragraphs)
    if doc.args:
        rows = []
        for name, text in doc.args:
            param = params.get(name.lstrip("*"))
            annotation = _annotation(param.annotation) if param is not None else ""
            default = ("" if param is None or param.default is inspect.Parameter.empty
                       else _short_repr(param.default))
            meta = f'<code class="t">{esc(annotation)}</code>' if annotation else ""
            if default:
                meta += f'<span class="default">default <code>{esc(default)}</code></span>'
            rows.append(f'<tr><td><code class="p">{esc(name)}</code></td><td>{meta}</td><td>{link(text)}</td></tr>')
        blocks.append(_table("Arguments", rows))
    if doc.attributes:
        rows = [f'<tr><td><code class="p">{esc(name)}</code></td><td>{link(text)}</td></tr>'
                for name, text in doc.attributes]
        blocks.append(_table("Attributes", rows, columns=2))
    if doc.returns or returns:
        chip = f'<code class="t">{esc(returns)}</code> ' if returns else ""
        label = "Yields" if entry.kind == "function" and inspect.isgeneratorfunction(entry.obj) else "Returns"
        blocks.append(f'<div class="field"><h4>{label}</h4><p>{chip}{link(doc.returns)}</p></div>')
    if doc.raises:
        items = "".join(f'<li><code class="t">{esc(name)}</code> {link(text)}</li>' for name, text in doc.raises)
        blocks.append(f'<div class="field"><h4>Raises</h4><ul>{items}</ul></div>')
    for name, text in doc.extra:
        blocks.append(f'<div class="field"><h4>{esc(name)}</h4><pre class="block">{esc(text)}</pre></div>')
    source = ""
    if entry.source:
        file, _, line = entry.source.partition(":")
        source = (f'<a class="source" href="{SOURCE_URL}/{esc(file)}#L{esc(line)}" target="_blank" '
                  f'rel="noopener">{esc(entry.source)}</a>')
    note = f'<span class="note">{esc(entry.note)}</span>' if entry.note else ""
    chip = f'<span class="returns" title="Return type">→ {esc(returns)}</span>' if returns else ""
    search = esc(f"{entry.title} {doc.summary}".lower())
    classes = "entry internal" if entry.internal else "entry"
    if entry.kind in ("class", "typeddict"):
        classes += " is-class"
    return f"""
<article class="{classes}" id="{esc(entry.anchor)}" data-search="{search}">
  <header>
    <span class="kind k-{entry.kind}">{entry.kind}</span>
    <h3><a href="#{esc(entry.anchor)}">{esc(entry.title)}</a></h3>
    {chip}
  </header>
  {render_signature(entry)}
  <div class="body">{''.join(blocks)}</div>
  <footer>{note}{source}</footer>
</article>"""


def _table(title: str, rows: list[str], columns: int = 3) -> str:
    """Return a titled definition table.

    Args:
        title: Field title.
        rows: Rendered table rows.
        columns: Number of columns, used for the header.

    Returns:
        HTML.
    """
    head = "<th>Name</th><th>Type</th><th>Description</th>" if columns == 3 else "<th>Name</th><th>Description</th>"
    return (f'<div class="field"><h4>{title}</h4><div class="table-wrap"><table class="params">'
            f'<thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div></div>')


def render_section(section: Section, link: Linker) -> str:
    """Return one module section.

    Args:
        section: Section to render.
        link: Docstring text renderer.

    Returns:
        HTML.
    """
    public = sum(not entry.internal for entry in section.entries)
    internal = len(section.entries) - public
    counts = f'<span class="count">{public} public</span>'
    if internal:
        counts += f'<span class="count internal-count">{internal} internal</span>'
    cards = "".join(render_entry(entry, link) for entry in section.entries)
    return f"""
<section class="module" id="{section.anchor}">
  <div class="module-head">
    <p class="eyebrow">{esc(PACKAGE)}.{esc(section.module)}</p>
    <h2>{esc(section.title)}</h2>
    <p class="lede">{link(section.lede)}</p>
    <div class="module-meta"><code class="import">{esc(section.import_line)}</code>{counts}</div>
  </div>
  {cards}
</section>"""


###################################### HAND-WRITTEN SECTIONS

def render_guides(sections: list[Section], link: Linker) -> str:
    """Return the overview, conventions, concepts and configuration sections.

    Args:
        sections: All module sections, for counts and entry points.
        link: Docstring text renderer.

    Returns:
        HTML.
    """
    entries = [entry for section in sections for entry in section.entries]
    public = sum(not entry.internal for entry in entries)
    internal = len(entries) - public
    starts = "".join(
        f'<a class="tile" href="#{esc(anchor)}"><code>{esc(anchor.split(".", 1)[1])}</code>'
        f'<span>{esc(next(e.doc.summary for e in entries if e.anchor == anchor))}</span></a>'
        for anchor in KEY_ENTRY_POINTS)
    config_rows = "".join(
        f'<tr><td><code class="p">{esc(key)}</code></td><td><code>{esc(repr(value))}</code></td>'
        f'<td>{link(CONFIG_HELP.get(key, ""))}</td></tr>' for key, value in DEFAULT_CONFIG.items())
    yaml = "\n".join(["synced_structure:"] + [f"  {key}: {value}" for key, value in DEFAULT_CONFIG.items()])
    return f"""
<section class="guide hero" id="overview">
  <p class="eyebrow">API reference · v{esc(synced_structure.__version__)}</p>
  <h1>{esc(PACKAGE)}</h1>
  <p class="lede">{esc(inspect.getdoc(synced_structure) or "")}</p>
  <p>A Biopython structure wrapper that keeps every observed residue linked to its position in the reference
  sequence, and builds contacts, contact maps, interface metrics, HADDOCK3 dockings, tables and 3D views on top
  of that link. Installation and a quick start are in the
  <a href="{REPO}#readme">README</a>.</p>
  <div class="stats">
    <div><strong>{public}</strong><span>public names</span></div>
    <div><strong>{internal}</strong><span>internal helpers</span></div>
    <div><strong>{len(sections)}</strong><span>modules</span></div>
  </div>
  <pre class="block"><span class="c"># install</span>
pip install "git+{REPO}.git"

<span class="c"># load a structure and look at it</span>
<span class="k">from</span> synced_structure <span class="k">import</span> SyncedStructure
s = SyncedStructure.from_rcsb(<span class="s">"1AAY"</span>)
s.chain_table()</pre>
  <h3 class="sub">Start here</h3>
  <div class="tiles">{starts}</div>
</section>

<section class="guide" id="conventions">
  <p class="eyebrow">Read once</p>
  <h2>Conventions</h2>
  <h3 class="sub">Names say what they return</h3>
  <p>A name starts with a verb only when the verb fixes the return type or a side effect. Otherwise the name is
  a noun whose last word is the type.</p>
  <div class="two-col">
    <div class="table-wrap"><table class="params plain">
      <thead><tr><th>Verb</th><th>Returns</th></tr></thead><tbody>
      <tr><td><code>set_</code>, <code>add_</code>, <code>remove_</code>, <code>rename_</code>, <code>remap_</code>, <code>apply_</code>, <code>clear_</code></td><td>the object itself (chainable)</td></tr>
      <tr><td><code>subset_</code>, <code>from_</code>, <code>load_</code></td><td>a new <code>SyncedStructure</code></td></tr>
      <tr><td><code>find_</code>, <code>fetch_</code></td><td>a DataFrame of matches / a DataFrame from the network</td></tr>
      <tr><td><code>search_</code></td><td>a list of IDs</td></tr>
      <tr><td><code>plot_</code>, <code>show_</code></td><td><code>(Figure, Axes)</code> / IPython HTML</td></tr>
      <tr><td><code>export_&lt;name&gt;</code></td><td>the path(s) written for <code>&lt;name&gt;</code></td></tr>
      <tr><td><code>run_</code></td><td>runs a haddock3 subprocess</td></tr>
      <tr><td><code>get_</code></td><td>a stored value, paired with <code>set_</code></td></tr>
    </tbody></table></div>
    <div class="table-wrap"><table class="params plain">
      <thead><tr><th>Noun ending</th><th>Type</th></tr></thead><tbody>
      <tr><td><code>_table</code></td><td><code>pandas.DataFrame</code></td></tr>
      <tr><td><code>_tables</code></td><td><code>dict[str, DataFrame]</code></td></tr>
      <tr><td><code>_matrix</code></td><td>dict holding a matrix and its labels</td></tr>
      <tr><td><code>_map</code></td><td><code>dict</code></td></tr>
      <tr><td><code>_ids</code></td><td><code>list[str]</code></td></tr>
      <tr><td><code>_paths</code></td><td>file or folder paths</td></tr>
      <tr><td><code>_string</code>, <code>_text</code>, <code>_toml</code></td><td><code>str</code></td></tr>
    </tbody></table></div>
  </div>
  <h3 class="sub">Numbers and units</h3>
  <ul class="facts">
    <li><strong>Indices are 1-based</strong> in every returned table.</li>
    <li><strong>nm</strong> for contacts: <code>contact_table</code>, contact matrices, contact maps and every contact
    <code>cutoff=</code>.</li>
    <li><strong>Å</strong> and <strong>Å²</strong> for interface metrics and docking; names say <code>_a</code> or
    <code>_A</code> where it matters.</li>
    <li>A chain can be named by its <strong>key</strong>, <strong>label ID</strong> or <strong>author ID</strong>; keys
    are matched first.</li>
  </ul>
</section>

<section class="guide" id="concepts">
  <p class="eyebrow">Vocabulary</p>
  <h2>Concepts</h2>
  <p>A deposited structure is an incomplete copy of a biological sequence: some residues are not modelled,
  numbering is chosen by the depositors, and every chain has two IDs. Every table reports all of these names.</p>
  <div class="numbers">
    <div><span class="big">125</span><code>pdb_resnum</code><span>author number in the file</span></div>
    <div><span class="big">1</span><code>sequence_index</code><span>position in the observed sequence</span></div>
    <div><span class="big">6</span><code>reference_index</code><span>position in the reference sequence</span></div>
  </div>
  <p class="caption">The first modelled residue of chain A of 5D39 has three numbers. Only
  <code>reference_index</code> means the same residue in every chain, structure and model.</p>
  <div class="table-wrap"><table class="params glossary"><thead><tr><th>Term</th><th>Meaning</th></tr></thead><tbody>
    <tr><td>chain key</td><td>The name a chain is stored under: its author ID by default (<code>key_by="auth"</code>) or its label ID (<code>key_by="label"</code>).</td></tr>
    <tr><td>label ID <code>chain_id</code></td><td>Chain ID assigned by the PDB (<code>label_asym_id</code>).</td></tr>
    <tr><td>author ID <code>auth_chain_id</code></td><td>Chain ID chosen by the depositors (<code>auth_asym_id</code>); the one most viewers show. The same letter can name a different chain in the other system.</td></tr>
    <tr><td>observed sequence</td><td>One-letter sequence of the residues present in the file.</td></tr>
    <tr><td>reference sequence</td><td>Full sequence of the chain's entity, from the mmCIF file or set by you.</td></tr>
    <tr><td>mapping</td><td>Observed → reference correspondence of a chain, built by global alignment and rebuilt whenever a chain, reference or numbering changes.</td></tr>
    <tr><td>moiety</td><td>Type of a residue: <code>protein</code>, <code>nucleic</code> or <code>other</code> (ligands, ions, water). Polymer residues are protein and nucleic.</td></tr>
    <tr><td>contact point</td><td>One point per residue for the <code>ca</code> scheme: Cα (protein), C4′ (nucleic), heavy-atom centroid (other).</td></tr>
    <tr><td>contact scheme</td><td><code>ca</code>: distance between contact points. <code>sidechain</code>: closest side-chain (protein) or base (nucleic) heavy atoms.</td></tr>
    <tr><td>contact type</td><td><code>inter</code> (between chains), <code>intra</code> (within a chain) or <code>all</code>.</td></tr>
    <tr><td>contact string</td><td>One character per residue: <code>1</code> contact, <code>0</code> none, <code>X</code> no contact point.</td></tr>
    <tr><td>interface</td><td>Residues of two chains that bury surface against each other, measured with heavy atoms.</td></tr>
    <tr><td>partner</td><td>One rigid body of a HADDOCK3 docking; it may hold several chains.</td></tr>
  </tbody></table></div>
</section>

<section class="guide" id="configuration">
  <p class="eyebrow">Optional</p>
  <h2>Configuration</h2>
  <p>Every loader accepts <code>config_path=</code>, a YAML file whose <code>synced_structure</code> section
  overrides these defaults. The file is never written.</p>
  <div class="table-wrap"><table class="params"><thead><tr><th>Key</th><th>Default</th><th>Meaning</th></tr></thead>
  <tbody>{config_rows}</tbody></table></div>
  <pre class="block">{esc(yaml)}</pre>
</section>"""


###################################### PAGE

def build_html() -> str:
    """Return the whole reference page.

    Returns:
        The HTML document.
    """
    groups = [(group, [collect_section(module, title) for module, title in modules]) for group, modules in LAYOUT]
    sections = [section for _, group_sections in groups for section in group_sections]
    link = Linker(sections)
    nav = ['<div class="nav-group"><p>Guide</p>']
    nav += [f'<a href="#{anchor}" data-target="{anchor}">{esc(title)}</a>' for anchor, title in GUIDE_SECTIONS]
    nav.append("</div>")
    for group, group_sections in groups:
        nav.append(f'<div class="nav-group"><p>{esc(group)}</p>')
        for section in group_sections:
            nav.append(f'<a href="#{section.anchor}" data-target="{section.anchor}">{esc(section.title)}'
                       f'<span>{sum(not e.internal for e in section.entries)}</span></a>')
        nav.append("</div>")
    body = render_guides(sections, link) + "".join(render_section(section, link) for section in sections)
    return PAGE.format(
        title=f"{PACKAGE} API",
        version=esc(synced_structure.__version__),
        repo=REPO,
        nav="\n".join(nav),
        body=body,
        css=CSS,
        js=JS,
        package=esc(PACKAGE),
    )


def main(argv: list[str] | None = None) -> Path:
    """Write the reference page.

    Args:
        argv: Command-line arguments, or None for sys.argv.

    Returns:
        The written path.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-o", "--output", type=Path, default=OUTPUT, help=f"output file (default {OUTPUT})")
    output = parser.parse_args(argv).output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_html(), encoding="utf-8")
    print(f"wrote {output}")
    return output


CSS = """
:root {
  --bg: #f6f7f9; --surface: #ffffff; --surface-2: #f0f2f5; --text: #15171c; --muted: #5d6370; --faint: #8a909c;
  --border: #e3e6eb; --accent: #4f46e5; --accent-soft: #eef0ff; --accent-text: #3730a3;
  --type: #0f766e; --param: #9a3412; --default: #6d28d9; --code-bg: #f3f4f7;
  --k-method: #4f46e5; --k-property: #0d9488; --k-classmethod: #b45309; --k-staticmethod: #64748b;
  --k-function: #2563eb; --k-class: #be185d; --k-typeddict: #be185d; --k-alias: #475569;
  --shadow: 0 1px 2px rgba(16, 24, 40, .04), 0 1px 3px rgba(16, 24, 40, .06);
  --radius: 12px; --header: 60px;
  --sans: "Inter", ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
  --mono: "JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  color-scheme: light;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #0d0f13; --surface: #151820; --surface-2: #1b1f29; --text: #e7e9ee; --muted: #a0a6b3; --faint: #737a88;
    --border: #262b36; --accent: #8b8cf9; --accent-soft: #1f2140; --accent-text: #c7c8ff;
    --type: #5eead4; --param: #fdba74; --default: #c4b5fd; --code-bg: #1b1f29;
    --k-method: #8b8cf9; --k-property: #2dd4bf; --k-classmethod: #fbbf24; --k-staticmethod: #94a3b8;
    --k-function: #60a5fa; --k-class: #f472b6; --k-typeddict: #f472b6; --k-alias: #94a3b8;
    --shadow: none; color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #0d0f13; --surface: #151820; --surface-2: #1b1f29; --text: #e7e9ee; --muted: #a0a6b3; --faint: #737a88;
  --border: #262b36; --accent: #8b8cf9; --accent-soft: #1f2140; --accent-text: #c7c8ff;
  --type: #5eead4; --param: #fdba74; --default: #c4b5fd; --code-bg: #1b1f29;
  --k-method: #8b8cf9; --k-property: #2dd4bf; --k-classmethod: #fbbf24; --k-staticmethod: #94a3b8;
  --k-function: #60a5fa; --k-class: #f472b6; --k-typeddict: #f472b6; --k-alias: #94a3b8;
  --shadow: none; color-scheme: dark;
}
* { box-sizing: border-box; }
html { scroll-padding-top: calc(var(--header) + 16px); }
body { margin: 0; background: var(--bg); color: var(--text); font: 15px/1.6 var(--sans);
  -webkit-font-smoothing: antialiased; }
a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }
code, pre { font-family: var(--mono); font-size: .86em; }
code { background: var(--code-bg); padding: .1em .35em; border-radius: 5px; }
a code { color: inherit; }

.topbar { position: sticky; top: 0; z-index: 20; height: var(--header); display: flex; align-items: center;
  gap: 16px; padding: 0 20px; background: color-mix(in srgb, var(--bg) 82%, transparent);
  backdrop-filter: saturate(160%) blur(12px); border-bottom: 1px solid var(--border); }
.brand { display: flex; align-items: center; gap: 10px; font-weight: 650; color: var(--text); white-space: nowrap; }
.brand:hover { text-decoration: none; }
.logo { width: 28px; height: 28px; border-radius: 8px; display: grid; place-items: center; color: #fff;
  background: linear-gradient(135deg, #6366f1, #14b8a6); font: 700 13px var(--mono); }
.version { font: 500 12px var(--mono); color: var(--muted); border: 1px solid var(--border); padding: 2px 8px;
  border-radius: 999px; }
.search { flex: 1; max-width: 520px; margin-left: auto; position: relative; }
.search input { width: 100%; height: 38px; padding: 0 64px 0 38px; border-radius: 10px; border: 1px solid var(--border);
  background: var(--surface); color: var(--text); font: inherit; outline: none; }
.search input:focus { border-color: var(--accent); box-shadow: 0 0 0 3px var(--accent-soft); }
.search svg { position: absolute; left: 12px; top: 11px; color: var(--faint); }
.search kbd { position: absolute; right: 10px; top: 9px; font: 11px var(--mono); color: var(--faint);
  border: 1px solid var(--border); border-radius: 5px; padding: 1px 6px; }
.controls { display: flex; gap: 8px; align-items: center; }
.btn { height: 38px; padding: 0 12px; border-radius: 10px; border: 1px solid var(--border); background: var(--surface);
  color: var(--text); font: 500 13px var(--sans); cursor: pointer; display: inline-flex; align-items: center; gap: 8px;
  white-space: nowrap; }
.btn:hover { border-color: var(--accent); }
.btn[aria-pressed="true"] { background: var(--accent-soft); border-color: var(--accent); color: var(--accent-text); }
.menu { display: none; }

.layout { display: grid; grid-template-columns: 270px minmax(0, 1fr); max-width: 1320px; margin: 0 auto; }
.sidebar { position: sticky; top: var(--header); height: calc(100vh - var(--header)); overflow-y: auto;
  padding: 24px 16px 48px 20px; border-right: 1px solid var(--border); }
.nav-group { margin-bottom: 18px; }
.nav-group p { margin: 0 0 6px 10px; font: 600 11px var(--sans); letter-spacing: .08em; text-transform: uppercase;
  color: var(--faint); }
.nav-group a { display: flex; justify-content: space-between; align-items: center; padding: 5px 10px;
  border-radius: 8px; color: var(--muted); font-size: 14px; }
.nav-group a span { font: 11px var(--mono); color: var(--faint); }
.nav-group a:hover { background: var(--surface-2); color: var(--text); text-decoration: none; }
.nav-group a.active { background: var(--accent-soft); color: var(--accent-text); font-weight: 550; }
.nav-group a.empty { opacity: .35; }

main { padding: 32px 40px 120px; min-width: 0; }
.content { max-width: 940px; }
section { margin-bottom: 64px; }
.eyebrow { margin: 0 0 6px; font: 600 12px var(--mono); color: var(--accent); letter-spacing: .02em; }
h1 { font-size: 44px; line-height: 1.1; letter-spacing: -.03em; margin: 0 0 12px; }
h2 { font-size: 28px; letter-spacing: -.02em; margin: 0 0 8px; }
h3.sub { font-size: 17px; margin: 32px 0 10px; }
.lede { color: var(--muted); font-size: 17px; margin: 0 0 16px; }
.hero h1 { background: linear-gradient(120deg, var(--text), var(--accent)); -webkit-background-clip: text;
  background-clip: text; color: transparent; }
.stats { display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; margin: 24px 0; }
.stats div { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius); padding: 14px 16px;
  box-shadow: var(--shadow); }
.stats strong { display: block; font-size: 26px; letter-spacing: -.02em; }
.stats span { color: var(--muted); font-size: 13px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 10px; }
.tile { display: flex; flex-direction: column; gap: 4px; padding: 12px 14px; border: 1px solid var(--border);
  border-radius: var(--radius); background: var(--surface); color: var(--muted); font-size: 13px; box-shadow: var(--shadow);
  transition: border-color .15s, transform .15s; }
.tile code { background: none; padding: 0; color: var(--accent-text); font-size: 13px; font-weight: 600; }
.tile:hover { border-color: var(--accent); transform: translateY(-1px); text-decoration: none; }
.two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
.facts { padding-left: 20px; }
.facts li { margin: 6px 0; }
.numbers { display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; margin: 20px 0 8px; }
.numbers div { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius);
  padding: 16px; display: flex; flex-direction: column; gap: 4px; box-shadow: var(--shadow); }
.numbers .big { font: 700 30px var(--mono); color: var(--accent); }
.numbers code { align-self: flex-start; }
.numbers span:last-child { color: var(--muted); font-size: 13px; }
.caption { color: var(--muted); font-size: 14px; }
pre.block { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius);
  padding: 14px 16px; overflow-x: auto; line-height: 1.55; }
pre.block .c { color: var(--faint); } pre.block .k { color: var(--accent); } pre.block .s { color: var(--type); }

.module-head { padding-bottom: 16px; margin-bottom: 20px; border-bottom: 1px solid var(--border); }
.module-meta { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
.import { font-size: 12.5px; }
.count { font: 500 12px var(--sans); color: var(--muted); background: var(--surface-2); padding: 2px 9px;
  border-radius: 999px; }
.internal-count { display: none; }
body.show-internal .internal-count { display: inline; }

.entry { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius);
  padding: 16px 20px; margin: 0 0 14px; box-shadow: var(--shadow); }
.entry:target { border-color: var(--accent); box-shadow: 0 0 0 3px var(--accent-soft); }
.entry.is-class { border-left: 3px solid var(--k-class); }
.entry.internal { display: none; background: var(--surface-2); border-style: dashed; }
body.show-internal .entry.internal { display: block; }
.entry.hidden { display: none !important; }
.entry header { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
.entry h3 { margin: 0; font: 600 15.5px var(--mono); word-break: break-word; }
.entry h3 a { color: var(--text); }
.kind { font: 600 10.5px var(--sans); text-transform: uppercase; letter-spacing: .06em; padding: 2px 8px;
  border-radius: 6px; color: var(--k); background: color-mix(in srgb, var(--k) 13%, transparent); }
.k-method { --k: var(--k-method); } .k-property { --k: var(--k-property); } .k-classmethod { --k: var(--k-classmethod); }
.k-staticmethod { --k: var(--k-staticmethod); } .k-function { --k: var(--k-function); } .k-class { --k: var(--k-class); }
.k-typeddict { --k: var(--k-typeddict); } .k-alias { --k: var(--k-alias); }
.returns { margin-left: auto; font: 500 12px var(--mono); color: var(--type); background: color-mix(in srgb, var(--type) 10%, transparent);
  padding: 2px 9px; border-radius: 999px; max-width: 100%; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
pre.sig { margin: 12px 0 4px; padding: 10px 14px; background: var(--code-bg); border-radius: 9px; overflow-x: auto;
  line-height: 1.55; font-size: 13px; }
pre.sig .n { color: var(--text); font-weight: 600; } pre.sig .p { color: var(--param); } pre.sig .t { color: var(--type); }
pre.sig .d { color: var(--default); } pre.sig .o { color: var(--faint); }
.body p { margin: 10px 0; }
.body .summary { font-size: 15.5px; }
.field h4 { margin: 16px 0 6px; font: 600 11px var(--sans); text-transform: uppercase; letter-spacing: .08em; color: var(--faint); }
.field ul { margin: 0; padding-left: 18px; }
.field p { margin: 0; }
.table-wrap { overflow-x: auto; }
table.params { width: 100%; border-collapse: collapse; font-size: 14px; }
table.params th { text-align: left; font: 600 11px var(--sans); text-transform: uppercase; letter-spacing: .06em;
  color: var(--faint); padding: 6px 10px; border-bottom: 1px solid var(--border); }
table.params td { padding: 8px 10px; border-bottom: 1px solid var(--border); vertical-align: top; }
table.params tr:last-child td { border-bottom: none; }
table.params td:first-child { white-space: nowrap; }
code.p { color: var(--param); } code.t { color: var(--type); }
.default { display: block; margin-top: 4px; color: var(--faint); font-size: 12px; white-space: nowrap; }
.default code { color: var(--default); }
table.params.plain td { white-space: normal; width: 50%; }
table.params.plain td code { white-space: nowrap; }
.glossary td:first-child { font-weight: 600; white-space: normal; min-width: 140px; }
.entry footer { display: flex; justify-content: space-between; gap: 12px; margin-top: 10px; font: 12px var(--mono); }
.entry footer:empty { display: none; }
.note { color: var(--faint); }
.source { color: var(--faint); margin-left: auto; }
.source:hover { color: var(--accent); }
a.xref code { color: var(--accent-text); background: var(--accent-soft); }
.no-results { display: none; padding: 40px; text-align: center; color: var(--muted); border: 1px dashed var(--border);
  border-radius: var(--radius); }
body.searching section.guide { display: none; }
body.searching .module.no-match { display: none; }

@media (max-width: 960px) {
  .layout { grid-template-columns: 1fr; }
  .sidebar { position: fixed; inset: var(--header) 0 0 0; z-index: 15; background: var(--bg); display: none; border: none; }
  body.nav-open .sidebar { display: block; }
  .menu { display: inline-flex; }
  main { padding: 24px 16px 96px; }
  .two-col, .numbers, .stats { grid-template-columns: 1fr; }
  h1 { font-size: 34px; }
}
@media (max-width: 640px) {
  .topbar { gap: 8px; padding: 0 12px; }
  .brand .name, .version, .btn .label, .search kbd { display: none; }
  .returns { margin-left: 0; }
}
"""

JS = """
(() => {
  const root = document.documentElement, body = document.body;
  const store = {
    get(key) { try { return localStorage.getItem(key); } catch (e) { return null; } },
    set(key, value) { try { localStorage.setItem(key, value); } catch (e) {} },
  };
  const theme = store.get("ss-theme");
  if (theme) root.dataset.theme = theme;
  document.getElementById("theme").addEventListener("click", () => {
    const dark = root.dataset.theme ? root.dataset.theme === "dark"
      : matchMedia("(prefers-color-scheme: dark)").matches;
    root.dataset.theme = dark ? "light" : "dark";
    store.set("ss-theme", root.dataset.theme);
  });

  const internal = document.getElementById("internal");
  const setInternal = (on) => {
    body.classList.toggle("show-internal", on);
    internal.setAttribute("aria-pressed", String(on));
    store.set("ss-internal", on ? "1" : "0");
    filter();
  };
  internal.addEventListener("click", () => setInternal(!body.classList.contains("show-internal")));

  document.getElementById("menu").addEventListener("click", () => body.classList.toggle("nav-open"));
  document.querySelectorAll(".sidebar a").forEach((a) => a.addEventListener("click", () => body.classList.remove("nav-open")));

  const input = document.getElementById("q");
  const entries = [...document.querySelectorAll(".entry")];
  const modules = [...document.querySelectorAll(".module")];
  const empty = document.getElementById("no-results");
  function filter() {
    const terms = input.value.toLowerCase().trim().split(/\\s+/).filter(Boolean);
    body.classList.toggle("searching", terms.length > 0);
    const showInternal = body.classList.contains("show-internal");
    let shown = 0;
    entries.forEach((el) => {
      const hit = terms.every((t) => el.dataset.search.includes(t));
      el.classList.toggle("hidden", terms.length > 0 && !hit);
      if (hit && (showInternal || !el.classList.contains("internal"))) shown++;
    });
    modules.forEach((m) => {
      const any = [...m.querySelectorAll(".entry")].some((el) => !el.classList.contains("hidden")
        && (showInternal || !el.classList.contains("internal")));
      m.classList.toggle("no-match", !any);
      const link = document.querySelector(`.sidebar a[data-target="${m.id}"]`);
      if (link) link.classList.toggle("empty", terms.length > 0 && !any);
    });
    empty.style.display = terms.length && !shown ? "block" : "none";
  }
  input.addEventListener("input", filter);
  input.addEventListener("keydown", (e) => { if (e.key === "Escape") { input.value = ""; filter(); input.blur(); } });
  document.addEventListener("keydown", (e) => {
    if ((e.key === "/" || (e.key === "k" && (e.metaKey || e.ctrlKey))) && document.activeElement !== input) {
      e.preventDefault(); input.focus(); input.select();
    }
  });
  if (location.hash) {
    const target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
    if (target && target.classList.contains("internal")) body.classList.add("show-internal");
  }
  setInternal(store.get("ss-internal") === "1" || body.classList.contains("show-internal"));

  const links = new Map([...document.querySelectorAll(".sidebar a[data-target]")].map((a) => [a.dataset.target, a]));
  const observer = new IntersectionObserver((items) => {
    items.forEach((item) => {
      if (!item.isIntersecting) return;
      links.forEach((a) => a.classList.remove("active"));
      const a = links.get(item.target.id);
      if (a) { a.classList.add("active"); a.scrollIntoView({ block: "nearest" }); }
    });
  }, { rootMargin: "-20% 0px -70% 0px" });
  document.querySelectorAll("main section[id]").forEach((s) => observer.observe(s));
})();
"""

PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<meta name="description" content="API reference of {package}: every function, method and property with its signature.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600;700&display=swap">
<style>{css}</style>
</head>
<body>
<header class="topbar">
  <button class="btn menu" id="menu" aria-label="Menu">☰</button>
  <a class="brand" href="#overview"><span class="logo">ss</span><span class="name">{package}</span></a>
  <span class="version">v{version}</span>
  <label class="search">
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true"><circle cx="11" cy="11" r="7"/><path d="m20 20-3.5-3.5"/></svg>
    <input id="q" type="search" placeholder="Search names and summaries…" autocomplete="off" aria-label="Search">
    <kbd>/</kbd>
  </label>
  <div class="controls">
    <button class="btn" id="internal" aria-pressed="false" title="Show internal helpers"><span>⚙</span><span class="label">Internal</span></button>
    <button class="btn" id="theme" title="Toggle theme" aria-label="Toggle theme">◐</button>
    <a class="btn" href="{repo}" target="_blank" rel="noopener"><span class="label">GitHub</span><span>↗</span></a>
  </div>
</header>
<div class="layout">
  <nav class="sidebar">{nav}</nav>
  <main><div class="content">
    {body}
    <div class="no-results" id="no-results">No name matches this search. Internal helpers are searched only when shown.</div>
  </div></main>
</div>
<script>{js}</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
