"""The documentation works: guide code runs, links resolve, the MkDocs site builds, and the API reference is complete.

Each guide page runs like a notebook: its ```python blocks are executed top to bottom in one namespace. A block
preceded by <!-- skip-test --> is not run. Structures downloaded with from_rcsb are loaded from tests/data
instead, so the test is offline. The API reference must document every public name exactly once.
"""

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from synced_structure import SyncedStructure  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"
DOCS = ROOT / "docs"
GUIDES = [ROOT / "README.md", DOCS / "index.md", DOCS / "quickstart.md", *sorted((DOCS / "guide").glob("*.md")),
          *sorted((DOCS / "tutorials").rglob("*.md"))]
PAGES = [ROOT / "README.md", ROOT / "CHANGELOG.md", *sorted(DOCS.rglob("*.md"))]
DIRECTIVE = re.compile(r"^::: (synced_structure\.[\w.]+)", re.M)
BLOCK = re.compile(r"(<!-- skip-test -->\s*)?```python\n(.*?)```", re.S)
LINK = re.compile(r"\]\(([^)\s]+)\)")


def _blocks(path):
    """Return the code of every ```python block of a page that is not marked skip-test."""
    return [code for skip, code in BLOCK.findall(path.read_text()) if not skip]


def _slug(heading):
    """Return the GitHub anchor of a Markdown heading."""
    text = re.sub(r"[`*_]", "", heading.strip().lower())
    return re.sub(r"[^\w\- ]", "", text).replace(" ", "-")


def _anchors(path):
    return {_slug(line.lstrip("#")) for line in path.read_text().splitlines() if line.startswith("#")}


@pytest.mark.parametrize("guide", GUIDES, ids=[str(path.relative_to(ROOT)) for path in GUIDES])
def test_guide_code_runs(guide, tmp_path, monkeypatch):
    local = classmethod(lambda cls, pdb_id, **kwargs: cls.from_path(DATA / f"{pdb_id.upper()}.cif", **kwargs))
    monkeypatch.setattr(SyncedStructure, "from_rcsb", local)
    monkeypatch.chdir(tmp_path)
    blocks = _blocks(guide)
    if not blocks:
        pytest.skip(f"{guide.name} has no runnable code")
    namespace = {"__name__": f"guide_{guide.stem}"}
    for number, code in enumerate(blocks, start=1):
        try:
            exec(compile(code, f"{guide.name} block {number}", "exec"), namespace)
        except Exception as error:
            raise AssertionError(f"{guide.name}, code block {number} failed:\n{code}") from error
        finally:
            plt.close("all")


def test_relative_links_resolve():
    broken = []
    for page in PAGES:
        for target in LINK.findall(page.read_text()):
            if re.match(r"[a-z]+:", target):
                continue
            file_part, _, anchor = target.partition("#")
            linked = (page.parent / file_part).resolve() if file_part else page
            if not linked.exists():
                broken.append(f"{page.relative_to(ROOT)}: {target} (missing file)")
            elif anchor and linked.suffix == ".md" and anchor not in _anchors(linked):
                broken.append(f"{page.relative_to(ROOT)}: {target} (missing anchor)")
    assert not broken, "\n".join(broken)


def _documented():
    """Return every object path named by a ::: directive in the API reference, with repeats."""
    return [path for page in sorted((DOCS / "reference").glob("*.md")) for path in DIRECTIVE.findall(page.read_text())]


def test_api_reference_lists_every_public_name_once():
    import synced_structure
    from synced_structure import HaddockDocking

    documented = _documented()
    repeated = sorted({path for path in documented if documented.count(path) > 1})
    assert not repeated, f"documented more than once: {repeated}"
    expected = {f"synced_structure.{name}" for name in synced_structure.__all__}
    for name in dir(SyncedStructure):
        member = getattr(SyncedStructure, name)
        if not name.startswith("_") and (callable(member) or isinstance(member, property)):
            expected.add(f"synced_structure.SyncedStructure.{name}")
    # HaddockDocking is documented as a whole class, members included
    expected -= {f"synced_structure.HaddockDocking.{name}" for name in dir(HaddockDocking)}
    missing = sorted(expected - set(documented))
    assert not missing, f"add these to a docs/reference/*.md page: {missing}"
    unknown = sorted(set(documented) - expected)
    assert not unknown, f"documented but not public: {unknown}"


def test_configuration_page_lists_every_setting():
    from synced_structure.constants import DEFAULT_CONFIG

    page = (DOCS / "guide" / "configuration.md").read_text()
    missing = [key for key in DEFAULT_CONFIG if f"| `{key}` | `{DEFAULT_CONFIG[key]!r}` |" not in page]
    assert not missing, f"missing or stale in docs/guide/configuration.md: {missing}"


@pytest.mark.skipif(importlib.util.find_spec("mkdocs") is None or importlib.util.find_spec("mkdocstrings") is None,
                    reason="mkdocs and mkdocstrings are not installed (pip install -e '.[docs]')")
def test_mkdocs_builds_strict(tmp_path):
    result = subprocess.run([sys.executable, "-m", "mkdocs", "build", "--strict", "--site-dir", str(tmp_path / "site")],
                            cwd=ROOT, capture_output=True, text=True, env={**__import__("os").environ,
                                                                          "NO_MKDOCS_2_WARNING": "1"})
    assert result.returncode == 0, result.stderr[-3000:]
