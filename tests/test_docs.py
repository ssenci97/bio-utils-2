"""The documentation works: the README's code runs, every relative link resolves, and the API page is current.

The README runs like a notebook: its ```python blocks are executed top to bottom in one namespace. A block
preceded by <!-- skip-test --> is not run. Structures downloaded with from_rcsb are loaded from tests/data
instead, so the test is offline. docs/index.html must equal what docs/build_docs.py generates now.
"""

import importlib.util
import re
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from synced_structure import SyncedStructure  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"
GUIDES = [ROOT / "README.md"]
PAGES = [ROOT / "README.md", ROOT / "CHANGELOG.md", *sorted((ROOT / "docs").glob("*.md"))]
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


@pytest.mark.parametrize("guide", GUIDES, ids=[path.stem for path in GUIDES])
def test_guide_code_runs(guide, tmp_path, monkeypatch):
    local = classmethod(lambda cls, pdb_id, **kwargs: cls.from_path(DATA / f"{pdb_id.upper()}.cif", **kwargs))
    monkeypatch.setattr(SyncedStructure, "from_rcsb", local)
    monkeypatch.chdir(tmp_path)
    blocks = _blocks(guide)
    assert blocks, f"{guide.name} has no runnable code"
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


def _build_docs():
    spec = importlib.util.spec_from_file_location("build_docs", ROOT / "docs" / "build_docs.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_api_page_is_current():
    page = _build_docs().build_html()
    assert page == (ROOT / "docs" / "index.html").read_text(encoding="utf-8"), (
        "docs/index.html is stale: run python docs/build_docs.py")


def test_api_page_lists_every_public_name():
    page = _build_docs().build_html()
    import synced_structure
    from synced_structure import HaddockDocking

    missing = [name for name in synced_structure.__all__ if f'id="{name}"' not in page and f'.{name}"' not in page]
    for cls, owner in ((SyncedStructure, "SyncedStructure"), (HaddockDocking, "HaddockDocking")):
        for name in dir(cls):
            member = getattr(cls, name)
            if not name.startswith("_") and (callable(member) or isinstance(member, property)):
                if f'id="{owner}.{name}"' not in page:
                    missing.append(f"{owner}.{name}")
    assert not missing, missing
