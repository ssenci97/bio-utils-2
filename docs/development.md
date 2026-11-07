# Development

How the code is written and checked. Read this before changing the package.

## Set up

```bash
git clone https://github.com/ssenci97/bio-utils-2.git
cd bio-utils-2
pip install -e ".[dev]"            # editable install + pytest, ruff, pydoclint, mypy, typeguard, type stubs, mkdocs
pip install -e ".[dev,haddock]"    # also haddock3, for the docking tests
```

Always use the editable install (`-e`): tests and scripts then import the code in `src/`, not an old copy in `site-packages`.

## Conventions

- **Docstrings:** [Google style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) on every function and method: a summary line, then `Args:`, `Returns:` and `Raises:` where they apply. Keep them short; types belong in the signature, not the docstring.
- **Type hints:** every parameter and return is annotated (not `self`/`cls`). Methods that change the structure return `Self`. Shared aliases live in `src/synced_structure/_types.py`.
- **Mixins:** `SyncedStructure` is assembled from one mixin per area. Every mixin inherits `StructureState` (`base.py`), which declares the shared attributes and, for type checkers only, the methods mixins call on each other. If a mixin starts calling a method of another mixin, add its signature there.
- **Tables first:** analyses return DataFrames (or dicts of them) with both chain IDs and 1-based indices; figures are drawn from those tables.
- **Units:** nm for contacts, Å for interface metrics and docking.
- **Names:** a name says what it returns. A verb prefix is used only when it fixes the return type or a side effect (`set_` → `self`, `export_` → paths, `plot_` → figure); otherwise the name is a noun whose last word is the type (`contact_table`, `chain_partner_map`, `contact_string`). Every `_table` is a method. New names must not be longer than needed. [Conventions](guide/conventions.md) states the full grammar; renames replace names outright, without deprecated aliases, and are listed in the changelog.
- **Columns:** tables that list chains carry `chain_key`, `chain_id`, `auth_chain_id` (`chain_key_a`, `chain_a`, `auth_chain_a` and `_b` for pairs), so any value read from a table can be passed back to a method.

## Tests

```bash
pytest                 # offline tests (default): library, docs, type and docstring checks
pytest -m network      # live RCSB services
pytest -m haddock      # real haddock3 runs at low throughput (a few minutes)
pytest -m ""           # everything — run this before every commit
```

| Layer | File | What it checks |
|---|---|---|
| Smoke and regressions | `tests/test_smoke.py`, `tests/test_regressions.py` | Every analysis on 1AAY, edge cases, bugs that were fixed |
| Docking | `tests/test_docking.py` | HADDOCK3 set up offline; real low-throughput runs with `-m haddock` |
| Network | `tests/test_network.py` | Live RCSB download, search and citations |
| Documentation | `tests/test_docs.py` | README and guide code runs; links resolve; the MkDocs site builds with `--strict`; the API reference lists every public name once; the configuration page lists every setting |
| Types and docstrings | `tests/test_typing.py` | See below |

`tests/test_typing.py` is a second layer for the type hints; each check is skipped when its tool is not installed:

| Check | Tool | What it verifies |
|---|---|---|
| Coverage | standard library (`ast`) | Every function has parameter and return types and a docstring |
| Style | [ruff](https://docs.astral.sh/ruff/) (`D`, `ANN` rules) | Google-style docstrings, no missing annotations |
| Docstrings vs signatures | [pydoclint](https://jsh9.github.io/pydoclint/) | Every argument documented, no stale ones, `Returns:`/`Raises:` present |
| Static types | [mypy](https://mypy.readthedocs.io/) | Annotations consistent with the code |
| Runtime types | [typeguard](https://typeguard.readthedocs.io/) | Re-runs the offline suite checking every annotated call and return |

Runtime types across the whole suite, network and haddock3 included:

```bash
pytest -m "" --typeguard-packages=synced_structure
```

The tools can also be run directly: `ruff check src`, `pydoclint src`, `mypy`. Their settings are in `pyproject.toml`. mypy targets Python 3.12 because numpy's type stubs need it; the package itself supports Python 3.10 and newer.

## Documentation

The site is built with [MkDocs](https://www.mkdocs.org/), the [Material](https://squidfunk.github.io/mkdocs-material/) theme and [mkdocstrings](https://mkdocstrings.github.io/python/), configured in `mkdocs.yml`. It is built locally and not published.

```bash
mkdocs serve                       # live preview at http://127.0.0.1:8000
mkdocs build --strict              # what the tests run; fails on any warning
```

- `README.md` holds the pitch, a short install and the quick start. Pages under `docs/` hold the rest: `quickstart.md`, `guide/*.md` (concepts, conventions, configuration, docking) and `reference/*.md` (API).
- The API reference is generated from signatures and docstrings. Each `reference/*.md` page lists the `SyncedStructure` members of one area under `members:`. A new public method must be added to one of these lists; `tests/test_docs.py` fails until it is.
- Python blocks in the README and in `docs/index.md`, `docs/quickstart.md` and `docs/guide/*.md` run in `tests/test_docs.py` (one namespace per page, with `from_rcsb` served from `tests/data`). Mark a block that must not run with `<!-- skip-test -->` on the line before it.

## Releasing

1. `pytest -m ""` passes, and `python examples/export_all.py` and `python examples/dock_1kxp.py` run.
2. Update `CHANGELOG.md` and the version in `pyproject.toml` and `src/synced_structure/__init__.py`.
3. Commit, then tag: `git tag v0.1.0 && git push --tags`.
