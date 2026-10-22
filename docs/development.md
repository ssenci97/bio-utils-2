# Development

How the code is written and checked. Read this before changing the package.

## Set up

```bash
git clone https://github.com/ssenci97/bio-utils-2.git
cd bio-utils-2
pip install -e ".[dev]"            # editable install + pytest, ruff, pydoclint, mypy, typeguard, type stubs
pip install -e ".[dev,haddock]"    # also haddock3, for the docking tests
```

Always use the editable install (`-e`): tests and scripts then import the code in `src/`, not an old copy in `site-packages`.

## Conventions

- **Docstrings:** [Google style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) on every function and method: a summary line, then `Args:`, `Returns:` and `Raises:` where they apply. Keep them short; types belong in the signature, not the docstring.
- **Type hints:** every parameter and return is annotated (not `self`/`cls`). Methods that change the structure return `Self`. Shared aliases live in `src/synced_structure/_types.py`.
- **Mixins:** `SyncedStructure` is assembled from one mixin per area. Every mixin inherits `StructureState` (`base.py`), which declares the shared attributes and, for type checkers only, the methods mixins call on each other. If a mixin starts calling a method of another mixin, add its signature there.
- **Tables first:** analyses return DataFrames (or dicts of them) with both chain IDs and 1-based indices; figures are drawn from those tables.
- **Units:** nm for contacts, Å for interface metrics and docking.
- **Names:** a name says its domain and what it returns (`interface_buried_area`, `docked_model_table`). The [API index](api_index.md) lists them all; renames replace names outright, without deprecated aliases.

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
| Documentation | `tests/test_docs.py` | Every guide's code runs; every relative link resolves |
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

- The README is the landing page; guides go in `docs/guides/`, definitions in `docs/concepts.md`, the reference in `docs/api_index.md`.
- Code in the guides runs in `tests/test_docs.py`, top to bottom per guide, with `from_rcsb` served from `tests/data`. Mark a block that must not run (a real haddock3 run, a network search) with `<!-- skip-test -->` on the line before it.
- When you add or rename a public name, update `docs/api_index.md`.

## Releasing

1. `pytest -m ""` passes, and `python examples/export_all.py` runs.
2. Update `CHANGELOG.md` and the version in `pyproject.toml` and `src/synced_structure/__init__.py`.
3. Commit, then tag: `git tag v0.1.0 && git push --tags`.
