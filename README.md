# synced_structure

A Biopython structure wrapper that keeps every observed residue linked to its position in the reference sequence, and builds contacts, contact maps, interface metrics, HADDOCK3 dockings, tables and 3D views on top of that link.

## 📖 Documentation: [`docs/`](docs/)

[Quick start](docs/quickstart.md) · [Concepts](docs/guide/concepts.md) · [API reference](docs/reference/index.md) · [HADDOCK3 docking](docs/guide/docking.md) · [Changelog](CHANGELOG.md) · [Development](docs/development.md)

The site can be built and previewed locally with `pip install -e ".[docs]" && mkdocs serve`.

## Installation

Python 3.10 or newer. Install `mdtraj` from conda-forge first if `pip` cannot build it on your system:

```bash
conda create -n synced python=3.12 && conda activate synced
conda install -c conda-forge mdtraj
pip install "git+https://github.com/ssenci97/bio-utils-2.git"
```

Extras: `[haddock]` for HADDOCK3 docking, `[docs]` to build the documentation, `[dev]` for tests and linters. The [installation page](docs/install.md) has the details.

## Quick start

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")          # download Zif268 zinc fingers bound to DNA
# s = SyncedStructure.from_path("model.cif")   # or load a local .pdb / .cif / .mmcif file

chains = s.chain_table()                       # one row per chain: chain key, label and author IDs, entity, coverage
residues = s.chain_residue_table("A")          # author number ↔ observed index ↔ reference index
contacts = s.contact_table(cutoff=0.8)         # inter-chain residue pairs closer than 0.8 nm
metrics = s.interface_complex("A", "B")        # buried area (Å²) and other interface metrics
print(chains[["chain_key", "chain_id", "auth_chain_id", "entity_id"]])
print(len(contacts), "contacts;", round(metrics["dASA_A2"]), "Å² buried between A and B")

fig, ax = s.plot_contact_map("A", "B", crop=True)        # contact map of one chain pair
path = s.export_chain_table(outdir="results")            # results/1AAY_chains.tsv
view = s.show_structure_3d()                             # interactive 3D view (display it in Jupyter)
```

Every analysis returns a `pandas.DataFrame` (or a dict of them). Contact distances are in **nm**; interface metrics and docking are in **Å**. Names say what they return: `_table` is a DataFrame, `_map` a dict, `export_` writes files and returns their paths.

A complete script that writes every table, contact map and a 3D page is in [`examples/export_all.py`](examples/export_all.py):

```bash
python examples/export_all.py      # writes to examples/output/
```

[`examples/dock_1kxp.py`](examples/dock_1kxp.py) re-docks actin and vitamin D-binding protein (1KXP) with HADDOCK3 and compares the models with the crystal (needs the `haddock` extra; writes to `examples/output/docking_1KXP/`). [`examples/dock_barnase_barstar.py`](examples/dock_barnase_barstar.py) docks two unbound structures from different PDB entries and evaluates against their complex. The [HADDOCK3 tutorials](docs/tutorials/haddock/index.md) walk through docking step by step.
