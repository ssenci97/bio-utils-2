# synced_structure

A Biopython structure wrapper that keeps every observed residue linked to its position in the reference sequence, and builds contacts, contact maps, interface metrics, HADDOCK3 dockings, tables and 3D views on top of that link.

**Documentation:** [`docs/index.html`](docs/index.html) is the API reference: every function and method with its signature, plus conventions, concepts and configuration. · [Changelog](CHANGELOG.md) · [Development](docs/development.md)

## Installation

### 1. Requirements
- Python **3.10 or newer**.
- A C compiler is only needed if `pip` has to build `mdtraj` from source; the conda route below avoids it.

### 2. Create an environment (recommended)

```bash
conda create -n synced python=3.12
conda activate synced
conda install -c conda-forge mdtraj      # prebuilt binary; skip if pip can install it on your system
```

A plain virtual environment works too: `python -m venv .venv && source .venv/bin/activate`.

### 3. Install the package

```bash
pip install "git+https://github.com/ssenci97/bio-utils-2.git"              # latest version
pip install "git+https://github.com/ssenci97/bio-utils-2.git@v0.1.0"       # a fixed release
```

This installs `biopython`, `numpy`, `pandas`, `mdtraj`, `matplotlib`, `pyyaml`, `requests`, `tqdm` and `ipython`.

### 4. Optional extras

| Need | Install |
|---|---|
| HADDOCK3 docking | `pip install "synced_structure[haddock] @ git+https://github.com/ssenci97/bio-utils-2.git"`, then check with `haddock3 -h` |
| Interactive 3D views | `pip install jupyterlab`; views load [3Dmol.js](https://3dmol.org) from the web, so the browser needs internet access |
| Development (tests, linters) | `git clone https://github.com/ssenci97/bio-utils-2.git && cd bio-utils-2 && pip install -e ".[dev]"` |

### 5. Check the installation

```bash
python -c "import synced_structure; print(synced_structure.__version__)"
```

## Quick start

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")          # download Zif268 zinc fingers bound to DNA
# s = SyncedStructure.from_path("model.cif")   # or load a local .pdb / .cif / .mmcif file

chains = s.chain_table()                       # one row per chain: structure and chain IDs, entity, references, coverage
residues = s.chain_residue_table("A")          # author number ↔ observed index ↔ reference index
contacts = s.contact_table(cutoff=0.8)         # inter-chain residue pairs closer than 0.8 nm
metrics = s.interface_complex("A", "B")        # buried area (Å²) and other interface metrics
print(chains[["chain_id", "auth_chain_id", "entity_id"]])
print(len(contacts), "contacts;", round(metrics["dASA_A2"]), "Å² buried between A and B")

fig, ax = s.plot_contact_map("A", "B", crop=True)        # contact map of one chain pair
paths = s.export_chain_table(outdir="results")           # results/1AAY_chains.tsv
view = s.show_structure_3d()                             # interactive 3D view (display it in Jupyter)
```

Every analysis returns a `pandas.DataFrame` (or a dict of them). Contact distances are in **nm**; interface metrics and docking are in **Å**. Names say what they return: `_table` is a DataFrame, `_map` a dict, `export_` writes files and returns their paths. The full list is in the [API reference](docs/index.html).

A complete script that writes every table, contact map and a 3D page is in [`examples/export_all.py`](examples/export_all.py):

```bash
python examples/export_all.py      # writes to examples/output/
```

## Rebuilding the documentation

```bash
python docs/build_docs.py          # regenerates docs/index.html from the signatures and docstrings
```

To publish it online, enable GitHub Pages for the `main` branch, folder `/docs`.
