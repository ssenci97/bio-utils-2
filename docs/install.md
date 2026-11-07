# Installation

## Requirements

- Python **3.10 or newer**.
- A C compiler only if `pip` has to build `mdtraj` from source; the conda route avoids it.

## Create an environment (recommended)

```bash
conda create -n synced python=3.12
conda activate synced
conda install -c conda-forge mdtraj      # prebuilt binary; skip if pip can install it on your system
```

A plain virtual environment works too: `python -m venv .venv && source .venv/bin/activate`.

## Install the package

```bash
pip install "git+https://github.com/ssenci97/bio-utils-2.git"              # latest version
pip install "git+https://github.com/ssenci97/bio-utils-2.git@v0.1.0"       # a fixed release
```

This installs `biopython`, `numpy`, `pandas`, `mdtraj`, `matplotlib`, `pyyaml`, `requests`, `tqdm` and `ipython`.

## Optional extras

| Need | Install |
|---|---|
| HADDOCK3 docking | `pip install "synced_structure[haddock] @ git+https://github.com/ssenci97/bio-utils-2.git"`, then check with `haddock3 -h` |
| Interactive 3D views | `pip install jupyterlab`; views load [3Dmol.js](https://3dmol.org) from the web, so the browser needs internet access |
| Building these docs | `pip install -e ".[docs]"`, then `mkdocs serve` |
| Development (tests, linters) | `git clone https://github.com/ssenci97/bio-utils-2.git && cd bio-utils-2 && pip install -e ".[dev]"` |

## Check the installation

```bash
python -c "import synced_structure; print(synced_structure.__version__)"
```
