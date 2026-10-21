# synced_structure

Synchronized Biopython structures: chain mapping, contacts, contact maps, interfaces and 3D views.

`synced_structure` wraps a `Bio.PDB.Structure` and links every observed residue to its **reference sequence** (the full mmCIF entity sequence, or one you supply). That observed→reference mapping stays consistent through renaming, replacing, renumbering and re-referencing chains, so contacts, interfaces, tables and plots can always be reported in reference coordinates.

See [PHILOSOPHY.md](PHILOSOPHY.md) for the design principles and the full API reference.

## Installation

Requires Python ≥ 3.10.

Install directly from GitHub:

```bash
pip install "git+https://github.com/ssenci97/bio-utils-2.git"
pip install "git+https://github.com/ssenci97/bio-utils-2.git@main"   # pin a branch, tag or commit
```

Or clone for development:

```bash
git clone https://github.com/ssenci97/bio-utils-2.git
cd bio-utils-2
pip install -e .            # library only
pip install -e ".[dev]"     # library + pytest
```

`mdtraj` ships compiled extensions; if `pip` cannot build it on your platform, install it from conda-forge first:

```bash
conda install -c conda-forge mdtraj
pip install -e .
```

### Dependencies

| Package | Used for |
|---|---|
| `biopython>=1.80` | Parsing PDB/mmCIF, structure objects, sequence alignment, Shrake–Rupley SASA |
| `numpy`, `pandas` | Numerics and the tabular outputs every analysis returns |
| `mdtraj` | Fast residue–residue distance and contact calculations |
| `matplotlib` | Contact maps, multi-chain contact plots, panel grids, PDF export |
| `pyyaml` | Optional YAML configuration |
| `requests` | Downloading entries and querying the RCSB Search and GraphQL APIs |
| `tqdm` | Progress bars |
| `ipython` | Inline 3D views in Jupyter |
| `pytest` *(dev)* | Test suite |

The 3D viewer loads [3Dmol.js](https://3dmol.org) from its CDN at display time, so it needs a browser with internet access; no Python package is required for it.

## Capabilities

- **Loading** from RCSB (`SyncedStructure.from_rcsb("4HHB")`), a local `.pdb`/`.cif`/`.mmcif` file (`from_path`), or an existing Biopython structure (`from_biopython`). Chains can be keyed by author (`key_by="auth"`, default) or label (`key_by="label"`) IDs; tables always report both.
- **Chain and reference management**: per-chain observed and reference sequences, external accessions (e.g. UniProt), custom references, adding/replacing/removing chains, renaming chains by mapping or by sequence identity, and finding chains that match a sequence. Every change re-syncs the mapping.
- **Alignment and residue mapping**: global alignment of observed to reference sequence, per-chain alignment statistics and column-by-column positions, residue tables (PDB number ↔ sequence index ↔ reference index), moiety classification (protein, nucleic, ligand, ion, water), and renumbering by reference or 1..N.
- **Contacts**: intra-/inter-chain contacts via MDTraj with Cα or closest-heavy-atom schemes, contact matrices and block matrices, contact statistics, and contact sequences exportable in observed or reference numbering.
- **Interfaces**: buried surface area, SASA, buried fraction, contact count, minimum distance, contact score, interface residues and per-chain partner summaries.
- **Visualization**: contact maps, multi-chain contact plots, combined contact panels, multi-panel figure grids exported to PDF, and interactive 3Dmol.js views in Jupyter (`show_structure_3d`, `side_by_side`).
- **Export**: CSV export of chains, residue mappings, moieties, alignment statistics/positions and interface tables; figure export in several formats.
- **RCSB search**: `search_rcsb` queries the RCSB Search API by organism, keyword, method, resolution, release date, UniProt ID, EC number and more; `search_rcsb_pubs` attaches publication metadata to the results.

## Quick start

```python
from synced_structure import SyncedStructure, search_rcsb

s = SyncedStructure.from_rcsb("4HHB")
s.summary()                                  # chains, entities, references, coverage
s.get_residue_mapping("A")                   # observed ↔ reference numbering
s.contact_matrix("A", "B")                   # inter-chain contact map data
s.interface_complex("A", "B")                # interface metrics
s.plot_contact_map("A", "B")
s.show_structure_3d(style="cartoon")              # 3D view in Jupyter

hits = search_rcsb(organisms=["Homo sapiens"], keywords=["hemoglobin"], resolution_max=2.0)
```

A full worked example, the STAT6-DNA complex 5D39, is in [`examples/stat6_dna_5d39.ipynb`](examples/stat6_dna_5d39.ipynb).

## Configuration

All constructors accept an optional `config_path` pointing to a YAML file. Values are read from the `synced_structure` section; missing keys fall back to defaults, and the file is never written.

```yaml
synced_structure:
  save_outdir: data/saved          # default output directory for exports
  aln_match_score: 1
  aln_mismatch_score: -1
  aln_open_gap_score: -2
  aln_extend_gap_score: -0.5
  interface_contact_cutoff_a: 4.5  # Å
  interface_score_cutoff_a: 12.0   # Å
  interface_block_size: 1000
```

## Tests

```bash
pytest                 # offline tests (default)
pytest -m network      # tests that call live RCSB services
```
