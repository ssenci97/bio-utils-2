# Changelog

## 0.1.0

First stable version.

### Structures and numbering
- `SyncedStructure` wraps a Biopython structure and keeps, for every chain, the observed → reference mapping, both chain IDs (label and author) and the entity, rebuilding the mapping whenever a chain, reference or numbering changes.
- Loaders: `from_rcsb`, `from_path` (PDB and mmCIF), `from_biopython`; chains keyed by author (default) or label IDs.
- Chain management: reference sequences, renumbering, renaming (by mapping or by sequence), subsetting, sequence search.

### Analyses
- Contacts with MDTraj (`ca` and `sidechain` schemes), contact strings on observed and reference sequences, pair and block matrices, contact statistics.
- Interface metrics: buried surface, interface residues and partners, contact scores.
- Contact maps (pair and multi-chain), lettered panel grids, 3Dmol.js views; TSV, PDF and HTML export.
- RCSB search and citation lookup.

### Docking
- `HaddockDocking` / `prepare_haddock_docking`: HADDOCK3 set up from your chains and numbering (ambiguous, unambiguous, patch and contact restraints), `ab-initio`, `oracle` and `eval` presets, runs in a subprocess, scores and docked models read back as `SyncedStructure` objects.

### Quality
- Google-style docstrings and type hints on every function, checked by ruff, pydoclint and mypy, and at runtime by typeguard.
- Tests: offline, network, low-throughput haddock3, documentation (every guide is executed) and typing layers.
- Documentation: [concepts](docs/concepts.md), six step-by-step [guides](docs/README.md) and an [API index](docs/api_index.md).

### Changed names
Several methods were renamed to say what they do; the old names were removed. See [Renamed names](docs/api_index.md#renamed-names) to update existing scripts.

### Fixed
- `rename_chains_by_sequence` and `find_chains_by_sequence` no longer accept a match made of a few aligned letters (for example a DNA strand matching a protein); both take `min_coverage` (default 50 %).
