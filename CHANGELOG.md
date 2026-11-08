# Changelog

## Unreleased

### Added
- `core_count()`: the CPU cores this process may use, the smallest of the affinity mask, the cgroup CPU quota
  (containers, Kubernetes, systemd, Slurm) and the batch-scheduler allocation (`SLURM_CPUS_PER_TASK`,
  `SLURM_CPUS_ON_NODE`, PBS `NCPUS`, LSF `LSB_DJOB_NUMPROC`, SGE `NSLOTS`). `os.cpu_count()` counts the whole node.
- Operations map (`docs/assets/operations.png`, at most 1920 x 1080 px) on the API reference page: the parts of a
  `SyncedStructure` (and of each `ChainRecord`), with every public method linked to the parts it changes or reads, as
  a Graphviz directed graph from `docs/scripts/operations_graph.py`. An mkdocs hook (`docs/hooks.py`) redraws
  it before each build. `tests/test_docs.py` fails when a public method is missing from it. `graphviz` joins the `docs`
  extra.
- `SyncedStructure.export_structure(file_format="cif" | "pdb", chain_ids=None)` writes the structure, or some chains,
  to a file; `mmcif_string()` and `pdb_string()` return the text. mmCIF keeps label and author chain IDs, entities,
  descriptions, accessions, and reference sequences, so `from_path()` reloads the same chains. PDB keeps author chain
  IDs and raises a clear error for chain IDs longer than one character, more than 99,999 atoms, or residue numbers
  outside -999..9999.
- `method="alascan"` for `HaddockDocking`: no rigid-body docking. The partners are written as one complex, refined
  (`mdref`), clustered (`rmsdmatrix`, `clustrmsd`), alanine-scanned (`alascan`) and evaluated (`caprieval`).
  `add_scan_residues()` restricts the scan in your numbering. `alascan_tables()` returns the per-model and per-cluster
  results mapped back to chain keys and reference positions.
- HADDOCK3 partners from different structures: `HaddockDocking.set_reference(complex, chain_map=...)` names a complex
  (for example the holo complex of two apo partners) that caprieval compares models with, written to
  `inputs/reference.pdb` in the partners' segment IDs and HADDOCK3 numbering, and that `add_contact_restraints()` reads
  contacts from by default, so `oracle` works on mixed partners. `reference=` also accepts a `SyncedStructure`.
- `HaddockDocking.plot_caprieval(x="irmsd", y="score", color="dockq")`: a scatter plot of caprieval metrics, returning
  `(fig, ax)`.
- `examples/docking/barnase_barstar_oracle.py`: docking of unbound barnase (1A2P) and barstar (1A19), evaluated against
  the complex (1BRS). `examples/docking/1kxp_oracle.py` also saves the caprieval plot.
- `run_haddock` logs progress: module starts and finishes, models written so far with the time left, and stops haddock3 when Python is interrupted. `progress=False` turns it off.

### Changed: core detection
- `HaddockDocking(ncores=None)` is the new default and uses `core_count()`; it was a fixed 4. The examples use
  `core_count()` instead of `min(4, os.cpu_count())`.
- Docking examples moved to `examples/docking/`: `1kxp_oracle.py`, `1kxp_ab_initio.py`, `barnase_barstar_oracle.py`
  and `barnase_barstar_ab_initio.py` (new), with every argument of the docking calls written out. Outputs go to
  `examples/output/docking_<system>_<method>/`.

### Changed: 3D views
- `show_structure_3d` and `export_structure_3d` build their mmCIF with the same writer as `mmcif_string()`, with chain
  keys as author IDs so colors and click labels still follow your chain keys.

### Changed: HADDOCK3 partners
- Residues transferred between structures (contacts from `source=`, the reference complex) are matched through a
  global alignment of the two reference sequences when they differ (tags, truncations, mutants), instead of assuming
  equal reference positions. A log line gives the identity of each aligned pair.
- `separate="auto"` also moves partners apart when they come from different structures but touch (heavy atoms within
  `clearance_a`), as superposed models or models centred at the origin do.

### Changed: HADDOCK3 presets
- Every preset clusters the rigid-body models and refines only the best clusters: `clustfcc` → `seletopclusts` (the top fifth of the clusters, at least 5, 4 models each) → `emref` → `caprieval`. Before, `oracle` and `eval` refined all 200 rigid-body models and `ab-initio` refined the top 200 before clustering.
- The `eval` method is renamed `targeted`: docking from patch-pair restraints you give. It no longer requires a reference; like every method, it ends with `caprieval` against the reference when there is one, otherwise against the best-scoring model.

### Documentation
- A documentation site built with MkDocs, Material and mkdocstrings: installation, a quick start, guides (concepts, conventions, configuration, HADDOCK3 docking) and an API reference organised by pipeline step. It replaces `docs/index.html` and `docs/build_docs.py`.
- `examples/docking/1kxp_oracle.py`: a HADDOCK3 re-docking of actin and vitamin D-binding protein (1KXP), compared with the crystal.
- HADDOCK3 tutorials (`docs/tutorials/haddock/`): re-docking a complex, restraints, partners from different
  structures, workflow and parameters, and analysis of the results. Their code blocks run in the tests.
- The guide code blocks run in the tests, like the README's, and `mkdocs build --strict` is part of the test suite.

### Fixed
- A chain's mapping was rebuilt by global alignment whenever any chain was added, replaced, renamed or renumbered. This silently discarded mappings set by `apply_reference_numbering` and empty mappings from `set_reference_sequence(remap=False)`. Mappings are now kept while the chain's observed and reference sequences are unchanged.

### Changed: methods
Each change makes a method follow the naming rules: a verb fixes the return type, and every `_table` is a method.

| Before | Now |
|---|---|
| `s.residue_table`, `s.moiety_table`, `s.alignment_table`, `s.alignment_column_table` (properties) | `s.residue_table()`, `s.moiety_table()`, `s.alignment_table()`, `s.alignment_column_table()` |
| `dock.residue_table`, `dock.restraint_table`, `dock.contact_patch_table` (properties) | `dock.residue_table()`, `dock.restraint_table()`, `dock.contact_patch_table()` |
| `rename_chains_by_sequence(...)` returned a report DataFrame | returns the structure, like `rename_chains`; unmatched targets are logged; inspect matches first with `find_chains_by_sequence` |
| `renumber_chain(chain, in_place=False)` returned a renumbered Biopython `Chain` copy | `renumber_chain(chain)` renumbers in place and returns the structure; for a copy, use `subset_chains([chain]).renumber_chain(chain)` |
| `export_contact_strings(output_path, ...)` | `export_contact_strings(outdir=None, ...)`, writing `{structure_id}_contact_strings.tsv` (`_reference.tsv` with `reference=True`) |
| `export_chain_table` wrote a reduced table (entity, description, residue count) | writes `chain_table()` |
| — | `export_contact_stat_tables(outdir=None, ...)` writes the four `contact_stat_tables` |
| `export_plot_grid(..., name)` wrote `name` | writes `{structure_id}_{name}`, like every other export |
| `chain_ids_by_moiety` returned label IDs | returns chain keys |
| `set_accession(chain_id, external_ref)`, `set_reference_id(chain_id, refseq_id)` | `set_accession(chain_id, accession)`, `set_reference_id(chain_id, reference_id)` |
| `plot_contact_map(..., save_path=)`, `plot_multi_contact_map(..., save_path=)` | removed: save with `export_plot(lambda: s.plot_contact_map(...), name)`, which applies the PDF settings |
| default chains `"A"`/`"B"` in `contact_matrix`, `block_contact_matrix`, `interface_*`, `plot_contact_map`, `plot_multi_contact_map` | chain arguments are required (`chain_b=None` still means intra-chain) |
| `chain_partner_structures(moieties=POLYMER_MOIETIES)` | same value, now the shared `DEFAULT_CONTACT_MOIETIES` default |

### Changed: columns
Every table that lists chains now has a `chain_key` column, or `chain_key_a`/`chain_key_b` for pairs, so values read from a table can be passed back to methods (label and author IDs can be different letters for the same chain). Column names are shared across tables:

| Table | Before | Now |
|---|---|---|
| `contact_table` | `chain_1`, `auth_chain_1`, `resnum_1`, `resname_1`, `sequence_index_1`, `moiety_1` (and `_2`) | `chain_key_a`, `chain_a`, `auth_chain_a`, `resnum_a`, `resname_a`, `sequence_index_a`, `moiety_a` (and `_b`) |
| interface complex rows (`interface_complex`, `interface_tables()["complex"]`, ...) | `structure` | `structure_id`, plus `chain_key_a`, `chain_key_b` |
| interface residue rows | `structure`, `chain`, `auth_chain` | `structure_id`, `chain_key`, `chain_id`, `auth_chain_id` |
| `contact_stat_tables()["pairs"]` | `chain_a`, `chain_b` (label IDs) | `chain_key_a`, `chain_a`, `auth_chain_a`, `chain_key_b`, `chain_b`, `auth_chain_b`; rows in chain order |
| `contact_stat_tables()["profile"]`, `["separation"]` | `chain_id` (separation), `chain_id`, `auth_chain_id` (profile) | `chain_key`, `chain_id`, `auth_chain_id` |
| `chain_table`, `residue_table`, `moiety_table`, `chain_residue_table`, `alignment_table`, `alignment_column_table`, `find_chains_by_sequence`, contact strings export | — | `chain_key` added |

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
- Tests: offline, network, low-throughput haddock3, documentation (README code runs, API page is current) and typing layers.
- Documentation: a short README with installation and quick start, and a single-page HTML API reference generated from the docstrings (replaced by the MkDocs site in the next version).

### Changed names
Methods were renamed so that each name says what it returns: a verb prefix only where it fixes the return type or a side effect, otherwise a noun ending in the type (`contact_table`, `chain_partner_map`, `interface_tables`, `step_paths`). The six scalar interface methods (`interface_buried_area` and the others) were removed; read the value from `interface_complex(a, b)`. The two interface exports were merged into `export_interface_tables`. The old names were removed; update existing scripts with this table.

| Old | New |
|---|---|
| `subset`, `match_chains`, `rename_by_sequence`, `summary` | `subset_chains`, `find_chains_by_sequence`, `rename_chains_by_sequence`, `chain_table` |
| `residue_map`, `moiety_summary`, `reindex_chain` | `residue_table`, `moiety_table`, `renumber_chain` |
| `alignment_stats`, `alignment_positions` | `alignment_table`, `alignment_column_table` |
| `to_ca_trajectory`, `get_contact_sequence`, `get_reference_contact_sequence`, `export_contact_sequences` | `contact_point_trajectory`, `contact_string`, `reference_contact_string`, `export_contact_strings` |
| `plot_contacts_multi` | `plot_multi_contact_map` |
| `interface`, `interface_all`, `chain_partner_frame` | `interface_metrics`, `interface_tables`, `chain_partner_table` |
| `buried_area`, `initial_sasa`, `buried_fraction`, `contact_count`, `min_distance`, `contact_score` and their `interface_*` versions | removed: read `interface_complex(a, b)["dASA_A2"]`, `["initial_sasa_A2"]`, `["frac_dASA"]`, `["n_contacts"]`, `["min_dist_A"]`, `["contact_score"]` |
| `export_chains`, `export_residues`, `export_moieties`, `export_alignment_stats`, `export_alignment_positions` | `export_chain_table`, `export_residue_table`, `export_moiety_table`, `export_alignment_table`, `export_alignment_column_table` |
| `search_rcsb_pubs`, `side_by_side` | `fetch_rcsb_citations`, `show_side_by_side_3d` |
| `jupy_viz`, `complete_contacts`, `partners`, `from_rcsb_id`, `residue_moiety` | `show_structure_3d`, `block_contact_matrix`, `interface_partner_map`, `from_rcsb`, `residue_class` |
| `chain_summary`, `chain_contact_partners`, `chain_partner_interfaces`, `chains_by_moiety` | `chain_table`, `chain_partner_map`, `chain_partner_metrics`, `chain_ids_by_moiety` |
| `get_external_ref` / `set_external_ref`, `get_refseq_id` / `set_refseq_id`, `get_residue_mapping` | `get_accession` / `set_accession`, `get_reference_id` / `set_reference_id`, `chain_residue_table` |
| `reference_alignment_stats`, `reference_alignment_positions` | `alignment_table`, `alignment_column_table` |
| `compute_contacts`, `contact_statistics` | `contact_table`, `contact_stat_tables` |
| `interface_metrics_all`, `interface_residues`, `interface_partner_labels` | `interface_tables`, `interface_residue_table`, `interface_partner_map` |
| `export_chain_summary`, `export_residue_mapping`, `export_reference_alignment_stats`, `export_reference_alignment_positions` | `export_chain_table`, `export_chain_residue_table`, `export_alignment_table`, `export_alignment_column_table` |
| `export_interface_complex` + `export_interface_residues` | `export_interface_tables` (writes both files) |
| `HaddockDocking`: `haddock_residue_map`, `run_steps`, `haddock_scores`, `partner_contact_pairs`, `restraint_file_text`, `write_inputs`, `run_caprieval_on_steps` | `residue_table`, `step_paths`, `score_table`, `partner_contact_map`, `restraint_text`, `export_inputs`, `run_caprieval` |
| `record["field"]` on a chain record | `record.field` |
| `HaddockSetup`, `prepare_haddock` and its methods | `HaddockDocking`, `prepare_haddock_docking` |

### Fixed
- `rename_chains_by_sequence` and `find_chains_by_sequence` no longer accept a match made of a few aligned letters (for example a DNA strand matching a protein); both take `min_coverage` (default 50 %).
