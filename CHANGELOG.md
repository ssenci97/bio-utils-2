# Changelog

## Unreleased

### Documentation
- A documentation site built with MkDocs, Material and mkdocstrings and published to GitHub Pages: installation, a quick start, guides (concepts, conventions, configuration, HADDOCK3 docking) and an API reference organised by pipeline step. It replaces `docs/index.html` and `docs/build_docs.py`.
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
