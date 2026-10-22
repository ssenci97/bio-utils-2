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
- Tests: offline, network, low-throughput haddock3, documentation (README code runs, API page is current) and typing layers.
- Documentation: a short [README](README.md) with installation and quick start, and a single-page HTML API reference generated from the docstrings by `docs/build_docs.py`.

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
