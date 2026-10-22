# API index

This page lists every public name in `synced_structure`, grouped by area, with a one-line description. Every name also has a full Google-style docstring: `help(SyncedStructure.compute_contacts)` shows it. Names changed in this version are listed in [Renamed names](#renamed-names); the old names no longer exist. For definitions, see [Concepts](concepts.md); for step-by-step use, the [guides](README.md).

Conventions:
- Indices and positions are 1-based in every table.
- Distances are in **nm** for contacts and in **Å** for interface metrics and docking (names say `_a` where it matters).
- A chain can be named by its key, label ID or author ID.
- Methods that change the structure return `self`.

## `SyncedStructure`

### Loading
| Name | Kind | Description |
|---|---|---|
| `SyncedStructure.from_rcsb(pdb_id, config_path=None, key_by="auth")` | classmethod | Download an mmCIF entry from RCSB |
| `SyncedStructure.from_path(path, config_path=None, key_by="auth")` | classmethod | Load a local `.pdb`, `.cif` or `.mmcif` file |
| `SyncedStructure.from_biopython(structure, ...)` | classmethod | Wrap an existing Biopython structure |
| `SyncedStructure(structure, structure_id="structure", mmcif_dict=None, file_format="cif", config_path=None, default_entity_id=None, key_by="auth")` | constructor | Build from a Biopython structure directly (the loaders call this) |
| `structure`, `model`, `structure_id`, `file_format`, `mmcif_dict`, `chains`, `key_by`, `config` | attributes | Underlying Biopython objects, ID, source format, mmCIF data, chain records, key mode, config |

### Chains and references
| Name | Kind | Description |
|---|---|---|
| `chain_ids` | property | Loaded chain keys |
| `chain_summary(partners=False, ...)` | table | One row per chain: IDs, entity, moiety, description, references, mapping coverage |
| `get_sequence(chain)` / `get_reference_sequence(chain)` | str | Observed / reference one-letter sequence |
| `set_reference_sequence(chain, seq, remap=True)` / `set_reference_sequences({chain: seq})` | mutator | Replace references and rebuild the mapping |
| `remap_chain(chain)` / `remap_all_chains()` | mutator | Rebuild observed→reference mappings |
| `get_external_ref` / `set_external_ref` | get / set | External accession (e.g. UniProt) |
| `get_refseq_id` / `set_refseq_id` | get / set | ID of the reference sequence |
| `set_entity_id`, `set_chain_metadata` | mutator | Annotate a chain |
| `add_chain`, `replace_chain`, `remove_chain` | mutator | Edit the chain set |
| `subset_chains(chains, structure_id=None)` | new structure | Copy holding only some chains |
| `find_chains_by_sequence(sequence, min_identity=90, min_coverage=50)` | table | Chains whose sequence matches a query |
| `rename_chains({old: new})` | mutator | Rename chain keys, keeping all maps in sync |
| `rename_chains_by_sequence({new: sequence}, min_identity=90, min_coverage=50)` | table + mutator | Rename chains matched by sequence |
| `chain_contact_partners(cutoff, moieties)` | dict | `{chain: [chains in contact]}` |

### Residues and reference alignment
| Name | Kind | Description |
|---|---|---|
| `get_residue_mapping(chain)` | table | Author number ↔ sequence index ↔ reference index for one chain |
| `residue_table` | property, table | Every polymer residue of every chain with its indices |
| `moiety_table` | property, table | Residue and atom counts per chain and moiety |
| `chains_by_moiety(moiety)` | list | Chains with `protein`, `nucleic` or `other` residues |
| `apply_reference_numbering(chains=None)` | mutator | Trim unmapped residues, renumber by reference position |
| `renumber_chain(chain, in_place=False)` | chain | Renumber polymer residues 1..N |
| `reference_alignment_stats` | property, table | Local alignment statistics of each chain against its reference |
| `reference_alignment_positions` | property, table | The same alignment column by column |

### Contacts
| Name | Kind | Description |
|---|---|---|
| `compute_contacts(cutoff, contact_type, chain_ids, scheme, moieties)` | table | Residue pairs closer than `cutoff` (nm) |
| `contact_point_trajectory(chain_ids, moieties)` | MDTraj | One contact point per residue: Cα, C4′ or centroid |
| `contact_string(chain, ...)` | str | Per residue `1` contact / `0` none / `X` no contact point |
| `reference_contact_string(chain, ...)` | str | The same, laid out on the reference sequence |
| `export_contact_strings(path, ...)` | path | TSV of sequences and contact strings |
| `contact_matrix(chain_a, chain_b, ...)` | dict | Residue-pair distance matrix of one chain pair |
| `block_contact_matrix(chain_a, chain_b, ...)` | dict | Block matrix with intra- and inter-chain masks |
| `contact_statistics(cutoff, moieties, scheme)` | dict of tables | Pairs, per-residue profile, sequence separation, moiety pairs |

### Interfaces (heavy atoms, Å)
| Name | Kind | Description |
|---|---|---|
| `interface_metrics(chain_a, chain_b)` | dict | `{"complex": metrics, "residues": table}` for one pair |
| `interface_metrics_all(chain_pairs=None)` | dict of tables | The same for many pairs |
| `interface_complex(a, b)` / `interface_residues(a, b)` | dict / table | Complex-level / residue-level part of `interface_metrics` |
| `interface_buried_area(a, b)` | float | Buried surface area, Å² |
| `interface_initial_sasa(a, b)` | float | Summed SASA of the isolated chains, Å² |
| `interface_buried_fraction(a, b)` | float | Buried area / initial SASA |
| `interface_atom_contact_count(a, b)` | int | Heavy-atom pairs within the interface cutoff |
| `interface_min_distance_a(a, b)` | float | Minimum heavy-atom distance, Å |
| `interface_contact_score(a, b)` | float | Sum of inverse distances within the score cutoff |
| `interface_partner_labels(a, b)` | dict | Interface residue → partner residues |
| `chain_partner_interfaces(chain, cutoff)` | dict | `interface_metrics` for every contact partner of a chain |
| `chain_partner_structures(chain, cutoff, moieties)` | dict | Two-chain structures for every contact partner |
| `chain_partner_table(chain, cutoff, level="complex")` | table | Interface rows of all partners of a chain |

### Plots and 3D
| Name | Kind | Description |
|---|---|---|
| `plot_contact_map(chain_a, chain_b=None, ...)` | figure | Contact map of one chain or a pair |
| `plot_multi_contact_map(chain_ids, ...)` | figure | Several chains on one block map |
| `show_structure_3d(...)` | HTML | Interactive 3Dmol.js view in Jupyter |

Axis ticks are author residue numbers. `crop=True` zooms a pair map onto the residues in contact, and maps more elongated than 2.5:1 (such as protein–DNA) switch to automatic aspect unless `aspect` is given. Plots default to 3 in (3.5 in for maps of three or more chains); 3D views are 500 px, or 350 px each side by side.

### Export (TSV / PDF / HTML, file names unchanged)
| Name | Writes |
|---|---|
| `export_chain_summary(outdir)` | `{id}_chains.tsv` |
| `export_residue_mapping(chain_ids, outdir)` | `{id}_residue_mapping.tsv` |
| `export_residue_table(outdir)` | `{id}_residues.tsv` |
| `export_moiety_table(outdir)` | `{id}_moieties.tsv` |
| `export_reference_alignment_stats(outdir)` | `{id}_alignment_stats.tsv` |
| `export_reference_alignment_positions(outdir)` | `{id}_alignment_positions.tsv` |
| `export_interface_complex(pairs, outdir)` / `export_interface_residues(pairs, outdir)` | `{id}_interface_complex.tsv` / `{id}_interface_residues.tsv` |
| `export_plot(plot_func, name, ...)` / `export_plot_grid(plot_funcs, name, ...)` | Figures, PDF by default |
| `export_structure_3d(name, ...)` | Standalone HTML 3D view |

### Docking entry point
| Name | Kind | Description |
|---|---|---|
| `prepare_haddock_docking(*partners, workdir, **kwargs)` | `HaddockDocking` | Set up a HADDOCK3 docking of this structure's chains |

## `HaddockDocking` (HADDOCK3)

Guide: [6. Docking with HADDOCK3](guides/06-docking.md).

| Name | Kind | Description |
|---|---|---|
| `HaddockDocking(partners, workdir, method, seed, ...)` | class | Partners, numbering, separation and configuration |
| `haddock_residue_map` | property, table | Written residue ↔ original chain, author number, sequence and reference index |
| `add_ambiguous_restraints(chain_a, active_a, chain_b, active_b, passive_a, passive_b, numbering)` | mutator | Active/passive AIRs |
| `add_unambiguous_restraint(chain_a, positions_a, chain_b, positions_b, numbering, distance)` | mutator | One patch-to-patch restraint |
| `add_patch_pair_restraints(pairs, kind, numbering)` | mutator | One restraint set per patch couple (used by `eval`) |
| `add_contact_restraints(kind, cutoff, scheme, source, chain_map, patch_gap)` | mutator | Restraints from inter-partner contacts |
| `partner_contact_pairs(cutoff, scheme, source, chain_map)` | dict | Contacts between partners, as sequence-index pairs |
| `clear_restraints()` | mutator | Remove all restraints |
| `restraint_table` / `contact_patch_table` | property, table | Restraints added / contact patches behind unambiguous restraints |
| `restraint_file_text(kind)` | str | CNS `.tbl` text of ambiguous or unambiguous restraints |
| `workflow_modules()` / `config_toml()` | list / str | Workflow with seed, restraints, overrides / its TOML |
| `input_paths` / `run_dir` | property | Files written by `write_inputs()` / the haddock3 run folder |
| `write_inputs()` | path | Write PDBs, reference, restraints, residue map, config |
| `validate_config()` | self | Check the configuration with `haddock3 --setup` |
| `run_haddock(overwrite=False)` | self | Run haddock3 in a subprocess |
| `run_steps()` | dict | Step folders of the run |
| `haddock_scores(step=None)` | table | HADDOCK score, clusters and energy terms per model of a step |
| `caprieval_tables(step=None)` | dict of tables | `capri_ss` / `capri_clt` of a caprieval or analysed step |
| `run_caprieval_on_steps(steps)` | dict | Run caprieval on finished steps (`haddock3-analyse`) |
| `docked_model_table(step=None)` | table | Final models: caprieval metrics + scores + energies + path |
| `docked_model_paths(step=None)` | list | Final model files, best first |
| `load_docked_model(i_or_path)` / `load_docked_models(top)` | structure / dict | Docked models as `SyncedStructure` with original chains |
| `cluster_contact_patches(points, gap)` | function | Single-linkage clustering of contact pairs into patches |

## Module-level functions
| Name | Description |
|---|---|
| `search_rcsb(organisms, keywords, title_contains, keywords_field, experimental_method, resolution_max, resolution_min, release_date_from, release_date_to, polymer_entity_count_min, polymer_entity_count_max, uniprot_ids, enzyme_classification, return_type="entry", max_rows=10000, ...)` | Search RCSB entries; values within one criterion are OR-ed, criteria are AND-ed |
| `fetch_rcsb_citations(df, id_col="pdb_id", ...)` | Citations of the PDB IDs in a DataFrame |
| `build_panel_grid(plot_funcs, ...)` / `export_panels_pdf(grids, outdir, ...)` | Lettered panel grids / multi-page PDF of grids |
| `show_side_by_side_3d(*views, ...)` | Several 3D views in a row |

## Configuration

All loaders accept `config_path=`, a YAML file whose `synced_structure` section overrides these defaults. The file is never written.

```yaml
synced_structure:
  save_outdir: data/saved          # default output folder for exports
  aln_match_score: 1               # local alignment scoring (reference_alignment_*, sequence matching)
  aln_mismatch_score: -1
  aln_open_gap_score: -2
  aln_extend_gap_score: -0.5
  interface_contact_cutoff_a: 4.5  # Å, heavy-atom contact cutoff of interface metrics
  interface_score_cutoff_a: 12.0   # Å, distance range of the interface contact score
  interface_block_size: 1000       # rows per block when computing interface atom distances (memory)
```

## Renamed names
The old names were removed; update existing scripts with this table.

| Old | New |
|---|---|
| `subset`, `match_chains`, `rename_by_sequence`, `summary` | `subset_chains`, `find_chains_by_sequence`, `rename_chains_by_sequence`, `chain_summary` |
| `residue_map`, `moiety_summary`, `reindex_chain` | `residue_table`, `moiety_table`, `renumber_chain` |
| `alignment_stats`, `alignment_positions` | `reference_alignment_stats`, `reference_alignment_positions` |
| `to_ca_trajectory`, `get_contact_sequence`, `get_reference_contact_sequence`, `export_contact_sequences` | `contact_point_trajectory`, `contact_string`, `reference_contact_string`, `export_contact_strings` |
| `plot_contacts_multi` | `plot_multi_contact_map` |
| `interface`, `interface_all`, `chain_partner_frame` | `interface_metrics`, `interface_metrics_all`, `chain_partner_table` |
| `buried_area`, `initial_sasa`, `buried_fraction`, `contact_count`, `min_distance`, `contact_score` | `interface_buried_area`, `interface_initial_sasa`, `interface_buried_fraction`, `interface_atom_contact_count`, `interface_min_distance_a`, `interface_contact_score` |
| `export_chains`, `export_residues`, `export_moieties`, `export_alignment_stats`, `export_alignment_positions` | `export_chain_summary`, `export_residue_table`, `export_moiety_table`, `export_reference_alignment_stats`, `export_reference_alignment_positions` |
| `search_rcsb_pubs`, `side_by_side` | `fetch_rcsb_citations`, `show_side_by_side_3d` |
| `jupy_viz`, `complete_contacts`, `partners`, `from_rcsb_id`, `residue_moiety` | `show_structure_3d`, `block_contact_matrix`, `interface_partner_labels`, `from_rcsb`, `residue_class` |
| `record["field"]` on a chain record | `record.field` |
| `HaddockSetup`, `prepare_haddock` and its methods | `HaddockDocking`, `prepare_haddock_docking`; see the HaddockDocking table |
