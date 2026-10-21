# synced_structure: Philosophy and API Overview

## Philosophy

`synced_structure` does not try to replace Biopython, MDTraj or the RCSB APIs. Those tools are more general and better tested, and they cover far more ground than this package. Biopython parses files and stores the atoms, MDTraj computes distances quickly, and RCSB provides the data. This package works on top of them and tackles one problem they leave open: **a structure on its own does not tell you what its residues mean.**

A deposited chain is a partial, renumbered and sometimes mislabelled copy of a biological sequence. It may have missing loops, insertion codes, author IDs that differ from label IDs, modified residues and ligands mixed in. Biopython stores what was observed. `synced_structure` links each observed residue to its **reference** sequence (the full mmCIF entity sequence, or one you supply) and keeps that link consistent through every later operation.

Four principles follow from this:

1. **Wrap, don't fork.** The underlying `Bio.PDB.Structure` stays reachable as `.structure` and is never hidden. Anything Biopython can do still works, and `SyncedStructure` only adds meaning on top.
2. **One mapping, kept in sync.** Each chain has a single observed→reference map, built by global alignment. Changing a reference, renaming a chain, replacing a chain or renumbering residues updates that map instead of leaving it stale. This is what *synced* means.
3. **Identity is plural.** A chain can be addressed by its key, its label (`chain_id`) ID or its author (`auth_chain_id`) ID. You choose the key with `key_by`, and the tables always report both IDs. Nothing silently assumes PDB-style author numbering.
4. **Tables first, figures second.** Every analysis returns a `pandas.DataFrame` or a dict and can be exported. Plots and 3D views are built from those same tables, so a figure always matches the numbers behind it.

The design is a single class assembled from focused mixins (loading, chains, alignment, residues, contacts, interfaces, plots, export and 3D). Each mixin handles one area, and they all share the same synced chain state. Mutating methods return `self`, so calls can be chained.

---

## Construction

```python
SyncedStructure(structure, structure_id="structure", mmcif_dict=None, file_format="cif",
                config_path=DEFAULT_CONFIG_PATH, default_entity_id=None, key_by="auth")
SyncedStructure.from_rcsb_id(pdb_id, config_path=..., key_by="auth")
SyncedStructure.from_path(path, config_path=..., key_by="auth")
SyncedStructure.from_biopython(structure, structure_id=None, mmcif_dict=None, file_format="pdb",
                               config_path=..., key_by="auth")
```

Attributes: `structure`, `model`, `structure_id`, `file_format`, `mmcif_dict`, `config`, `key_by`, `chains`.

## Chains and references (remapping core)

| Method | Purpose |
|---|---|
| `chain_ids` *(property)* | Loaded chain keys |
| `summary()` | Label/author IDs, entity, description, external ref, refseq ID, coverage |
| `get_sequence(chain_id)` / `get_reference_sequence(chain_id)` | Observed / reference sequence |
| `set_reference_sequence(chain_id, sequence, remap=True)` | Replace one reference |
| `set_reference_sequences(sequences, remap=True)` | Replace several references from `{chain: seq}` |
| `remap_chain(chain_id)` / `remap_all_chains()` | Rebuild mappings |
| `get_external_ref(chain_id)` / `set_external_ref(chain_id, ref)` | Accession (for example UniProt) |
| `get_refseq_id(chain_id)` / `set_refseq_id(chain_id, ref)` | ID of the reference sequence; falls back to the external ref |
| `set_entity_id(chain_id, entity_id)` / `set_chain_metadata(chain_id, key, value)` | Annotate a chain |
| `add_chain(chain, reference_sequence=None, entity_id=None, metadata=None, replace=False)` | Add a chain |
| `replace_chain(chain_id, chain, ...)` / `remove_chain(chain_id)` | Replace or remove a chain |
| `match_chains(sequence, min_identity=90.0)` | Find chains matching a sequence |
| `rename_chains(mapping)` / `rename_by_sequence(targets, min_identity=90.0)` | Rename while keeping state in sync |
| `chain_contact_partners(cutoff=..., moieties=...)` | `{chain: [partners]}` |

## Alignment

`alignment_stats` *(property)*: local alignment statistics per chain. `alignment_positions` *(property)*: the column-by-column alignment.

## Residues

| Method | Purpose |
|---|---|
| `get_residue_mapping(chain_id)` | pdb_resnum ↔ sequence_index ↔ reference_index |
| `residue_map` / `moiety_summary` *(properties)* | All residues; counts per chain and moiety |
| `chains_by_moiety(moiety)` | Chains containing `protein` / `nucleic` / `ligand` / `ion` / `water` |
| `apply_reference_numbering(chain_ids=None)` | Trim unmapped residues and renumber by reference, in place |
| `reindex_chain(chain_id, in_place=False)` | Renumber 1..N |

## Contacts

```python
to_ca_trajectory(chain_ids=None, moieties=...)
compute_contacts(cutoff=..., contact_type=..., chain_ids=None, scheme="ca", moieties=...)
get_contact_sequence(chain_id, ...) / get_reference_contact_sequence(chain_id, ...)
export_contact_sequences(output_path, ..., reference=False, ...)
contact_matrix(chain_a="A", chain_b=None, cutoff=..., moieties=..., scheme="ca")
complete_contacts(chain_a="A", chain_b="B", ..., min_separation=...)
contact_statistics(cutoff=..., moieties=..., scheme="ca")
```

## Interfaces

`interface(a, b)`, `interface_all(chain_pairs=None)`, `interface_complex(a, b)`, `interface_residues(a, b)`, `chain_partner_interfaces(chain_id, cutoff=...)`, `chain_partner_frame(chain_id, cutoff=..., level="complex")`, plus the scalar metrics `buried_area`, `initial_sasa`, `buried_fraction`, `contact_count`, `min_distance`, `contact_score` and `partners`, which all take `(chain_a="A", chain_b="B")`.

## Visualization

`plot_contact_map(chain_a="A", chain_b=None, ...)`, `plot_contacts_multi(chain_ids=("A","B"), ..., plddt=False, ...)`, `plot_contact_panel(chain_a="A", chain_b=None, ...)` and `jupy_viz(angles=(0,0,0), chain_id=None, colors=None, title=None, style="cartoon", ...)`.

## Export

All of these take `outdir=None`: `export_chains`, `export_residue_mapping(chain_ids=None)`, `export_residues`, `export_moieties`, `export_alignment_stats`, `export_alignment_positions`, `export_interface_complex(chain_pairs=None)`, `export_interface_residues(chain_pairs=None)`, `export_plot(plot_func, name, formats=..., dpi=..., raster=True)` and `export_plot_grid(plot_funcs, name, ...)`.

## Module-level functions (`from synced_structure import ...`)

```python
search_rcsb(organisms=None, keywords=None, title_contains=None, keywords_field=None,
            experimental_method=None, resolution_max=None, resolution_min=None,
            release_date_from=None, release_date_to=None, polymer_entity_count_min=None,
            polymer_entity_count_max=None, uniprot_ids=None, enzyme_classification=None,
            return_type="entry", max_rows=10000, return_all_hits=False,
            include_computational=False, sort_by="score", sort_direction="desc", timeout=30)
search_rcsb_pubs(df, id_col="pdb_id", primary_only=True, merge=True, batch_size=200, timeout=30)
build_panel_grid(plot_funcs, sizes=3.5, ncols=2, panel_labels=True, main_title=None, wspace=0.5, hspace=0.5)
export_panels_pdf(grids, outdir, filename=..., raster=True, dpi=..., **grid_kwargs)
side_by_side(*outputs, gap_px=...)
```
