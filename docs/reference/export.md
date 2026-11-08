# Export

Every `export_` method writes `{structure_id}_{name}` into `outdir` (default: the `save_outdir` setting, see [Configuration](../guide/configuration.md)) and returns the written path or paths. `export_structure` writes `{structure_id}.cif` or `.pdb` instead.

## Structure files

mmCIF reloads with `from_path()` into the same chains: keys, label and author IDs, entities, descriptions, accessions, and reference sequences. PDB keeps coordinates, residue IDs, and author chain IDs only.

::: synced_structure.SyncedStructure.export_structure

::: synced_structure.SyncedStructure.mmcif_string

::: synced_structure.SyncedStructure.pdb_string

## Tables, figures, and pages

::: synced_structure.SyncedStructure.export_chain_table

::: synced_structure.SyncedStructure.export_chain_residue_table

::: synced_structure.SyncedStructure.export_residue_table

::: synced_structure.SyncedStructure.export_moiety_table

::: synced_structure.SyncedStructure.export_alignment_table

::: synced_structure.SyncedStructure.export_alignment_column_table

::: synced_structure.SyncedStructure.export_contact_strings

::: synced_structure.SyncedStructure.export_contact_stat_tables

::: synced_structure.SyncedStructure.export_interface_tables

::: synced_structure.SyncedStructure.export_plot

::: synced_structure.SyncedStructure.export_plot_grid

::: synced_structure.SyncedStructure.export_structure_3d

## Module functions

::: synced_structure.export_panels_pdf
