# Comparison

Compare a structure (the model) with a reference: matched residues, RMSD after an optimal (Kabsch) superposition, and the CAPRI quality metrics of docking (i-RMSD, l-RMSD, Fnat, DockQ). Distances are in **Å**. Residues match through their reference positions, through a global alignment when the reference sequences differ, and through the observed sequences when a chain has no reference (PDB files without SEQRES). See the [comparison guide](../guide/comparison.md).

::: synced_structure.SyncedStructure.residue_match_table

::: synced_structure.SyncedStructure.rmsd_metrics

::: synced_structure.SyncedStructure.apply_superposition

::: synced_structure.SyncedStructure.capri_metrics

::: synced_structure.SyncedStructure.capri_table
