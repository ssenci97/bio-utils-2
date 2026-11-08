# Quick start

This walks through the pipeline on 1AAY: the three zinc fingers of Zif268 bound to DNA. Every step returns a table you can inspect; figures and files are made from those tables.

## Load

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")          # download from the RCSB
# s = SyncedStructure.from_path("model.cif")   # or load a local .pdb / .cif / .mmcif file
```

## Chains

```python
chains = s.chain_table()
print(chains[["chain_key", "chain_id", "auth_chain_id", "moiety", "description"]])
```

Chains are stored under their author IDs (`chain_key`). In 1AAY, author chain `A` is the protein and `B`, `C` are the two DNA strands, while the label IDs (`chain_id`) are different letters. Every method accepts any of the three names; see [Concepts](guide/concepts.md#chain-ids).

## Numbering

```python
residues = s.chain_residue_table("A")          # author number ↔ observed index ↔ reference index
print(residues.head(3))
alignment = s.alignment_table()                # identity and coverage of each chain against its reference
```

To renumber chains by reference position and drop unmodelled residues, call `s.apply_reference_numbering()`.

## Contacts

```python
contacts = s.contact_table(cutoff=0.8)          # inter-chain residue pairs closer than 0.8 nm
stats = s.contact_stat_tables()                  # per pair, per residue, by separation, by moiety pair
print(stats["pairs"][["chain_key_a", "chain_key_b", "contact_type", "contact_pairs"]])
print(s.contact_string("A"))                     # 1 = contact, 0 = none, X = no contact point
print(s.masked_sequence("A"))                    # residues in contact, others masked with '-'
```

## Interfaces

```python
metrics = s.interface_complex("A", "B")          # buried area (Å²), interface residues, contact score
print(round(metrics["dASA_A2"]), "Å² buried between A and B")
iface = s.interface_residue_table("A", "B")      # one row per interface residue
```

## Figures

```python
fig, ax = s.plot_contact_map("A", "B", crop=True)           # contact map of one chain pair
fig, ax = s.plot_multi_contact_map(["A", "B", "C"])         # all three chains on one map
view = s.show_structure_3d()                                # interactive 3D view; display it in Jupyter
```

## Files

```python
s.export_structure(outdir="results")                        # results/1AAY.cif; reloads with the same chains
s.export_structure("pdb", chain_ids=["A"], outdir="results")   # results/1AAY_A.pdb
s.export_chain_table(outdir="results")                      # results/1AAY_chains.tsv
s.export_contact_stat_tables(outdir="results")              # results/1AAY_contact_stats_*.tsv
s.export_plot(lambda: s.plot_contact_map("A", "B"), "map_A_B", outdir="results")   # results/1AAY_map_A_B.pdf
s.export_structure_3d(outdir="results")                     # results/1AAY_structure_3d.html
```

A script that writes every table, contact map and a 3D page for a larger complex is in [`examples/export_all.py`](https://github.com/ssenci97/bio-utils-2/blob/main/examples/export_all.py).

## Next

- [Concepts](guide/concepts.md): numbers, chain IDs and the vocabulary used in every table.
- [Conventions](guide/conventions.md): how names, columns and units are chosen.
- [HADDOCK3 docking](guide/docking.md): dock chains with restraints in your own numbering.
- [API reference](reference/index.md): every method with its signature.
