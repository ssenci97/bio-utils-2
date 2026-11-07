# synced_structure

**A Biopython structure wrapper that keeps every observed residue linked to its position in the reference sequence**, and builds contacts, contact maps, interface metrics, HADDOCK3 dockings, tables and 3D views on top of that link.

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")            # Zif268 zinc fingers bound to DNA
s.chain_table()                                  # chains with label and author IDs, entity, coverage
s.contact_table(cutoff=0.8)                      # inter-chain residue contacts, in nm
s.interface_complex("A", "B")                    # buried surface and interface metrics, in Å
fig, ax = s.plot_contact_map("A", "B")           # contact map
```

<div class="grid cards" markdown>

-   **Get started**

    [Install](install.md) the package and follow the [quick start](quickstart.md) through the whole pipeline.

-   **Understand the tables**

    [Concepts](guide/concepts.md) explains the three residue numbers and the two chain IDs every table reports.

-   **Find a method**

    The [API reference](reference/index.md) lists every function by pipeline step, with signatures and examples.

-   **Dock**

    The [HADDOCK3 guide](guide/docking.md) sets up, runs and reads back a docking in your own numbering.

</div>

## Why

A deposited structure is an incomplete copy of its sequence: residues are missing, numbering is the depositors' choice, and each chain has two IDs. `synced_structure` aligns every chain to its reference sequence once, keeps that mapping in sync through every edit, and reports reference positions next to author numbers in every table. Contacts, interfaces and docking restraints from different structures then refer to the same residues.

## At a glance

- **Tables first**: every analysis returns a `pandas.DataFrame` (or a dict of them), with chain keys, label and author IDs, and 1-based indices.
- **Names say what they return**: `_table` is a DataFrame, `_map` a dict, `export_` writes files and returns their paths ([Conventions](guide/conventions.md)).
- **Units**: nm for contacts, Å for interfaces and docking.
