# HADDOCK3 docking

`HaddockDocking` sets up a [HADDOCK3](https://www.bonvinlab.org/haddock3/) run from your chains, writes restraints in the numbering you choose, runs `haddock3` in a subprocess, and reads scores and models back as `SyncedStructure` objects. Install the extra first: `pip install "synced_structure[haddock] @ git+https://github.com/ssenci97/bio-utils-2.git"`.

## 1. Prepare

Each positional argument is one **partner** (rigid body): one chain, a list of chains, or `(structure, chains)` for a chain of another structure. Partners can come from different PDB entries or models. Their residues are matched through reference positions, and `set_reference(complex, chain_map=...)` names a complex to take contact restraints from and to evaluate against ([tutorial](../tutorials/haddock/03-different-structures.md)).

The [HADDOCK3 tutorials](../tutorials/haddock/index.md) go through every step in detail.

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")
dock = s.prepare_haddock_docking("A", ["B", "C"], workdir="dock", method="oracle", seed=7)
dock.residue_table().head()          # partner, segid, haddock_resid ↔ chain_key, pdb_resnum, reference_index
```

| `method` | Workflow |
|---|---|
| `ab-initio` | 10,000 rigid-body models, with random and center-of-mass restraints unless you add your own; the default. |
| `oracle` | 200 rigid-body models with restraints from the contacts of the input complex; checks the set up recovers a known complex. |
| `targeted` | 200 rigid-body models with restraints you give as patch pairs (`add_patch_pair_restraints`). |

Every method then clusters the rigid-body models (`clustfcc`), keeps the best fifth of the clusters (at least 5, or all when there are fewer) with 4 models each (`seletopclusts`), and refines only those (`emref`). It ends with `caprieval`, against the reference when there is one and otherwise against the best-scoring model. The number of clusters is known only after clustering, so `run_haddock` runs haddock3 up to `clustfcc`, counts the clusters, and restarts the run from `seletopclusts`. Setting `module_params={"seletopclusts": {"top_clusters": 10}}` fixes the number and runs haddock3 once.

Partners that come from the same structure, or that touch, are moved apart with random rotations (`separate="auto"`), so the docking does not start from the answer. Override the preset of any module with `module_params={"rigidbody": {"sampling": 1000}}`.

## 2. Add restraints

Restraints are given on the original chains, in `reference` (default), `sequence` or `author` numbering, and are translated to HADDOCK3 numbering through the residue map.

```python
dock.add_contact_restraints(kind="ambiguous")    # active residues from the contacts of the input complex
dock.restraint_table()
```

The other options are `add_ambiguous_restraints` (active/passive residues), `add_unambiguous_restraint` (residue-to-residue distances) and `add_patch_pair_restraints`. `clear_restraints` starts over.

## 3. Write, check and run

```python
config = dock.export_inputs()        # partner PDBs, restraint .tbl files, residue map and config.toml
print(dock.config_toml()[:200])
```

<!-- skip-test -->
```python
dock.validate_config()               # haddock3 --setup in a scratch folder
dock.run_haddock()                   # runs haddock3; the run folder is dock/run
```

`run_haddock` logs progress while haddock3 runs: each module's start and finish, and while a module writes models, the count so far with an estimate of the time left. Interrupting Python (Ctrl+C, or stopping a notebook cell) stops haddock3. Pass `progress=False` to log only the end of the run.

```text
HADDOCK3 [2/5] rigidbody started | run 14s
HADDOCK3 [2/5] rigidbody: 90/200 models | 11m40s | ~14m16s left
HADDOCK3 [2/5] rigidbody finished in 25m53s
```

## 4. Read the results

<!-- skip-test -->
```python
scores = dock.score_table()                  # every model of the last step, best first
final = dock.docked_model_table()            # final models with caprieval metrics, if caprieval ran
fig, ax = dock.plot_caprieval()              # HADDOCK score vs i-RMSD, colored by DockQ
models = dock.load_docked_models(top=5)      # {name: SyncedStructure} with your chain keys and residue IDs
best = dock.load_docked_model(0)
best.interface_complex("A", "B")             # analyse it like any other structure
```

Docked models keep the reference sequences and numbering of the input chains, so contacts and interfaces of a model compare directly with the input.

A complete run is in [`examples/dock_1kxp.py`](https://github.com/ssenci97/bio-utils-2/blob/main/examples/dock_1kxp.py): it re-docks actin and vitamin D-binding protein (1KXP) with the `oracle` preset and writes the restraints, scores, caprieval metrics, the interfaces of the top models next to the crystal's, and a 3D page of the best model.
