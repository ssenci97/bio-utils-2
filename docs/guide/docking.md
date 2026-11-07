# HADDOCK3 docking

`HaddockDocking` sets up a [HADDOCK3](https://www.bonvinlab.org/haddock3/) run from your chains, writes restraints in the numbering you choose, runs `haddock3` in a subprocess, and reads scores and models back as `SyncedStructure` objects. Install the extra first: `pip install "synced_structure[haddock] @ git+https://github.com/ssenci97/bio-utils-2.git"`.

## 1. Prepare

Each positional argument is one **partner** (rigid body): one chain, a list of chains, or `(structure, chains)` for a chain of another structure.

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")
dock = s.prepare_haddock_docking("A", ["B", "C"], workdir="dock", method="oracle", seed=7)
dock.residue_table().head()          # partner, segid, haddock_resid ↔ chain_key, pdb_resnum, reference_index
```

| `method` | Workflow |
|---|---|
| `ab-initio` | Full sampling (10,000 rigid-body models), refinement, clustering; the default. |
| `oracle` | A short run (200 models) that checks restraints recover a known complex. |
| `eval` | Like `oracle`, always ending with `caprieval` against the reference. |

Partners that come from the same structure are moved apart with random rotations (`separate="auto"`), so the docking does not start from the answer. Override the preset of any module with `module_params={"rigidbody": {"sampling": 1000}}`.

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

## 4. Read the results

<!-- skip-test -->
```python
scores = dock.score_table()                  # every model of the last step, best first
final = dock.docked_model_table()            # final models with caprieval metrics, if caprieval ran
models = dock.load_docked_models(top=5)      # {name: SyncedStructure} with your chain keys and residue IDs
best = dock.load_docked_model(0)
best.interface_complex("A", "B")             # analyse it like any other structure
```

Docked models keep the reference sequences and numbering of the input chains, so contacts and interfaces of a model compare directly with the input.
