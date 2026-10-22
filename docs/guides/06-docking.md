# 6. Docking with HADDOCK3

**Goal:** dock chains with [HADDOCK3](https://www.bonvinlab.org/haddock3/), driving it with restraints written on your own chains and numbering, and get the docked models back as `SyncedStructure` objects with scores.

**You need:**
- [guides 1–3](01-getting-started.md), so that chain keys and contacts are familiar;
- haddock3 for the run itself: `pip install haddock3` (its wheels include the CNS engine on Linux and macOS). Steps 1–4 work without it.

haddock3 runs in a separate process. Its command line calls `sys.exit` on errors, which would end a notebook session; here an error raises `RuntimeError` with the end of the log instead.

## Step 1: choose the partners

Each **partner** is one rigid body. Pass one argument per partner: a chain, or a list of chains that move together.

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")
dock = s.prepare_haddock_docking(
    "A", ["B", "C"],                 # partner 1: the protein; partner 2: the DNA duplex (two strands)
    workdir="runs/1aay",
    method="oracle",
    seed=7,
    ncores=2,
    module_params={"rigidbody": {"sampling": 4}, "seletop": {"select": 2}},   # tiny run for a laptop
)
print(dock.haddock_residue_map.groupby(["segid", "chain_key"])["haddock_resid"].agg(["min", "max"]))
```

```text
                 min  max
segid chain_key
A     A          103  187
B     B            1   11
      C           51   61
```

- Each partner gets one segment ID (`A`, `B`, …) shared by all its chains: haddock3 moves every segment ID separately, so a duplex written as two IDs would be pulled apart.
- Author numbers are kept when they are unique within a partner; otherwise the partner is renumbered. `haddock_residue_map` (also saved as `haddock_residue_map.tsv`) links every written residue back to its chain and numbers.
- Partners from the same structure are rotated randomly and moved at least `clearance_a` = 10 Å apart, so the docking does not start from the answer. `seed` makes this reproducible and is also HADDOCK's `iniseed`.
- To dock chains from different files, pass `(other_structure, chains)` as a partner.

## Step 2: add restraints

Restraints tell HADDOCK which residues should meet. Write positions on your chains; `numbering` says what they mean: `"reference"` (default), `"sequence"` or `"author"`.

```python
# from the contacts of the input complex, grouped into patches
dock.add_contact_restraints(kind="unambiguous", scheme="sidechain", cutoff=0.45)
print(dock.contact_patch_table[["patch", "chain_a", "residues_a", "chain_b", "residues_b"]].head(3))

# or by hand: residues 44 and 47 of the protein's reference sequence should touch nucleotides 7-8 of strand B
dock.add_ambiguous_restraints("A", [44, 47], "B", [7, 8], numbering="reference")
print(dock.restraint_table[["kind", "origin", "segid_a", "residues_a", "segid_b", "residues_b"]].tail(1))
```

```text
   patch chain_a            residues_a chain_b         residues_b
0      1       A  [118, 120, 121, 124]       B  [7, 8, 9, 10, 11]
1      2       A  [146, 148, 149, 152]       B       [5, 6, 7, 8]
2      3       A       [174, 177, 180]       B    [1, 2, 3, 4, 5]
        kind  origin segid_a  residues_a segid_b residues_b
6  ambiguous  custom       A  [144, 147]       B     [7, 8]
```

Residue numbers in these tables are HADDOCK numbers (see `haddock_residue_map`).

| Method | Adds |
|---|---|
| `add_ambiguous_restraints(chain_a, active_a, chain_b, active_b, passive_a=(), passive_b=())` | Ambiguous interaction restraints (AIRs): each active residue must touch any active or passive residue of the other side |
| `add_unambiguous_restraint(chain_a, positions_a, chain_b, positions_b, distance=(2.0, 2.0, 0.0))` | One restraint joining two patches |
| `add_patch_pair_restraints([((chain_a, positions), (chain_b, positions)), ...], kind="ambiguous")` | One restraint set per couple of patches |
| `add_contact_restraints(kind, cutoff, scheme, source=None, chain_map=None, patch_gap=4)` | Restraints from contacts; `"unambiguous"` groups contacts into patches whose positions differ by at most `patch_gap` on both sequences. A different `source` complex is transferred through reference positions |

Positions with no observed residue are logged and dropped. `clear_restraints()` starts over.

## Step 3: pick a method

| `method` | Workflow (HADDOCK modules) | Restraints |
|---|---|---|
| `ab-initio` | topoaa → rigidbody (10 000 models, `cmrest` + `ranair`) → seletop 200 → emref (`contactairs`) → clustfcc → seletopclusts (10 clusters) → caprieval | none needed; your ambiguous restraints replace the random ones (`ranair`) |
| `oracle` | topoaa → rigidbody (200) → emref → seletop 10 → caprieval | contact patches are added automatically if you add none |
| `eval` | as `oracle` | your `add_patch_pair_restraints` couples, as ambiguous restraints; needs a reference |

The presets are sized for real dockings. On a laptop, shrink them with `module_params` as in Step 1 (for `ab-initio`, also `{"clustfcc": {"min_population": 1}}`). Inspect the final workflow with:

```python
print([name for name, _ in dock.workflow_modules()])
```

```text
['topoaa', 'rigidbody', 'emref', 'seletop', 'caprieval']
```

`caprieval` compares the models with the input complex when all partners come from one structure, otherwise with the best-scoring model.

## Step 4: write the inputs

```python
config = dock.write_inputs()      # runs/1aay/config.toml, inputs/*.pdb, restraints/*.tbl, haddock_residue_map.tsv
print(config.name)
```

```text
config.toml
```

Everything HADDOCK needs is now in `runs/1aay`. Restraint files are checked with HADDOCK's own validator.

## Step 5: run HADDOCK3

<!-- skip-test -->
```python
dock.validate_config()            # haddock3 --setup: checks every parameter without running (a few seconds)
dock.run_haddock()                # about a minute with the tiny settings above
```

The run goes to `runs/1aay/run`, with the haddock3 log in `runs/1aay/haddock3.log`. `run_haddock(overwrite=True)` replaces a previous run, but only one this class created.

## Step 6: read scores and models

<!-- skip-test -->
```python
table = dock.docked_model_table()          # one row per final model, best first
print(table[["model", "score", "irmsd", "dockq", "path"]])

best = dock.load_docked_model(0)           # the best model, with the original chains and numbering
print(best.chain_ids)

rigid = dock.haddock_scores("rigidbody")   # scores and energy terms of every model of any step
dock.run_caprieval_on_steps(["rigidbody"]) # CAPRI metrics for an intermediate step...
early = dock.caprieval_tables("rigidbody") # ...then read like the final ones
```

| Method | Returns |
|---|---|
| `docked_model_table(step=None)` | Final models: CAPRI metrics (iRMSD, fnat, DockQ, …), HADDOCK score, energy terms, cluster and file path |
| `docked_model_paths()` / `load_docked_models(top=None)` | The same models as file paths / as `SyncedStructure` objects |
| `haddock_scores(step)` | HADDOCK score and energies for every model of a step |
| `caprieval_tables(step)` | caprieval's `capri_ss` and `capri_clt` tables |
| `run_caprieval_on_steps(steps)` | Runs caprieval on finished intermediate steps |
| `run_steps()` | The step folders of the run |

A `step` is a folder name (`"3_emref"`), a module name (`"emref"`, its last step) or a step number (`3`).

## Common mistakes

- **Running the default presets on a laptop.** `ab-initio` samples 10 000 models; shrink `sampling` and `select` first.
- **Restraint positions in the wrong numbering.** The default is the *reference* sequence; pass `numbering="author"` for the numbers you see in a viewer.
- **Splitting a duplex.** Pass both DNA strands as one partner (`["B", "C"]`), otherwise HADDOCK docks them apart.
- **Reading HADDOCK numbers as author numbers.** Restraint and patch tables use HADDOCK numbering; `haddock_residue_map` translates.

## Next

- [API index](../api_index.md#haddockdocking-haddock3) for every `HaddockDocking` method and argument.
