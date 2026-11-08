"""Dock actin and vitamin D-binding protein (1KXP) with the 'ab-initio' protocol and compare the models with the crystal.

Needs the haddock extra and the haddock3 executable on PATH. Run from anywhere with:
    python examples/docking/1kxp_ab_initio.py

The partners are pulled apart with random rotations and docked without any knowledge of the interface: rigidbody
draws random ambiguous restraints between surface patches (ranair) and keeps the partners together with
center-of-mass restraints (cmrest), emref refines with restraints from the contacts it finds (contactairs). The
crystal is used only by caprieval. The run is downsized to 1,000 rigid-body models for a laptop; set
MODULE_PARAMS = {} for the full preset (10,000). Every argument of the docking calls is written out, at its default
unless the comment says otherwise.

Outputs land in examples/output/docking_1KXP_ab-initio/:
    run/        the HADDOCK3 inputs, config.toml and run folder
    tables/     residue map, model scores, caprieval models and clusters vs the crystal, and the interface of the
                crystal vs the top models (TSV)
    figures/    HADDOCK score vs i-RMSD colored by DockQ (PDF)
    structure/  interactive 3D view of the best model as standalone HTML
"""

from pathlib import Path

import pandas as pd

from synced_structure import SyncedStructure, core_count

PDB_ID = "1KXP"
ACTIN, DBP = "A", "D"
METHOD = "ab-initio"
SEED = 7
# laptop-sized run: 1,000 rigid-body models instead of 10,000; {} keeps the preset
MODULE_PARAMS = {"rigidbody": {"sampling": 1000}, "clustfcc": {"min_population": 1}}
NCORES = core_count()           # every core granted to this job (affinity, cgroup quota, scheduler), not the whole node
TOP = 5                         # best models compared with the crystal
CAPRI_COLUMNS = ["model", "score", "irmsd", "fnat", "lrmsd", "dockq", "cluster_id"]
OUTDIR = Path(__file__).resolve().parents[1] / "output" / f"docking_{PDB_ID}_{METHOD}"
RUN, TABLES, FIGURES, STRUCTURE = (OUTDIR / name for name in ("run", "tables", "figures", "structure"))

###################################### PREPARE

s = SyncedStructure.from_rcsb(PDB_ID)
for folder in (TABLES, FIGURES, STRUCTURE):
    folder.mkdir(parents=True, exist_ok=True)

dock = s.prepare_haddock_docking(
    ACTIN, DBP,                         # one partner per argument: a chain, a list of chains, or (structure, chains)
    workdir=RUN,
    method=METHOD,                      # 'ab-initio' (default), 'oracle', 'targeted' or 'alascan'
    seed=SEED,                          # default 917
    ncores=NCORES,                      # default None: core_count()
    separate="auto",                    # move partners apart when they share a structure or touch
    clearance_a=10.0,                   # smallest gap in Å between separated partners
    module_params=MODULE_PARAMS,        # default None: the preset
    global_params=None,                 # top-level haddock3 parameters, such as {"max_cpus": True}
    reference="auto",                   # the input complex, since both partners come from 1KXP
    caprieval=True,                     # end with caprieval against the reference
    haddock3_exe=None,                  # search PATH and this Python's bin directory
)
# no restraints: adding ambiguous ones (add_ambiguous_restraints, add_contact_restraints) turns ranair off
dock.export_inputs()
dock.validate_config()                  # haddock3 --setup in a scratch folder: catches unknown module_params

###################################### RUN

dock.run_haddock(overwrite=True, progress=True)

###################################### TABLES (TSV)

models = dock.load_docked_models(top=TOP, step=None)   # step=None: the last step
interfaces = pd.DataFrame(
    [{"model": PDB_ID, **s.interface_complex(ACTIN, DBP)}]
    + [{"model": name, **model.interface_complex(ACTIN, DBP)} for name, model in models.items()])

tables = {
    "residue_map": dock.residue_table(),
    "scores": dock.score_table(step=None),
    "docked_models": dock.docked_model_table(step=None),   # final models with irmsd, fnat, lrmsd, dockq vs the crystal
    "interfaces": interfaces,
}
capri = dock.caprieval_tables(step=None)    # {"models": per model, "clusters": per cluster} vs the crystal
tables.update({f"caprieval_{name}": table for name, table in capri.items()})
table_paths = []
for name, table in tables.items():
    path = TABLES / f"{PDB_ID}_{METHOD}_{name}.tsv"
    table.to_csv(path, sep="\t", index=False)
    table_paths.append(path)

###################################### FIGURES (PDF) AND 3D VIEW (HTML)

figure_paths = s.export_plot(
    lambda: dock.plot_caprieval(x="irmsd", y="score", color="dockq", step=None,
                                title=f"Actin – DBP vs {PDB_ID} ({METHOD})", ax=None),
    "caprieval", outdir=FIGURES)

best_name, best = next(iter(models.items()))
structure_paths = [best.export_structure_3d(outdir=STRUCTURE, title=f"Actin – DBP, best docked model ({best_name})")]

print(f"caprieval of the top {TOP} models against the crystal:")
print(capri["models"].head(TOP)[[column for column in CAPRI_COLUMNS if column in capri["models"]]].to_string(index=False))
for folder, paths in ((TABLES, table_paths), (FIGURES, figure_paths), (STRUCTURE, structure_paths)):
    print(f"{folder.relative_to(OUTDIR.parent)}/")
    for path in paths:
        print(f"  {Path(path).name}")
