"""Dock unbound barnase (1A2P) and barstar (1A19) with the 'ab-initio' protocol and evaluate the models against 1BRS.

Needs the haddock extra and the haddock3 executable on PATH. Run from anywhere with:
    python examples/docking/barnase_barstar_ab_initio.py

The partners come from two different PDB entries and are docked without any knowledge of the interface: rigidbody
draws random ambiguous restraints between surface patches (ranair) and keeps the partners together with
center-of-mass restraints (cmrest), emref refines with restraints from the contacts it finds (contactairs). The bound
complex 1BRS is used only by caprieval, written in the partners' HADDOCK3 numbering (1A19 barstar has an extra Met
and Cys40, so its sequence is aligned to 1BRS's first). The run is downsized to 1,000 rigid-body models for a laptop;
set MODULE_PARAMS = {} for the full preset (10,000). Every argument of the docking calls is written out, at its
default unless the comment says otherwise.

Outputs land in examples/output/docking_barnase_barstar_ab-initio/:
    run/        the HADDOCK3 inputs, config.toml and run folder
    tables/     residue map, model scores, caprieval models and clusters vs 1BRS, and the interface of 1BRS vs the
                top models (TSV)
    figures/    HADDOCK score vs i-RMSD colored by DockQ (PDF)
    structure/  interactive 3D view of the best model as standalone HTML
"""

from pathlib import Path

import pandas as pd

from synced_structure import HaddockDocking, SyncedStructure, core_count

NAME = "barnase_barstar"
APO_BARNASE, APO_BARSTAR, HOLO = "1A2P", "1A19", "1BRS"
# partner chains (partner number, chain) → chains of the complex: barnase A–A, barstar A–D
CHAIN_MAP = {(1, "A"): "A", (2, "A"): "D"}
METHOD = "ab-initio"
SEED = 11
# laptop-sized run: 1,000 rigid-body models instead of 10,000; {} keeps the preset
MODULE_PARAMS = {"rigidbody": {"sampling": 1000}, "clustfcc": {"min_population": 1}}
NCORES = core_count()           # every core granted to this job (affinity, cgroup quota, scheduler), not the whole node
TOP = 5                         # best models compared with the complex
CAPRI_COLUMNS = ["model", "score", "irmsd", "fnat", "lrmsd", "dockq", "cluster_id"]
OUTDIR = Path(__file__).resolve().parents[1] / "output" / f"docking_{NAME}_{METHOD}"
RUN, TABLES, FIGURES, STRUCTURE = (OUTDIR / name for name in ("run", "tables", "figures", "structure"))

###################################### PREPARE

barnase = SyncedStructure.from_rcsb(APO_BARNASE)
barstar = SyncedStructure.from_rcsb(APO_BARSTAR)
holo = SyncedStructure.from_rcsb(HOLO)
for folder in (TABLES, FIGURES, STRUCTURE):
    folder.mkdir(parents=True, exist_ok=True)

dock = HaddockDocking(
    [(barnase, "A"), (barstar, "A")],   # one (structure, chain or list of chains) per partner
    workdir=RUN,
    method=METHOD,                      # 'ab-initio' (default), 'oracle', 'targeted' or 'alascan'
    seed=SEED,                          # default 917
    ncores=NCORES,                      # default None: core_count()
    separate="auto",                    # move partners apart when they share a structure or touch
    clearance_a=10.0,                   # smallest gap in Å between separated partners
    module_params=MODULE_PARAMS,        # default None: the preset
    global_params=None,                 # top-level haddock3 parameters, such as {"max_cpus": True}
    reference=None,                     # default 'auto'; the complex needs a chain map, so set_reference() below
    caprieval=True,                     # end with caprieval against the reference
    haddock3_exe=None,                  # search PATH and this Python's bin directory
)
dock.set_reference(holo, chain_map=CHAIN_MAP)   # caprieval reference only: no restraints are taken from it
# no restraints: adding ambiguous ones (add_ambiguous_restraints, add_contact_restraints) turns ranair off
dock.export_inputs()
dock.validate_config()                  # haddock3 --setup in a scratch folder: catches unknown module_params

###################################### RUN

dock.run_haddock(overwrite=True, progress=True)

###################################### TABLES (TSV)

models = dock.load_docked_models(top=TOP, step=None)   # barnase is chain A, barstar comes back as AB
interfaces = pd.DataFrame(
    [{"model": HOLO, **holo.interface_complex("A", "D")}]
    + [{"model": name, **model.interface_complex("A", "AB")} for name, model in models.items()])

tables = {
    "residue_map": dock.residue_table(),
    "scores": dock.score_table(step=None),
    "docked_models": dock.docked_model_table(step=None),   # final models with irmsd, fnat, lrmsd, dockq vs 1BRS
    "interfaces": interfaces,
}
capri = dock.caprieval_tables(step=None)
tables.update({f"caprieval_{name}": table for name, table in capri.items()})
table_paths = []
for name, table in tables.items():
    path = TABLES / f"{NAME}_{METHOD}_{name}.tsv"
    table.to_csv(path, sep="\t", index=False)
    table_paths.append(path)

###################################### FIGURES (PDF) AND 3D VIEW (HTML)

figure_paths = holo.export_plot(
    lambda: dock.plot_caprieval(x="irmsd", y="score", color="dockq", step=None,
                                title=f"Barnase–barstar vs {HOLO} ({METHOD})", ax=None),
    "caprieval", outdir=FIGURES)
best_name, best = next(iter(models.items()))
structure_paths = [best.export_structure_3d(outdir=STRUCTURE, title=f"Barnase–barstar, best docked model ({best_name})")]

print(f"caprieval of the top {TOP} models against {HOLO}:")
print(capri["models"].head(TOP)[[column for column in CAPRI_COLUMNS if column in capri["models"]]].to_string(index=False))
for folder, paths in ((TABLES, table_paths), (FIGURES, figure_paths), (STRUCTURE, structure_paths)):
    print(f"{folder.relative_to(OUTDIR.parent)}/")
    for path in paths:
        print(f"  {Path(path).name}")
