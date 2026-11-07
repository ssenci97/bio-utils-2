"""Dock unbound barnase (1A2P) and unbound barstar (1A19) with HADDOCK3 and evaluate the models against the complex (1BRS).

Needs the haddock extra and the haddock3 executable on PATH. Run from anywhere with:  python examples/dock_barnase_barstar.py

The partners come from two different PDB entries. The bound complex 1BRS gives both the restraints (its interface
contacts, transferred onto the apo chains through reference positions; 1A19 barstar has an extra Met and Cys40, so
its sequence is aligned to 1BRS's first) and the caprieval reference, written in the partners' HADDOCK3 numbering. The
'oracle' preset is used, downsized to 40 rigid-body models for a laptop: clustering, refinement of the best clusters, caprieval against 1BRS.

Outputs land in examples/output/docking_barnase_barstar/:
    run/        the HADDOCK3 inputs, config.toml and run folder
    tables/     restraints, residue map, model scores, caprieval models and clusters vs 1BRS, and the interface of
                1BRS vs the top models (TSV)
    figures/    HADDOCK score vs i-RMSD colored by DockQ (PDF)
    structure/  interactive 3D view of the best model as standalone HTML
"""

import os
from pathlib import Path

import pandas as pd

from synced_structure import HaddockDocking, SyncedStructure

NAME = "barnase_barstar"
APO_BARNASE, APO_BARSTAR, HOLO = "1A2P", "1A19", "1BRS"
# partner chains (partner number, chain) → chains of the complex: barnase A–A, barstar A–D
CHAIN_MAP = {(1, "A"): "A", (2, "A"): "D"}
METHOD = "oracle"
SEED = 11
# laptop-sized run: 40 rigid-body models instead of 200, at most 4 cores; set MODULE_PARAMS = {} for the full preset
MODULE_PARAMS = {"rigidbody": {"sampling": 40}, "clustfcc": {"min_population": 1}}
NCORES = min(4, os.cpu_count() or 2)
TOP = 5                         # best models compared with the complex
CAPRI_COLUMNS = ["model", "score", "irmsd", "fnat", "lrmsd", "dockq", "cluster_id"]
OUTDIR = Path(__file__).resolve().parent / "output" / f"docking_{NAME}"
RUN, TABLES, FIGURES, STRUCTURE = (OUTDIR / name for name in ("run", "tables", "figures", "structure"))

###################################### PREPARE

barnase = SyncedStructure.from_rcsb(APO_BARNASE)
barstar = SyncedStructure.from_rcsb(APO_BARSTAR)
holo = SyncedStructure.from_rcsb(HOLO)
for folder in (TABLES, FIGURES, STRUCTURE):
    folder.mkdir(parents=True, exist_ok=True)

dock = HaddockDocking([(barnase, "A"), (barstar, "A")], workdir=RUN, method=METHOD, seed=SEED,
                      ncores=NCORES, module_params=MODULE_PARAMS)
dock.set_reference(holo, chain_map=CHAIN_MAP)           # caprieval reference and contact source
dock.add_contact_restraints(kind="ambiguous")           # interface residues of 1BRS become active
dock.export_inputs()

###################################### RUN

dock.run_haddock(overwrite=True)

###################################### TABLES (TSV)

models = dock.load_docked_models(top=TOP)               # barnase is chain A, barstar comes back as AB
interfaces = pd.DataFrame(
    [{"model": HOLO, **holo.interface_complex("A", "D")}]
    + [{"model": name, **model.interface_complex("A", "AB")} for name, model in models.items()])

tables = {
    "restraints": dock.restraint_table(),
    "residue_map": dock.residue_table(),
    "scores": dock.score_table(),
    "docked_models": dock.docked_model_table(),       # final models with irmsd, fnat, lrmsd, dockq vs 1BRS
    "interfaces": interfaces,
}
capri = dock.caprieval_tables()
tables.update({f"caprieval_{name}": table for name, table in capri.items()})
table_paths = []
for name, table in tables.items():
    path = TABLES / f"{NAME}_{name}.tsv"
    table.to_csv(path, sep="\t", index=False)
    table_paths.append(path)

###################################### FIGURES (PDF) AND 3D VIEW (HTML)

figure_paths = holo.export_plot(lambda: dock.plot_caprieval(title=f"Barnase–barstar vs {HOLO}"),
                                "caprieval", outdir=FIGURES)
best_name, best = next(iter(models.items()))
structure_paths = [best.export_structure_3d(outdir=STRUCTURE, title=f"Barnase–barstar, best docked model ({best_name})")]

print(f"caprieval of the top {TOP} models against {HOLO}:")
print(capri["models"].head(TOP)[[column for column in CAPRI_COLUMNS if column in capri["models"]]].to_string(index=False))
for folder, paths in ((TABLES, table_paths), (FIGURES, figure_paths), (STRUCTURE, structure_paths)):
    print(f"{folder.relative_to(OUTDIR.parent)}/")
    for path in paths:
        print(f"  {Path(path).name}")
