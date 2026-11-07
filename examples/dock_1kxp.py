"""Re-dock actin and vitamin D-binding protein (1KXP) with HADDOCK3 and compare the models with the crystal.

Needs the haddock extra and the haddock3 executable on PATH. Run from anywhere with:  python examples/dock_1kxp.py

The partners are pulled apart with random rotations, then docked back with ambiguous restraints taken from the
crystal contacts (the 'oracle' preset: 200 rigid-body models, clustering, refinement of the best clusters only,
caprieval against the crystal). Progress is logged while haddock3 runs.

Outputs land in examples/output/docking_1KXP/:
    run/        the HADDOCK3 inputs, config.toml and run folder
    tables/     restraints, model scores, caprieval models and clusters vs the crystal, and the interface of the
                crystal vs the top models (TSV)
    structure/  interactive 3D view of the best model as standalone HTML
"""

import os
from pathlib import Path

import pandas as pd

from synced_structure import SyncedStructure

PDB_ID = "1KXP"
ACTIN, DBP = "A", "D"
METHOD = "oracle"
SEED = 7
TOP = 5                         # best models compared with the crystal
CAPRI_COLUMNS = ["model", "score", "irmsd", "fnat", "lrmsd", "dockq", "cluster_id"]
OUTDIR = Path(__file__).resolve().parent / "output" / f"docking_{PDB_ID}"
RUN, TABLES, STRUCTURE = (OUTDIR / name for name in ("run", "tables", "structure"))

###################################### PREPARE

s = SyncedStructure.from_rcsb(PDB_ID)
for folder in (TABLES, STRUCTURE):
    folder.mkdir(parents=True, exist_ok=True)

dock = s.prepare_haddock_docking(ACTIN, DBP, workdir=RUN, method=METHOD, seed=SEED, ncores=os.cpu_count() or 4)
dock.add_contact_restraints(kind="ambiguous")           # interface residues of the crystal become active
dock.export_inputs()

###################################### RUN

dock.run_haddock(overwrite=True)

###################################### TABLES (TSV)

models = dock.load_docked_models(top=TOP)
interfaces = pd.DataFrame(
    [{"model": PDB_ID, **s.interface_complex(ACTIN, DBP)}]
    + [{"model": name, **model.interface_complex(ACTIN, DBP)} for name, model in models.items()])

tables = {
    "restraints": dock.restraint_table(),
    "residue_map": dock.residue_table(),
    "scores": dock.score_table(),
    "docked_models": dock.docked_model_table(),       # final models with irmsd, fnat, lrmsd, dockq vs the crystal
    "interfaces": interfaces,
}
capri = dock.caprieval_tables()                         # {"models": per model, "clusters": per cluster} vs the crystal
tables.update({f"caprieval_{name}": table for name, table in capri.items()})
table_paths = []
for name, table in tables.items():
    path = TABLES / f"{PDB_ID}_{name}.tsv"
    table.to_csv(path, sep="\t", index=False)
    table_paths.append(path)

###################################### 3D VIEW (HTML)

best_name, best = next(iter(models.items()))
structure_paths = [best.export_structure_3d(outdir=STRUCTURE, title=f"Actin – DBP, best docked model ({best_name})")]

print(f"caprieval of the top {TOP} models against the crystal:")
print(capri["models"].head(TOP)[[column for column in CAPRI_COLUMNS if column in capri["models"]]].to_string(index=False))
for folder, paths in ((TABLES, table_paths), (STRUCTURE, structure_paths)):
    print(f"{folder.relative_to(OUTDIR.parent)}/")
    for path in paths:
        print(f"  {Path(path).name}")
