"""Re-dock the zinc-finger protein of 1AAY against its two DNA strands, which form one multi-chain partner.

Needs the haddock extra and the haddock3 executable on PATH. Run from anywhere with:
    python examples/docking/1aay_two_chains_oracle.py

Partner 1 is the protein (chain A); partner 2 is the DNA duplex made of two chains (B and C) that move together as
one rigid body, so HADDOCK3 sees two molecules while the structure has three chains. Restraints come from the
crystal contacts. The run is downsized to 20 rigid-body models; set MODULE_PARAMS = {} for the full preset (200).

Outputs land in examples/output/docking_1AAY_oracle/:
    run/        the HADDOCK3 inputs, config.toml and run folder
    tables/     residue map, restraints, model scores and caprieval vs the crystal (TSV)
    figures/    HADDOCK score vs i-RMSD colored by DockQ (PDF)
    structure/  3D view of the best model with the restrained residues highlighted (HTML)
"""

from pathlib import Path

from synced_structure import SyncedStructure, core_count

PDB_ID = "1AAY"
PROTEIN, DUPLEX = "A", ["B", "C"]       # one partner per argument: a chain or a list of chains
METHOD = "oracle"
SEED = 7
# laptop-sized run: 20 rigid-body models instead of 200; {} keeps the preset
MODULE_PARAMS = {"rigidbody": {"sampling": 20}, "clustfcc": {"min_population": 1}}
NCORES = core_count()
TOP = 5
DIM_ALPHA = 0.25                        # opacity of the residues that are not highlighted
OUTDIR = Path(__file__).resolve().parents[1] / "output" / f"docking_{PDB_ID}_{METHOD}"
RUN, TABLES, FIGURES, STRUCTURE = (OUTDIR / name for name in ("run", "tables", "figures", "structure"))

###################################### PREPARE

s = SyncedStructure.from_rcsb(PDB_ID)
for folder in (TABLES, FIGURES, STRUCTURE):
    folder.mkdir(parents=True, exist_ok=True)

dock = s.prepare_haddock_docking(
    PROTEIN, DUPLEX,
    workdir=RUN,
    method=METHOD,
    seed=SEED,
    ncores=NCORES,
    module_params=MODULE_PARAMS,
)
dock.add_contact_restraints(kind="ambiguous", cutoff=0.8, scheme="ca")
dock.export_inputs()
dock.validate_config()

###################################### RUN

dock.run_haddock(overwrite=True, progress=True)

###################################### TABLES (TSV)

tables = {
    "residue_map": dock.residue_table(),    # partner and segid of every chain, HADDOCK3 numbering vs yours
    "restraints": dock.restraint_table(),
    "scores": dock.score_table(step=None),
    "docked_models": dock.docked_model_table(step=None),
}
tables.update({f"caprieval_{name}": table for name, table in dock.caprieval_tables(step=None).items()})
table_paths = []
for name, table in tables.items():
    path = TABLES / f"{PDB_ID}_{METHOD}_{name}.tsv"
    table.to_csv(path, sep="\t", index=False)
    table_paths.append(path)

###################################### FIGURES (PDF) AND 3D VIEW (HTML)

figure_paths = s.export_plot(
    lambda: dock.plot_caprieval(x="irmsd", y="score", color="dockq", step=None,
                                title=f"Zif268 – DNA duplex vs {PDB_ID} ({METHOD})", ax=None),
    "caprieval", outdir=FIGURES)

models = dock.load_docked_models(top=TOP, step=None)
best_name, best = next(iter(models.items()))
restraints, residues = dock.restraint_table(), dock.residue_table()
restrained = {(row[f"segid_{side}"], resid) for _, row in restraints.iterrows()
              for side in "ab" for resid in row[f"residues_{side}"]}
mapped = residues[[(segid, resid) in restrained for segid, resid in zip(residues["segid"], residues["haddock_resid"])]]
highlight = {chain: group["pdb_resnum"].tolist() for chain, group in mapped.groupby("chain_key")}
structure_paths = [best.export_structure_3d(
    outdir=STRUCTURE, title=f"Best docked model ({best_name}), restrained residues highlighted",
    highlight=highlight,                # {chain: residue numbers or "10-25,40"}; None draws everything normally
    dim_alpha=DIM_ALPHA,                # opacity of the other residues
    highlight_sticks=True)]             # also show the side chains of the highlighted residues

print(dock.residue_table().groupby(["segid", "chain_key"]).size().to_string())
for folder, paths in ((TABLES, table_paths), (FIGURES, figure_paths), (STRUCTURE, structure_paths)):
    print(f"{folder.relative_to(OUTDIR.parent)}/")
    for path in paths:
        print(f"  {Path(path).name}")
