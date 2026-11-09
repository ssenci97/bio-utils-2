"""Compare structures: chain contact graph, RMSD, CAPRI metrics and a 3D superposition, on barnase-barstar.

1BRS holds three barnase-barstar complexes (barnase A, B, C; barstar D, E, F); 1A2P is barnase from another entry.

Run from anywhere with:  python examples/compare_structures.py

Outputs land in examples/output/compare/:
    1BRS_chain_graph.pdf          chains as nodes, edges as thick as their residue contacts
    1BRS_chain_graph.tsv          the contact counts behind the graph
    1A2P_on_1BRS_residues.tsv     residues of 1A2P matched to 1BRS chain A
    1BRS_capri.tsv                each barnase-barstar copy scored against the A-D complex
    1A2P_on_1BRS_3d.html          1A2P superposed on 1BRS, in one interactive 3D view
"""

from pathlib import Path

import pandas as pd

from synced_structure import SyncedStructure

OUTDIR = Path(__file__).resolve().parent / "output" / "compare"
OUTDIR.mkdir(parents=True, exist_ok=True)

complex_ = SyncedStructure.from_rcsb("1BRS")
barnase = SyncedStructure.from_rcsb("1A2P")

###################################### CHAIN CONTACT GRAPH

graph = complex_.chain_graph_table()
graph.to_csv(OUTDIR / "1BRS_chain_graph.tsv", sep="\t", index=False)
paths = complex_.export_plot(complex_.plot_chain_graph, "chain_graph", outdir=OUTDIR)

###################################### RMSD OF ANOTHER ENTRY

matches = barnase.residue_match_table(complex_, chain_map={"A": "A"})
matches.to_csv(OUTDIR / "1A2P_on_1BRS_residues.tsv", sep="\t", index=False)
for atoms in ("ca", "backbone", "heavy"):
    rmsd = barnase.rmsd_metrics(complex_, chain_map={"A": "A"}, atoms=atoms)
    print(f"1A2P A on 1BRS A, {atoms:<8}: RMSD {rmsd['rmsd_A']:.2f} Å over {rmsd['n_atoms']} atoms")

###################################### CAPRI METRICS OF EACH COPY AGAINST A-D

copies = {f"{r}-{l}": complex_.subset_chains([r, l], structure_id=f"1BRS_{r}{l}") for r, l in (("B", "E"), ("C", "F"))}
capri = [complex_.capri_table({name: model}, receptor=name[0], ligand=name[2],
                              chain_map={name[0]: "A", name[2]: "D"}) for name, model in copies.items()]
table = pd.concat(capri, ignore_index=True)
table.to_csv(OUTDIR / "1BRS_capri.tsv", sep="\t", index=False)
print(table[["structure_id", "irmsd_A", "lrmsd_A", "fnat", "dockq", "capri_class"]].to_string(index=False))

###################################### 3D SUPERPOSITION (HTML)

view = barnase.show_superposition_3d(complex_, chain_map={"A": "A"})
page = OUTDIR / "1A2P_on_1BRS_3d.html"
page.write_text(f'<!doctype html><html><head><meta charset="utf-8"><title>1A2P on 1BRS</title></head>'
                f'<body style="margin:16px">{view.data}</body></html>', encoding="utf-8")

for path in [*paths, page]:
    print(Path(path).relative_to(OUTDIR.parent.parent))
