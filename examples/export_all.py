"""Download adenylyl cyclase 9 bound to Gαs (6R3Q) and export every table (TSV), contact maps (PDF), and a 3D view (HTML).

Run from anywhere with:  python examples/export_all.py

Outputs land in examples/output/:
    tables/     every table as TSV
    figures/    contact maps as PDF, each page exactly the size it was drawn at
    structure/  interactive 3D view as standalone HTML
"""

from pathlib import Path

from synced_structure import SyncedStructure, export_panels_pdf

PDB_ID = "6R3Q"
OUTDIR = Path(__file__).resolve().parent / "output"
TABLES, FIGURES, STRUCTURE = (OUTDIR / name for name in ("tables", "figures", "structure"))
CUTOFF_NM = 0.8                 # C-alpha contacts drawn in color
SIDECHAIN_CUTOFF_NM = 0.45      # closest heavy-atom contacts drawn in color
GRAY_FACTOR = 2.0               # residue pairs up to this many cutoffs drawn in gray behind the contacts
GRID_CUTOFFS_NM = (0.6, 0.8, 1.0)

###################################### LOAD

s = SyncedStructure.from_rcsb(PDB_ID)
for folder in (TABLES, FIGURES, STRUCTURE):
    folder.mkdir(parents=True, exist_ok=True)

# Contact statistics report label chain IDs; plots take chain keys (author IDs by default)
key_of = {record.chain_id: key for key, record in s.chains.items()}
summary = s.summary()
stats = s.contact_statistics()
inter = stats["pairs"].query("contact_type == 'inter'")
inter_pairs = [(key_of[a], key_of[b]) for a, b in zip(inter["chain_a"], inter["chain_b"])]
polymers = [key_of[row.chain_id] for row in summary.itertuples() if row.moiety in ("protein", "nucleic")]

###################################### TABLES (TSV)

summary_path = TABLES / f"{PDB_ID}_summary.tsv"
summary.to_csv(summary_path, sep="\t", index=False)
tables = [
    summary_path,
    s.export_chains(outdir=TABLES),
    s.export_residue_mapping(outdir=TABLES),
    s.export_residues(outdir=TABLES),
    s.export_moieties(outdir=TABLES),
    s.export_alignment_stats(outdir=TABLES),
    s.export_alignment_positions(outdir=TABLES),
    s.export_interface_complex(inter_pairs, outdir=TABLES),
    s.export_interface_residues(inter_pairs, outdir=TABLES),
    s.export_contact_sequences(TABLES / f"{PDB_ID}_contact_sequences.tsv"),
    s.export_contact_sequences(TABLES / f"{PDB_ID}_contact_sequences_reference.tsv", reference=True),
]
for key, frame in stats.items():
    path = TABLES / f"{PDB_ID}_contact_stats_{key}.tsv"
    frame.to_csv(path, sep="\t", index=False)
    tables.append(path)

###################################### FIGURES (PDF)

chains = "_".join(polymers)
figures = {
    # all chains on one map: intra-chain contacts in gray, the interface in color
    f"contacts_{chains}": lambda: s.plot_contacts_multi(
        polymers, cutoff=CUTOFF_NM, gray_cutoff=CUTOFF_NM * GRAY_FACTOR, size_in=3.5,
        show_color_legend=True, show_distance_legend=True),
    # same map at side-chain resolution
    f"contacts_{chains}_sidechain": lambda: s.plot_contacts_multi(
        polymers, cutoff=SIDECHAIN_CUTOFF_NM, gray_cutoff=SIDECHAIN_CUTOFF_NM * GRAY_FACTOR, scheme="sidechain",
        size_in=3.5, show_color_legend=True, show_distance_legend=True),
}
figure_paths = []
for name, plot in figures.items():
    figure_paths += s.export_plot(plot, name, outdir=FIGURES)

# How the interface grows with the cutoff, one panel per cutoff
maps = [lambda ax, c=c: s.plot_contacts_multi(polymers, cutoff=c, gray_cutoff=c * GRAY_FACTOR, size_in=3.5, ax=ax)
        for c in GRID_CUTOFFS_NM]
figure_paths.append(export_panels_pdf([(maps, f"{PDB_ID} contacts by cutoff")], FIGURES,
                                      filename=f"{PDB_ID}_contacts_by_cutoff.pdf",
                                      sizes=3.5, ncols=len(GRID_CUTOFFS_NM)))

###################################### 3D VIEW (HTML)

structure_paths = [s.export_structure_3d(outdir=STRUCTURE, title=f"Adenylyl cyclase 9 – Gαs ({PDB_ID})")]

for folder, paths in ((TABLES, tables), (FIGURES, figure_paths), (STRUCTURE, structure_paths)):
    print(f"{folder.relative_to(OUTDIR.parent)}/")
    for path in paths:
        print(f"  {Path(path).name}")
