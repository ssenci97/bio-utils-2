"""Download the STAT6-DNA complex 5D39 and export every table (TSV), figure (PDF, PNG, SVG), and 3D view (HTML).

Run from anywhere with:  python examples/export_all.py
Outputs land in examples/output/.
"""

from pathlib import Path

from synced_structure import SyncedStructure, export_panels_pdf

PDB_ID = "5D39"
TITLE = "STAT6-DNA complex 5D39"
OUTDIR = Path(__file__).resolve().parent / "output"
FORMATS = ("pdf", "png", "svg")

###################################### LOAD

s = SyncedStructure.from_rcsb(PDB_ID)
OUTDIR.mkdir(parents=True, exist_ok=True)
paths = []

summary_path = OUTDIR / f"{PDB_ID}_summary.txt"
summary_path.write_text(s.summary().to_string() + "\n")
paths.append(summary_path)

# Interacting chain pairs, each listed once
partners = s.chain_contact_partners()
pairs = sorted({tuple(sorted((a, b))) for a, chains in partners.items() for b in chains})

###################################### TABLES (TSV)

paths += [
    s.export_chains(outdir=OUTDIR),
    s.export_residue_mapping(outdir=OUTDIR),
    s.export_residues(outdir=OUTDIR),
    s.export_moieties(outdir=OUTDIR),
    s.export_alignment_stats(outdir=OUTDIR),
    s.export_alignment_positions(outdir=OUTDIR),
    s.export_interface_complex(pairs, outdir=OUTDIR),
    s.export_interface_residues(pairs, outdir=OUTDIR),
    s.export_contact_sequences(OUTDIR / f"{PDB_ID}_contact_sequences.tsv"),
    s.export_contact_sequences(OUTDIR / f"{PDB_ID}_contact_sequences_reference.tsv", reference=True),
]
for key, frame in s.contact_statistics().items():
    path = OUTDIR / f"{PDB_ID}_contact_stats_{key}.tsv"
    frame.to_csv(path, sep="\t", index=False)
    paths.append(path)

###################################### FIGURES (PDF, PNG, SVG)

figures = {
    "contacts_multi_A_C": lambda: s.plot_contacts_multi(("A", "C"), show_color_legend=True),
    "contact_map_N_M_sidechain": lambda: s.plot_contact_map("N", "M", scheme="sidechain", n_ticks=None, tick_step=5),
    "contact_map_A": lambda: s.plot_contact_map("A"),
    "contact_panel_A_N": lambda: s.plot_contact_panel("A", "N"),
}
for name, plot in figures.items():
    paths += s.export_plot(plot, name, outdir=OUTDIR, formats=FORMATS)

# One multi-page PDF: a contact map per interacting pair, four per page
maps = [lambda ax, a=a, b=b: s.plot_contact_map(a, b, ax=ax) for a, b in pairs]
pages = [(maps[i:i + 4], f"{TITLE}: inter-chain contacts") for i in range(0, len(maps), 4)]
paths.append(export_panels_pdf(pages, OUTDIR, filename=f"{PDB_ID}_panels.pdf"))

###################################### 3D VIEWS (HTML)

paths.append(s.export_structure_3d(outdir=OUTDIR, title=TITLE))
paths.append(s.export_structure_3d("structure_3d_stick", outdir=OUTDIR, title=f"{TITLE} (sticks)", style="stick"))

print(f"Wrote {len(paths)} files to {OUTDIR}:")
for path in paths:
    print(f"  {Path(path).name}")
