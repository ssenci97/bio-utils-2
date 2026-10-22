# 5. Export

**Goal:** write tables (TSV), figures (PDF/PNG) and an interactive 3D page to a folder, ready to share.

**You need:** [guide 1](01-getting-started.md). For a complete script that exports everything for one entry, see [`examples/export_all.py`](../../examples/export_all.py).

## Step 1: tables

Every `export_*` method writes `{structure_id}_{name}.tsv` into `outdir` (created if missing) and returns the path. Without `outdir`, files go to `data/saved` (the `save_outdir` setting).

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")
out = "results"

tables = [
    s.export_chain_summary(outdir=out),                      # 1AAY_chains.tsv
    s.export_residue_mapping(outdir=out),                    # 1AAY_residue_mapping.tsv
    s.export_residue_table(outdir=out),                      # 1AAY_residues.tsv
    s.export_moiety_table(outdir=out),                       # 1AAY_moieties.tsv
    s.export_reference_alignment_stats(outdir=out),          # 1AAY_alignment_stats.tsv
    s.export_interface_complex([("A", "B"), ("A", "C")], outdir=out),   # 1AAY_interface_complex.tsv
    s.export_interface_residues([("A", "B")], outdir=out),              # 1AAY_interface_residues.tsv
    s.export_contact_strings(f"{out}/1AAY_contact_strings.tsv", reference=True),
]
print([path.name for path in tables])
```

`export_contact_strings` takes a file path rather than a folder, because its name depends on your cutoff and scheme. Any table you compute yourself can be saved the same way with pandas: `s.compute_contacts().to_csv("contacts.tsv", sep="\t", index=False)`.

## Step 2: one figure

`export_plot` takes a function that *returns* `(figure, axes)`, so the figure is drawn, saved in every format you ask for, and closed:

```python
paths = s.export_plot(lambda: s.plot_contact_map("A", "B", crop=True), "map_A_B", outdir=out, formats=("pdf", "png"))
print([path.name for path in paths])
```

```text
['1AAY_map_A_B.pdf', '1AAY_map_A_B.png']
```

Figures are saved at exactly the size they were drawn at, with editable text in PDFs.

## Step 3: several panels on one page

Panel functions take an `ax` to draw on. `export_plot_grid` lays them out lettered A, B, …:

```python
panels = [
    lambda ax: s.plot_contact_map("A", "B", ax=ax),
    lambda ax: s.plot_contact_map("A", "C", ax=ax),
]
s.export_plot_grid(panels, "maps_grid", outdir=out)          # results/maps_grid.pdf
```

Note that the grid file is named exactly `maps_grid.pdf`, without the structure ID. To put several grids into one multi-page PDF, use `export_panels_pdf`:

```python
from synced_structure import export_panels_pdf

export_panels_pdf([(panels, "Zif268 on both DNA strands")], out, filename="panels.pdf")
```

## Step 4: an interactive 3D page

```python
page = s.export_structure_3d(outdir=out, title="Zif268 on DNA")   # results/1AAY_structure_3d.html
print(page.name)
```

```text
1AAY_structure_3d.html
```

Open the HTML file in any browser; it needs internet access to load 3Dmol.js. Extra arguments (`style=`, `colors=`, …) are passed to `show_structure_3d`.

## Common mistakes

- **Passing a figure instead of a function.** `export_plot(fig, ...)` fails; wrap the call: `export_plot(lambda: s.plot_contact_map(...), ...)`.
- **Panel functions without `ax`.** Grid panels are called as `func(ax=ax)`; write `lambda ax: s.plot_contact_map("A", "B", ax=ax)`.
- **Looking for files in the wrong folder.** Without `outdir`, exports go to `data/saved` relative to where Python runs.

## Next

- [6. Docking](06-docking.md)
- [API index](../api_index.md) for every export method.
