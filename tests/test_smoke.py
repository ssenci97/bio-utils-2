"""Offline smoke tests on 1AAY: zif268 zinc fingers (protein) bound to two DNA strands, with Zn ions and water."""

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from synced_structure import SyncedStructure, build_panel_grid, export_panels_pdf, show_side_by_side_3d  # noqa: E402

CIF = Path(__file__).parent / "data" / "1AAY.cif"


@pytest.fixture
def sync():
    """Return a fresh SyncedStructure of 1AAY keyed by author chain IDs (A protein, B and C DNA)."""
    return SyncedStructure.from_path(CIF)


def test_chains_and_moieties(sync):
    assert sync.chain_ids == ["B", "C", "A"]
    assert sync.get_sequence("B") == "AGCGTGGGCGT"
    assert set(sync.moiety_table["moiety"]) == {"protein", "nucleic", "other"}
    assert sync.chain_ids_by_moiety("nucleic") == ["A", "B"]
    summary = sync.chain_table()
    assert summary.columns[0] == "structure_id" and summary["structure_id"].eq(sync.structure_id).all()
    assert summary["mapping_fraction"].eq(1.0).all()
    assert dict(zip(summary["chain_id"], summary["moiety"])) == {"B": "nucleic", "C": "protein", "A": "nucleic"}
    polymer = sync.moiety_table.query("moiety != 'other'").groupby("chain_id")["moiety"].agg("+".join)
    assert dict(zip(summary["chain_id"], summary["moiety"])) == polymer.to_dict()
    assert "partners" not in summary
    with_partners = sync.chain_table(partners=True)
    assert dict(zip(with_partners["auth_chain_id"], with_partners["partners"])) == {"A": "B;C", "B": "A", "C": "A"}
    assert sync.get_accession("B") == "1AAY"


def test_indices_are_one_based_and_consistent(sync):
    residues = sync.residue_table
    contacts = sync.contact_table(contact_type="all")
    assert residues["sequence_index"].min() == 1
    protein = residues[residues["chain_id"] == "C"]
    assert protein["reference_index"].iloc[0] == 3
    merged = contacts.merge(residues, left_on=["chain_1", "sequence_index_1"], right_on=["chain_id", "sequence_index"])
    assert (merged["resnum_1"] == merged["pdb_resnum"]).all()


def test_contacts_and_strings(sync, tmp_path):
    contacts = sync.contact_table(contact_type="inter", moieties=("protein", "nucleic", "other"))
    assert {"protein", "nucleic"} <= set(contacts["moiety_1"]) | set(contacts["moiety_2"])
    assert contacts.loc[contacts["moiety_2"] == "other", "sequence_index_2"].isna().all()
    text = sync.contact_string("B")
    assert "1" in text
    masked = sync.masked_sequence("B")
    assert len(masked) == len(text)
    assert all((char != "-") == (bit == "1") for char, bit in zip(masked, text))
    assert sync.masked_sequence("B", mask=".").replace(".", "") == masked.replace("-", "")
    with pytest.raises(ValueError):
        sync.masked_sequence("B", mask="--")
    reference = sync.reference_masked_sequence("B")
    assert len(reference) == len(sync.get_reference_sequence("B"))
    assert all(char in ("-", ref) for char, ref in zip(reference, sync.get_reference_sequence("B")))
    path = sync.export_contact_strings(tmp_path / "out" / "strings.tsv", reference=True)
    assert path.is_file()
    assert "\t" in path.read_text().splitlines()[0]
    assert sync.chain_partner_map() == {"A": ["B", "C"], "B": ["A"], "C": ["A"]}


def test_matrices_and_statistics(sync):
    pair = sync.contact_matrix("A", "B", scheme="sidechain")
    assert pair["matrix"].shape == (85, 11)
    block = sync.block_contact_matrix("A", "B")
    assert block["inter"].any()
    stats = sync.contact_stat_tables()
    assert {"pairs", "profile", "separation", "moiety_pairs"} == set(stats)
    assert "nucleic/protein" in set(stats["moiety_pairs"]["moiety_pair"])


def test_interfaces(sync):
    complex_row = sync.interface_complex("A", "B")
    assert complex_row["dASA_A2"] > 0
    assert len(sync.chain_partner_table("A")) == 2
    residues = sync.interface_residue_table("A", "B")
    assert (residues["energy_per_partner"].abs() <= residues["contact_energy"].abs()).all()


def test_partner_structures(sync):
    sync.set_reference_sequence("B", "TT" + sync.get_sequence("B"))
    pairs = sync.chain_partner_structures("A")
    assert list(pairs) == ["1AAY_A_B", "1AAY_A_C"]
    child = pairs["1AAY_A_B"]
    assert child.structure_id == "1AAY_A_B"
    assert sorted(child.chain_ids) == ["A", "B"]
    assert sync.chain_ids == ["B", "C", "A"]
    assert child.get_reference_sequence("B") == sync.get_reference_sequence("B")
    assert child.mmcif_dict is sync.mmcif_dict
    child_row, parent_row = child.interface_complex("A", "B"), sync.interface_complex("A", "B")
    assert child_row.pop("structure") == "1AAY_A_B"
    parent_row.pop("structure")
    assert child_row == parent_row


def test_reference_numbering_keeps_hetero(sync):
    before = sync.moiety_table
    sync.apply_reference_numbering()
    after = sync.moiety_table
    assert before[before["moiety"] == "other"].equals(after[after["moiety"] == "other"])
    assert sync["A"].residues[0].id[1] == 3


def test_plots_and_exports(sync, tmp_path):
    fig, _ = sync.plot_contact_map("A", "B")
    for chains in (("A", "B"), ("A", "B", "C")):
        sync.plot_multi_contact_map(chains, plddt=True, show_color_legend=True)
    sync.plot_contact_map("A", "B", crop=True, value="distance", style="heatmap")
    with pytest.raises(ValueError):
        sync.plot_multi_contact_map(("A", "A"))
    paths = sync.export_plot(lambda: sync.plot_contact_map("A"), "map", outdir=tmp_path, formats=("png", "pdf"))
    assert all(p.is_file() for p in paths)
    pdf = export_panels_pdf([([lambda ax: sync.plot_contact_map("A", ax=ax)] * 2, "page")], tmp_path)
    assert pdf.is_file()
    assert build_panel_grid([lambda ax: ax.plot([0, 1])]).axes
    for export in (sync.export_chain_table, sync.export_residue_table, sync.export_moiety_table, sync.export_alignment_table):
        path = export(outdir=tmp_path)
        assert path.parent == tmp_path and path.suffix == ".tsv"
        assert "\t" in path.read_text().splitlines()[0]


def test_viz_and_config(sync, tmp_path):
    views = [sync.show_structure_3d(title="all"), sync.show_structure_3d(style="stick")]
    assert views[0].data != views[1].data
    assert "<div" in show_side_by_side_3d(*views).data
    page = sync.export_structure_3d(outdir=tmp_path, title="all").read_text()
    assert page.startswith("<!doctype html>") and "3Dmol" in page and "<title>all</title>" in page
    config = tmp_path / "config.yaml"
    config.write_text("synced_structure:\n  save_outdir: elsewhere\n")
    loaded = SyncedStructure.from_path(CIF, config_path=config)
    assert loaded.config["save_outdir"] == "elsewhere"
    assert config.read_text() == "synced_structure:\n  save_outdir: elsewhere\n"
