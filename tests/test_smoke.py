"""Offline smoke tests on 1AAY: zif268 zinc fingers (protein) bound to two DNA strands, with Zn ions and water."""

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from synced_structure import SyncedStructure, build_panel_grid, export_panels_pdf, side_by_side  # noqa: E402

CIF = Path(__file__).parent / "data" / "1AAY.cif"


@pytest.fixture
def sync():
    """Return a fresh SyncedStructure of 1AAY keyed by author chain IDs (A protein, B and C DNA)."""
    return SyncedStructure.from_path(CIF)


def test_chains_and_moieties(sync):
    assert sync.chain_ids == ["B", "C", "A"]
    assert sync.get_sequence("B") == "AGCGTGGGCGT"
    assert set(sync.moiety_summary["moiety"]) == {"protein", "nucleic", "ion", "water"}
    assert sync.chains_by_moiety("nucleic") == ["A", "B"]
    assert sync.summary()["mapping_fraction"].eq(1.0).all()
    assert sync.get_external_ref("B") == "1AAY"


def test_indices_are_one_based_and_consistent(sync):
    residues = sync.residue_map
    contacts = sync.compute_contacts(contact_type="all")
    assert residues["sequence_index"].min() == 1
    protein = residues[residues["chain_id"] == "C"]
    assert protein["reference_index"].iloc[0] == 3
    merged = contacts.merge(residues, left_on=["chain_1", "sequence_index_1"], right_on=["chain_id", "sequence_index"])
    assert (merged["resnum_1"] == merged["pdb_resnum"]).all()


def test_contacts_and_strings(sync, tmp_path):
    contacts = sync.compute_contacts(contact_type="inter", moieties=("protein", "nucleic", "ion"))
    assert {"protein", "nucleic"} <= set(contacts["moiety_1"]) | set(contacts["moiety_2"])
    assert contacts.loc[contacts["moiety_2"] == "ion", "sequence_index_2"].isna().all()
    assert "1" in sync.get_contact_sequence("B")
    path = sync.export_contact_sequences(tmp_path / "out" / "strings.csv", reference=True)
    assert path.is_file()
    assert sync.chain_contact_partners() == {"A": ["B", "C"], "B": ["A"], "C": ["A"]}


def test_matrices_and_statistics(sync):
    pair = sync.contact_matrix("A", "B", scheme="sidechain")
    assert pair["matrix"].shape == (85, 11)
    block = sync.block_contact_matrix("A", "B")
    assert block["inter"].any()
    stats = sync.contact_statistics()
    assert {"pairs", "profile", "separation", "moiety_pairs"} == set(stats)
    assert "nucleic/protein" in set(stats["moiety_pairs"]["moiety_pair"])


def test_interfaces(sync):
    complex_row = sync.interface_complex("A", "B")
    assert complex_row["dASA_A2"] > 0
    assert len(sync.chain_partner_frame("A")) == 2
    residues = sync.interface_residues("A", "B")
    assert (residues["energy_per_partner"].abs() <= residues["contact_energy"].abs()).all()


def test_reference_numbering_keeps_hetero(sync):
    before = sync.moiety_summary
    sync.apply_reference_numbering()
    after = sync.moiety_summary
    hetero = ["ion", "water"]
    assert before[before["moiety"].isin(hetero)].equals(after[after["moiety"].isin(hetero)])
    assert sync["A"].residues[0].id[1] == 3


def test_plots_and_exports(sync, tmp_path):
    fig, _ = sync.plot_contact_map("A", "B")
    for chains in (("A", "B"), ("A", "B", "C")):
        sync.plot_contacts_multi(chains, plddt=True, show_color_legend=True)
    sync.plot_contact_panel("A", "B")
    with pytest.raises(ValueError):
        sync.plot_contacts_multi(("A", "A"))
    paths = sync.export_plot(lambda: sync.plot_contact_map("A"), "map", outdir=tmp_path, formats=("png", "pdf"))
    assert all(p.is_file() for p in paths)
    pdf = export_panels_pdf([([lambda ax: sync.plot_contact_map("A", ax=ax)] * 2, "page")], tmp_path)
    assert pdf.is_file()
    assert build_panel_grid([lambda ax: ax.plot([0, 1])]).axes
    for export in (sync.export_chains, sync.export_residues, sync.export_moieties, sync.export_alignment_stats):
        assert export(outdir=tmp_path).parent == tmp_path


def test_viz_and_config(sync, tmp_path):
    views = [sync.show_structure_3d(title="all"), sync.show_structure_3d(style="stick")]
    assert views[0].data != views[1].data
    assert "<div" in side_by_side(*views).data
    config = tmp_path / "config.yaml"
    config.write_text("synced_structure:\n  save_outdir: elsewhere\n")
    loaded = SyncedStructure.from_path(CIF, config_path=config)
    assert loaded.config["save_outdir"] == "elsewhere"
    assert config.read_text() == "synced_structure:\n  save_outdir: elsewhere\n"
