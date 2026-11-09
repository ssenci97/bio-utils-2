"""Offline tests of structure comparison (superposition, RMSD, CAPRI/DockQ) and chain contact graphs on 1BRS.

1BRS holds three barnase-barstar complexes: barnase A, B, C bound to barstar D, E, F.
"""

from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from synced_structure import SyncedStructure  # noqa: E402
from synced_structure.comparison import _kabsch, capri_class, dockq_score  # noqa: E402

DATA = Path(__file__).parent / "data"


@pytest.fixture
def ref():
    """Return 1BRS keyed by author chain IDs."""
    return SyncedStructure.from_path(DATA / "1BRS.cif")


def _rotation(axis, degrees):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    t = np.radians(degrees)
    return np.eye(3) + np.sin(t) * k + (1 - np.cos(t)) * k @ k


def _move(structure, chains, rotation, shift):
    atoms = [a for key in chains for a in structure[key].chain.get_atoms()]
    center = np.mean([a.coord for a in atoms], axis=0)
    for atom in atoms:
        atom.set_coord((atom.coord - center) @ rotation.T + center + shift)


def test_kabsch_recovers_rigid_motion():
    rng = np.random.default_rng(1)
    mobile = rng.normal(size=(30, 3))
    rotation = _rotation([1, 2, 3], 70)
    target = mobile @ rotation.T + [4, -2, 1]
    fit, shift = _kabsch(mobile, target)
    assert np.allclose(mobile @ fit + shift, target, atol=1e-8)
    assert np.isclose(np.linalg.det(fit), 1)
    with pytest.raises(ValueError, match="at least 3"):
        _kabsch(mobile[:2], target[:2])


def test_rmsd_and_superposition(ref):
    model = ref.subset_chains(ref.chain_ids, structure_id="model")
    _move(model, model.chain_ids, _rotation([0, 1, 1], 50), [8, 0, -3])
    assert model.rmsd_metrics(ref)["rmsd_A"] < 1e-3
    assert model.rmsd_metrics(ref, superpose=False)["rmsd_A"] > 5
    assert model.apply_superposition(ref, chains="A") is model
    assert model.rmsd_metrics(ref, superpose=False)["rmsd_A"] < 1e-3
    for atoms in ("ca", "backbone", "heavy"):
        metrics = model.rmsd_metrics(ref, chains=["A", "D"], atoms=atoms)
        assert metrics["n_residues"] == len(model["A"].residues) + len(model["D"].residues)
    with pytest.raises(ValueError, match="atoms must be"):
        model.rmsd_metrics(ref, atoms="side")


def test_fit_chains_gives_ligand_rmsd(ref):
    model = ref.subset_chains(["A", "D"], structure_id="model")
    _move(model, ["D"], np.eye(3), [3, 0, 0])
    assert model.rmsd_metrics(ref, chains="D", fit_chains="A")["rmsd_A"] == pytest.approx(3, abs=1e-3)
    assert model.rmsd_metrics(ref, chains="D")["rmsd_A"] < 1e-3


def test_residue_matching_across_entries(ref):
    barnase = SyncedStructure.from_path(DATA / "1A2P.cif")
    table = barnase.residue_match_table(ref, {"A": "A"})
    assert len(table) > 100 and (table["resname"] == table["reference_resname"]).mean() > 0.95
    assert table.columns[0] == "structure_id" and table["sequence_index"].min() >= 1
    with pytest.raises(ValueError, match="chain_map"):
        barnase.residue_match_table(ref, {})


def test_capri_metrics(ref):
    model = ref.subset_chains(["A", "D"], structure_id="model")
    same = model.capri_metrics(ref, "A", "D")
    assert (same["irmsd_A"], same["lrmsd_A"], same["fnat"], same["dockq"]) == (0, 0, 1, 1)
    assert same["capri_class"] == "high"
    scores = []
    for degrees, shift in ((10, 2), (30, 6)):
        moved = ref.subset_chains(["A", "D"], structure_id=f"m{degrees}")
        _move(moved, ["D"], _rotation([1, 0, 1], degrees), [shift, 0, 0])
        scores.append(moved.capri_metrics(ref, ["A"], ["D"]))
    assert scores[0]["dockq"] > scores[1]["dockq"] and scores[0]["lrmsd_A"] < scores[1]["lrmsd_A"]
    table = ref.capri_table({"same": model, "far": moved}, "A", "D")
    assert list(table["structure_id"]) == ["same", "far"]
    with pytest.raises(ValueError, match="share"):
        model.capri_metrics(ref, "A", "A")
    with pytest.raises(ValueError, match="no contact"):
        ref.capri_metrics(ref, "A", "E")


def test_capri_class_and_dockq():
    assert capri_class(0.6, 0.8, 3) == "high"
    assert capri_class(0.4, 3, 4) == "medium"
    assert capri_class(0.2, 5, 9) == "acceptable"
    assert capri_class(0.05, 0.5, 0.5) == "incorrect"
    assert dockq_score(1, 0, 0) == 1


@pytest.mark.haddock
def test_capri_matches_haddock_caprieval(ref, tmp_path):
    libcapri = pytest.importorskip("haddock.libs.libcapri")
    native = ref.subset_chains(["A", "D"], structure_id="native")
    native_pdb = native.export_structure("pdb", outdir=tmp_path)
    model = native.subset_chains(["A", "D"], structure_id="model")
    _move(model, ["D"], _rotation([1, 2, 0], 15), [0, 3, 0])
    params = {"allatoms": False, "receptor_chain": "A", "ligand_chains": ["D"], "keep_hetatm": False}
    capri = libcapri.CAPRI(1, model.export_structure("pdb", outdir=tmp_path), tmp_path, native_pdb, params)
    capri.calc_fnat(5.0)
    capri.calc_irmsd(10.0)
    capri.calc_lrmsd()
    ours = model.capri_metrics(native, "A", "D")
    assert ours["lrmsd_A"] == pytest.approx(capri.lrmsd, abs=0.01)
    assert ours["irmsd_A"] == pytest.approx(capri.irmsd, abs=0.05)
    assert ours["fnat"] == pytest.approx(capri.fnat, abs=0.02)


def test_chain_graph(ref):
    table = ref.chain_graph_table()
    assert table.iloc[0]["contact_pairs"] >= table.iloc[-1]["contact_pairs"]
    slow = ref.chain_graph_table(contacts=ref.contact_table(contact_type="inter"))
    key = ["chain_key_a", "chain_key_b"]
    assert table.sort_values(key).reset_index(drop=True).equals(slow.sort_values(key).reset_index(drop=True))
    fig, ax = ref.plot_chain_graph(min_contacts=2)
    widths = sorted(line.get_linewidth() for line in ax.lines)
    assert len(ax.lines) == int((table["contact_pairs"] >= 2).sum()) and widths[0] < widths[-1]
    ref.plot_chain_graph(layout="circular", labels="auth", describe=False, show_counts=False, ax=ax)
    with pytest.raises(ValueError, match="layout"):
        ref.plot_chain_graph(layout="tree")


def test_superposition_view(ref):
    model = ref.subset_chains(ref.chain_ids, structure_id="model")
    _move(model, model.chain_ids, _rotation([1, 1, 0], 30), [5, 5, 5])
    before = model["A"].residues[0]["CA"].coord.copy()
    html = model.show_superposition_3d(ref, chains=["A", "D"]).data
    assert np.array_equal(model["A"].residues[0]["CA"].coord, before)
    assert "RMSD 0.00" in html and html.count("extra=[{") == 1
    assert "extra=[]" in ref.show_structure_3d().data
