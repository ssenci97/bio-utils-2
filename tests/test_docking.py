"""HADDOCK3 set up on 1AAY: protein chain A docked to the DNA duplex B+C.

Offline tests need no haddock3. Tests marked 'haddock' call the haddock3 executable with low-throughput settings
(run with: pytest -m haddock).
"""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.PDB import PDBParser

from synced_structure import HaddockDocking, SyncedStructure, core_count
from synced_structure.constants import HADDOCK_PRESETS, HADDOCK_SEEDED_MODULES
from synced_structure.docking import cluster_contact_patches

DATA = Path(__file__).parent / "data"
CIF = DATA / "1AAY.cif"
LOW_THROUGHPUT = {"rigidbody": {"sampling": 8}, "clustfcc": {"min_population": 1}}


@pytest.fixture
def sync():
    """Return a fresh SyncedStructure of 1AAY keyed by author chain IDs (A protein, B and C DNA)."""
    return SyncedStructure.from_path(CIF)


@pytest.fixture
def setup(sync, tmp_path):
    """Return an oracle set up docking chain A against the B+C duplex."""
    return sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "dock", method="oracle", seed=7)


def _config(setup):
    tomllib = pytest.importorskip("tomllib")
    return tomllib.loads(setup.config_toml())


def _heavy_coords(structure, keys):
    return np.array([atom.coord for key in keys for residue in structure.chains[key].residues
                     for atom in residue.get_atoms() if atom.element != "H"], dtype=np.float64)


# partners, numbering, and placement

def test_segids_and_author_numbering(setup):
    frame = setup.residue_table()
    ranges = frame.groupby(["segid", "chain_key"])["haddock_resid"].agg(["min", "max", "count"])
    assert ranges.loc[("A", "A")].tolist() == [103, 187, 85]
    assert ranges.loc[("B", "B")].tolist() == [1, 11, 11]
    assert ranges.loc[("B", "C")].tolist() == [51, 61, 11]
    assert (frame["haddock_resid"] == frame["pdb_resnum"]).all()
    assert frame.loc[frame["chain_key"] == "A", "reference_index"].iloc[0] == 3


def test_partner_renumbered_when_author_numbers_collide(sync, tmp_path):
    sync.renumber_chain("B")
    sync.renumber_chain("C")
    frame = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path).residue_table()
    duplex = frame[frame["segid"] == "B"].groupby("chain_key")["haddock_resid"].agg(["min", "max"])
    assert duplex.loc["B"].tolist() == [1, 11]
    assert duplex.loc["C"].tolist() == [22, 32]
    assert frame.loc[frame["segid"] == "A", "haddock_resid"].min() == 103


def test_partners_validated(sync, tmp_path):
    with pytest.raises(ValueError):
        sync.prepare_haddock_docking("A", workdir=tmp_path)
    with pytest.raises(ValueError):
        sync.prepare_haddock_docking("A", ["A", "B"], workdir=tmp_path)
    with pytest.raises(KeyError):
        sync.prepare_haddock_docking("A", "Z", workdir=tmp_path)
    with pytest.raises(ValueError):
        sync.prepare_haddock_docking("A", "B", workdir=tmp_path, method="docking")


def test_same_structure_partners_are_separated(sync, setup, tmp_path):
    before = _heavy_coords(sync, ["A", "B", "C"]).copy()
    assert setup.separated and list(setup.placements) == [1]
    moved = setup._place_coords(1, _heavy_coords(sync, ["B", "C"]))
    distances = np.linalg.norm(_heavy_coords(sync, ["A"])[:, None] - moved[None], axis=-1)
    assert distances.min() > setup.clearance_a
    rotation = setup.placements[1]["rotation"]
    assert np.allclose(rotation @ rotation.T, np.eye(3)) and np.isclose(np.linalg.det(rotation), 1.0)
    assert np.array_equal(before, _heavy_coords(sync, ["A", "B", "C"]))
    same = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "same", seed=7)
    other = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "other", seed=8)
    assert np.allclose(same.placements[1]["center"], setup.placements[1]["center"])
    assert not np.allclose(other.placements[1]["center"], setup.placements[1]["center"])


def test_separation_options(sync, tmp_path):
    assert not sync.prepare_haddock_docking("A", "B", workdir=tmp_path, separate=False).placements
    touching = sync.subset_chains(["B", "C"])         # another structure ID, same coordinates: still moved apart
    assert touching.structure_id != sync.structure_id
    assert sync.prepare_haddock_docking("A", (touching, ["B", "C"]), workdir=tmp_path).separated
    apart = sync.subset_chains(["B", "C"])
    for atom in apart.structure.get_atoms():
        atom.coord = atom.coord + np.array([500.0, 0.0, 0.0], dtype=atom.coord.dtype)
    assert not sync.prepare_haddock_docking("A", (apart, ["B", "C"]), workdir=tmp_path).separated
    assert sync.prepare_haddock_docking("A", (apart, ["B", "C"]), workdir=tmp_path, separate=True).separated


# restraints

def test_positions_resolve_in_every_numbering(setup):
    assert setup._haddock_numbers("A", [3, 4], "reference") == (0, [103, 104])
    assert setup._haddock_numbers("A", [1, 2], "sequence") == (0, [103, 104])
    assert setup._haddock_numbers("A", [103, 104], "author") == (0, [103, 104])
    assert setup._haddock_numbers("C", [1], "reference") == (1, [51])
    assert setup._haddock_numbers("A", [1, 3], "reference") == (0, [103])
    with pytest.raises(ValueError):
        setup._haddock_numbers("A", [1, 2], "reference")
    with pytest.raises(ValueError):
        setup._haddock_numbers("A", [3], "pdb")


def test_ambiguous_restraints(setup):
    setup.add_ambiguous_restraints("A", [3, 4], "B", [2], passive_b=[3])
    text = setup.restraint_text("ambiguous")
    assert text.count("assign") == 3
    assert "assign (resi 103 and segid A)" in text and "assign (resi 104 and segid A)" in text
    assert "assign (resi 2 and segid B)" in text and "assign (resi 3 and segid B)" not in text
    first = text.split("\n\n")[0]
    assert "(resi 2 and segid B)" in first and "(resi 3 and segid B)" in first
    assert text.splitlines()[-1].strip().endswith("2.0 2.0 0.0")
    _validate(text)


def test_unambiguous_restraints_and_wrapping(setup):
    setup.add_unambiguous_restraint("A", [3, 4, 5, 9], "C", [1], numbering="reference", distance=(3.0, 3.0, 1.0))
    text = setup.restraint_text("unambiguous")
    assert "(segid A and (resi 103:105 or resi 109))" in text
    assert "(segid B and resi 51) 3.0 3.0 1.0" in text
    setup.add_unambiguous_restraint("A", list(range(1, 90, 2)), "B", [2, 5, 8], numbering="sequence")
    text = setup.restraint_text("unambiguous")
    assert max(len(line) for line in text.splitlines()) < 132
    _validate(text)


def test_restraints_need_two_partners(setup, sync, tmp_path):
    with pytest.raises(ValueError):
        setup.add_unambiguous_restraint("B", [1], "C", [1], numbering="sequence")
    twin = SyncedStructure.from_path(CIF)
    pair = HaddockDocking([(sync, "A"), (twin, "A")], tmp_path / "twin")
    with pytest.raises(ValueError):
        pair.add_unambiguous_restraint("A", [3], "A", [3])
    pair.add_unambiguous_restraint((1, "A"), [3], (2, "A"), [3])
    assert pair.restraint_table()[["segid_a", "segid_b"]].iloc[0].tolist() == ["A", "B"]
    duplex = sync.subset_chains(["B", "C"])
    assert duplex["B"].chain_id == "A"
    keyed = sync.prepare_haddock_docking("A", (duplex, ["B", "C"]), workdir=tmp_path / "keyed")
    assert keyed._locate_chain("A") == (0, "A")


def test_patch_pairs(setup):
    setup.add_patch_pair_restraints([(("A", [20, 23]), ("B", [3, 4])), (("A", [50]), ("C", [2]))])
    table = setup.restraint_table()
    assert table["kind"].eq("ambiguous").all() and table["origin"].eq("patch_pair").all()
    assert table["residues_a"].tolist() == [[120, 123], [150]]
    assert table["residues_b"].tolist() == [[3, 4], [52]]
    setup.add_patch_pair_restraints([(("A", [20]), ("B", [3]))], kind="unambiguous")
    assert setup.restraint_table()["kind"].iloc[-1] == "unambiguous"
    with pytest.raises(ValueError):
        setup.add_patch_pair_restraints([])


def test_cluster_patches():
    points = [(1, 1), (2, 3), (10, 10), (11, 14), (30, 1)]
    assert cluster_contact_patches(points, gap=4) == [[0, 1], [2, 3], [4]]
    assert cluster_contact_patches(points, gap=0) == [[0], [1], [2], [3], [4]]


def test_contact_restraints(setup, sync):
    contacts = sync.contact_table(contact_type="inter", chain_ids=["A", "B", "C"])
    protein = contacts[(contacts["auth_chain_a"] == "A") | (contacts["auth_chain_b"] == "A")]
    expected = set()
    for row in protein.itertuples():
        side = "a" if row.auth_chain_a == "A" else "b"
        expected.add(getattr(row, f"resnum_{side}"))
    setup.add_contact_restraints(kind="ambiguous")
    table = setup.restraint_table()
    assert set().union(*table["residues_a"]) == expected
    setup.clear_restraints().add_contact_restraints(kind="unambiguous", scheme="sidechain", cutoff=0.45)
    patches = setup.contact_patch_table()
    assert len(patches) == len(setup.restraint_table()) and patches["n_pairs"].sum() > len(patches)
    assert set(patches["chain_b"]) == {"B", "C"}


def test_contacts_transfer_through_reference_positions(sync, tmp_path):
    direct = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "direct")
    direct.add_contact_restraints(kind="ambiguous")
    apo = SyncedStructure.from_path(CIF)
    apo.renumber_chain("A")
    moved = apo.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "moved")
    with pytest.raises(ValueError):
        HaddockDocking([(apo, "A"), (sync, ["B", "C"])], tmp_path / "mixed").add_contact_restraints()
    moved.add_contact_restraints(kind="ambiguous", source=sync, chain_map={"A": "A", "B": "B", "C": "C"})
    shift = lambda residues: [resid - 102 for resid in residues]
    assert moved.restraint_table()["residues_a"].tolist() == [shift(r) for r in direct.restraint_table()["residues_a"]]
    assert moved.restraint_table()["residues_b"].tolist() == direct.restraint_table()["residues_b"].tolist()


# partners from different structures: apo barnase (1A2P) and apo barstar (1A19), complex 1BRS (A barnase, D barstar)

BARNASE_BARSTAR = {(1, "A"): "A", (2, "A"): "D"}


@pytest.fixture
def apo():
    """Return apo barnase, apo barstar, and the barnase-barstar complex."""
    return tuple(SyncedStructure.from_path(DATA / f"{pdb_id}.cif") for pdb_id in ("1A2P", "1A19", "1BRS"))


def test_partners_from_different_structures(apo, tmp_path):
    barnase, barstar, holo = apo
    dock = HaddockDocking([(barnase, "A"), (barstar, "A")], tmp_path, method="oracle")
    frame = dock.residue_table()
    assert frame.groupby("segid")["structure_id"].unique().map(list).to_dict() == {"A": ["1A2P"], "B": ["1A19"]}
    assert not dock.separated and dock._reference_path() is None
    with pytest.raises(ValueError, match="several partners"):
        dock._locate_chain("A")
    with pytest.raises(ValueError, match="set_reference"):
        dock.workflow_modules()                         # oracle needs contacts, and there is no complex to read
    with pytest.raises(ValueError, match="chain_map"):
        dock.set_reference(holo)                        # both partners have a chain 'A'


def test_reference_structure_written_in_haddock_numbering(apo, tmp_path):
    barnase, barstar, holo = apo
    dock = HaddockDocking([(barnase, "A"), (barstar, "A")], tmp_path, method="oracle")
    dock.set_reference(holo, chain_map=BARNASE_BARSTAR).export_inputs()
    reference = dock.input_paths["reference"]
    assert dock._caprieval_params()["reference_fname"] == str(reference)
    parsed = PDBParser(QUIET=True).get_structure("ref", reference)[0]
    assert sorted(chain.id for chain in parsed) == ["A", "B"]
    rows = dock._reference_rows()
    assert sum(1 for _ in parsed.get_residues()) == len(rows) > 190
    written = {(chain.id, residue.id[1]) for chain in parsed for residue in chain}
    assert written <= set(zip(dock.residue_table()["segid"], dock.residue_table()["haddock_resid"]))
    # 1A19 barstar has an extra N-terminal Met and Cys40, 1BRS has Ala40: aligned, so the same residue lands on the same number
    barstar_rows = {entry["resid"]: residue for index, entry, residue in rows if index == 1}
    assert barstar_rows[40].get_resname() == "ALA" and barstar_rows[39].get_resname() == "ASP"
    assert barstar_rows[38].id[1] == 38 and barstar_rows[38].get_resname() == barstar.chains["A"].residues[37].get_resname()


def test_contacts_from_reference_structure(apo, tmp_path):
    barnase, barstar, holo = apo
    dock = HaddockDocking([(barnase, "A"), (barstar, "A")], tmp_path, method="oracle", reference=None)
    dock.set_reference(holo, chain_map=BARNASE_BARSTAR)
    modules = dict(dock.workflow_modules())             # oracle adds contact patches from the reference
    assert "unambig_fname" in modules["rigidbody"] and len(dock.restraint_table()) > 0
    explicit = HaddockDocking([(barnase, "A"), (barstar, "A")], tmp_path / "explicit")
    explicit.add_contact_restraints(kind="unambiguous", source=holo, chain_map=BARNASE_BARSTAR)
    assert explicit.restraint_table()["residues_b"].tolist() == dock.restraint_table()["residues_b"].tolist()
    barstar_contacts = sorted({resid for resids in dock.restraint_table()["residues_b"] for resid in resids})
    assert 35 in barstar_contacts                                  # Asp35 binds the barnase active site


def test_contacts_transfer_across_different_reference_sequences(sync, tmp_path):
    direct = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "direct")
    direct.add_contact_restraints(kind="ambiguous")
    tagged = SyncedStructure.from_path(CIF)
    tagged.set_reference_sequence("A", "MGSSHHHHHH" + sync.chains["A"].reference_sequence)
    moved = tagged.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "tagged")
    moved.add_contact_restraints(kind="ambiguous", source=sync)
    assert moved.restraint_table()["residues_a"].tolist() == direct.restraint_table()["residues_a"].tolist()


def test_reference_path_and_constructor_structure(sync, tmp_path):
    duplex = sync.subset_chains(["B", "C"])
    dock = HaddockDocking([(sync, "A"), (duplex, ["B", "C"])], tmp_path, reference=sync)
    assert dock._reference_path() == dock.input_paths["reference"]
    dock.set_reference(None)
    assert dock._reference_path() is None and "reference_fname" not in dock._caprieval_params()
    with pytest.raises(FileNotFoundError):
        dock.set_reference(tmp_path / "missing.pdb")._reference_path()


def _validate(text):
    librestraints = pytest.importorskip("haddock.libs.librestraints")
    librestraints.validate_tbldata(text)


# configuration

def test_preset_modules_seed_and_restraints(sync, tmp_path):
    for method in HADDOCK_PRESETS:
        setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / method, method=method, seed=42)
        if method == "targeted":
            setup.add_patch_pair_restraints([(("A", [20, 23]), ("B", [3, 4]))])
        config = _config(setup)
        modules = [name for name, _ in HADDOCK_PRESETS[method]] + ["caprieval"]
        assert [key for key in config if isinstance(config[key], dict)] == modules
        for name in modules:
            if name in HADDOCK_SEEDED_MODULES:
                assert config[name]["iniseed"] == 42
        assert config["caprieval"]["receptor_chain"] == "A" and config["caprieval"]["ligand_chains"] == ["B"]
        assert len(config["molecules"]) == (1 if method == "alascan" else 2) and config["run_dir"].endswith("run")
    abinitio = _config(sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "a"))
    assert abinitio["rigidbody"] == {"sampling": 10000, "cmrest": True, "ranair": True, "iniseed": 917}
    assert abinitio["emref"]["contactairs"] is True and abinitio["seletopclusts"] == {"top_models": 4}
    oracle = _config(sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "o", method="oracle"))
    assert oracle["rigidbody"]["sampling"] == 200 and "top_clusters" not in oracle["seletopclusts"]
    for config in (abinitio, oracle):
        assert [key for key in config if isinstance(config[key], dict)][-4:] == [
            "clustfcc", "seletopclusts", "emref", "caprieval"]
    assert oracle["rigidbody"]["unambig_fname"] == oracle["emref"]["unambig_fname"]
    assert "ambig_fname" not in oracle["rigidbody"]


def test_targeted_uses_patch_pairs_as_ambiguous(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="targeted")
    with pytest.raises(ValueError):
        setup.config_toml()
    setup.add_patch_pair_restraints([(("A", [20, 23]), ("B", [3, 4]))])
    config = _config(setup)
    assert config["rigidbody"]["ambig_fname"].endswith("ambig.tbl") and "unambig_fname" not in config["rigidbody"]
    no_reference = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="targeted", reference=None)
    no_reference.add_patch_pair_restraints([(("A", [20]), ("B", [3]))])
    capri = _config(no_reference)["caprieval"]
    assert "reference_fname" not in capri and capri["receptor_chain"] == "A"


def test_overrides_and_guards(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="ab-initio", ncores=2,
                                 module_params=LOW_THROUGHPUT, global_params={"debug": True}, caprieval=False)
    setup.add_ambiguous_restraints("A", [20], "B", [3])
    config = _config(setup)
    assert config["rigidbody"]["sampling"] == 8 and config["clustfcc"]["min_population"] == 1
    assert config["rigidbody"]["ranair"] is False and config["rigidbody"]["cmrest"] is True
    assert config["emref"]["ambig_fname"] == config["rigidbody"]["ambig_fname"]
    assert config["debug"] is True and config["ncores"] == 2 and "caprieval" not in config
    assert sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "default").ncores == core_count()
    with pytest.raises(ValueError):
        sync.prepare_haddock_docking("A", "B", workdir=tmp_path, module_params={"rigidbody": {"iniseed": 1}})
    with pytest.raises(ValueError):
        sync.prepare_haddock_docking("A", "B", workdir=tmp_path, module_params={"flexref": {}})
    with pytest.raises(ValueError):
        sync.prepare_haddock_docking("A", "B", workdir=tmp_path, global_params={"run_dir": "x"})
    apart = sync.subset_chains(["B", "C"])
    unreferenced = _config(sync.prepare_haddock_docking("A", (apart, ["B", "C"]), workdir=tmp_path))
    assert "reference_fname" not in unreferenced["caprieval"]
    assert unreferenced["caprieval"]["receptor_chain"] == "A"


# files

def test_written_files(sync, setup):
    config = setup.export_inputs()
    paths = setup.input_paths
    assert config.is_file() and paths["unambiguous"].is_file() and not paths["ambiguous"].exists()
    written = sorted(path.name for path in setup.workdir.iterdir() if path.is_file())
    assert written == ["config.toml", "haddock_residue_map.tsv"]
    assert len(paths["haddock_residue_map"].read_text().splitlines()) == 85 + 11 + 11 + 1
    parser = PDBParser(QUIET=True)
    receptor = parser.get_structure("a", paths["molecules"][0])
    duplex = parser.get_structure("b", paths["molecules"][1])
    assert [chain.id for chain in receptor[0]] == ["A"] and [chain.id for chain in duplex[0]] == ["B"]
    assert {atom.get_parent().segid.strip() for atom in duplex.get_atoms()} == {"B"}
    assert len(list(duplex.get_residues())) == 22
    written = np.array([atom.coord for atom in receptor.get_atoms()])
    assert np.allclose(written, _heavy_coords(sync, ["A"]), atol=1e-3)
    reference = parser.get_structure("r", paths["reference"])
    assert np.allclose([atom.coord for atom in reference[0]["B"].get_atoms()], _heavy_coords(sync, ["B", "C"]),
                       atol=1e-3)
    assert all(atom.element != "H" for atom in reference.get_atoms())


def test_load_model_restores_chains(sync, setup):
    setup.export_inputs()
    model = setup.load_docked_model(setup.input_paths["reference"])
    assert sorted(model.chain_ids) == ["A", "B", "C"]
    for key in ("A", "B", "C"):
        assert model.get_sequence(key) == sync.get_sequence(key)
        assert model.get_reference_sequence(key) == sync.get_reference_sequence(key)
        assert [r.id for r in model[key].residues] == [r.id for r in sync[key].residues]
    assert model.chain_residue_table("A")["reference_index"].iloc[0] == 3


def test_run_guards(setup):
    setup.export_inputs()
    setup.run_dir.mkdir()
    (setup.run_dir / "keep.txt").write_text("not ours")
    with pytest.raises(FileExistsError):
        setup.run_haddock()
    with pytest.raises(FileExistsError):
        setup.run_haddock(overwrite=True)
    assert (setup.run_dir / "keep.txt").is_file()
    partners = [(partner["structure"], partner["keys"]) for partner in setup.partners]
    missing = HaddockDocking(partners, setup.workdir / "x", haddock3_exe="/nonexistent/haddock3")
    with pytest.raises(FileNotFoundError):
        missing.validate_config()



def test_model_counts(setup):
    modules = ["topoaa", "rigidbody", "clustfcc", "seletopclusts", "emref", "caprieval"]
    assert [name for name, _ in setup.workflow_modules()] == modules
    assert setup._model_counts() == [2, 200, 0, 0, None, 0]       # emref is unknown until clusters are counted
    setup._top_clusters = 5
    assert setup._model_counts() == [2, 200, 0, 0, 20, 0]
    assert setup._cluster_split() == 3


def test_top_cluster_count(setup):
    folder = setup.run_dir / "2_clustfcc"
    folder.mkdir(parents=True)
    with pytest.raises(RuntimeError, match="no clusters"):             # clustfcc writes no table without clusters
        setup._top_cluster_count(2)
    for found, top in ((1, 1), (3, 3), (5, 5), (12, 5), (25, 5), (26, 6), (100, 20)):
        rows = [f"{rank}\trigidbody_{rank}.pdb\t-1.00\t{cluster}"
                for rank, cluster in enumerate([*range(1, found + 1), 1, "-"], start=1)]
        (folder / "clustfcc.tsv").write_text("rank\tmodel_name\tscore\tcluster_id\n" + "\n".join(rows) + "\n\n")
        assert setup._top_cluster_count(2) == top
    (folder / "clustfcc.tsv").write_text("rank\tmodel_name\tscore\tcluster_id\n1\trigidbody_1.pdb\t-1.00\t-\n")
    with pytest.raises(RuntimeError, match="no clusters"):
        setup._top_cluster_count(2)
    fixed = setup.__class__([(p["structure"], p["keys"]) for p in setup.partners], setup.workdir / "f",
                            method="oracle", module_params={"seletopclusts": {"top_clusters": 3}})
    assert fixed._cluster_split() is None and fixed._model_counts()[4] == 12


FAKE_HADDOCK3 = """#!/bin/sh
# stands in for haddock3: logs like it and writes model files step by step
cd "$(dirname "$1")"
for step in 0_topoaa 1_rigidbody 2_clustfcc 3_seletopclusts 4_emref 5_caprieval; do
    name=${step#*_}
    echo "[INFO] Running [$name] module"
    mkdir -p run/$step
    case $name in
        topoaa) n=2 ;; rigidbody) n=8 ;; emref) n=4 ;; *) n=0 ;;
    esac
    i=1
    while [ $i -le $n ]; do touch run/$step/${name}_$i.pdb; i=$((i + 1)); sleep 0.05; done
    if [ $name = clustfcc ]; then
        printf 'rank\tmodel_name\tscore\tcluster_id\n1\trigidbody_1.pdb\t-9.00\t1\n2\trigidbody_2.pdb\t-8.00\t2\n3\trigidbody_3.pdb\t-7.00\t-\n' > run/$step/clustfcc.tsv
    fi
    echo "[INFO] Module [$name] finished."
done
exit ${FAKE_EXIT:-0}
"""


def test_run_progress(sync, tmp_path, caplog, monkeypatch):
    import synced_structure.docking as docking

    monkeypatch.setattr(docking, "HADDOCK_POLL_SECONDS", 0.02)
    monkeypatch.setattr(docking, "HADDOCK_PROGRESS_SECONDS", 0.0)
    exe = tmp_path / "haddock3"
    exe.write_text(FAKE_HADDOCK3)
    exe.chmod(0o755)
    params = {**LOW_THROUGHPUT, "seletopclusts": {"top_clusters": 1}}       # one haddock3 call, no restart
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "dock", method="oracle",
                                         module_params=params, haddock3_exe=exe)
    with caplog.at_level("INFO"):
        setup.run_haddock()
    lines = [record.getMessage() for record in caplog.records if record.getMessage().startswith("HADDOCK3 [")]
    modules = ["topoaa", "rigidbody", "clustfcc", "seletopclusts", "emref", "caprieval"]
    for index, name in enumerate(modules, start=1):
        assert any(line.startswith(f"HADDOCK3 [{index}/6] {name} started") for line in lines)
        assert any(line.startswith(f"HADDOCK3 [{index}/6] {name} finished in") for line in lines)
    assert any(line.startswith("HADDOCK3 [2/6] rigidbody: 8/8 models") for line in lines)
    assert any(line.startswith("HADDOCK3 [5/6] emref: 4/4 models") for line in lines)
    assert not any("clustfcc:" in line or "seletopclusts:" in line for line in lines)
    monkeypatch.setenv("FAKE_EXIT", "3")
    with pytest.raises(RuntimeError, match="exited with code 3"):
        setup.run_haddock(overwrite=True, progress=False)


def test_run_cluster_split_without_haddock_import(sync, tmp_path, monkeypatch):
    import importlib.util

    find_spec = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda name, *args: None if name.split(".")[0] in
                        ("haddock", "jsonpickle") else find_spec(name, *args))
    exe = tmp_path / "haddock3"
    exe.write_text(FAKE_HADDOCK3)
    exe.chmod(0o755)
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "dock", method="oracle",
                                         module_params=LOW_THROUGHPUT, haddock3_exe=exe)
    assert setup._cluster_split() == 3
    setup.run_haddock(progress=False)
    assert setup._top_clusters == 2                                     # clustfcc.tsv has clusters 1, 2 and '-'
    assert "top_clusters = 2" in (setup.workdir / "config.toml").read_text()


def test_interrupt_stops_haddock3(sync, tmp_path, monkeypatch):
    import os
    import time

    exe = tmp_path / "haddock3"
    exe.write_text('#!/bin/sh\necho $$ > "$(dirname "$1")/pid"\nexec sleep 60\n')
    exe.chmod(0o755)
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "dock", method="oracle",
                                         module_params={"seletopclusts": {"top_clusters": 1}}, haddock3_exe=exe)

    def interrupt(process, log_path, offset=0):
        while not (setup.workdir / "pid").is_file():
            time.sleep(0.02)
        raise KeyboardInterrupt

    monkeypatch.setattr(setup, "_watch_run", interrupt)
    with pytest.raises(KeyboardInterrupt):
        setup.run_haddock()
    pid = int((setup.workdir / "pid").read_text())
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)

def test_duration():
    from synced_structure.docking import _duration

    assert [_duration(s) for s in (4.4, 185, 3725)] == ["4s", "3m05s", "1h02m"]


# results of a finished run, read from a run directory laid out like haddock3's

REMARKS = """REMARK FILENAME="{name}"
REMARK ===============================================================
REMARK HADDOCK score: {score}
REMARK            total,bonds,angles,improper,dihe,vdw,elec,air,cdih,coup,rdcs,vean,dani,xpcs,rg
REMARK energies: -363.8, 28.9, 106.6, 25.7, 463.3, {vdw}, -336.6, 0, 0, 0, 0, 0, 0, 0, 0
REMARK            bonds,angles,impropers,dihe,air,cdih,coup,rdcs,vean,dani,xpcs
REMARK rms-dev.: 4.1E-03,0.50,0.39,42.3,42.3,0,0, 0, 0, 0, 0
REMARK Desolvation energy: 10.8
REMARK buried surface area: {bsa}
"""
CAPRI_SS = ("model\tmd5\tcaprieval_rank\tscore\tirmsd\tfnat\tlrmsd\tilrmsd\tdockq\tcluster_id\n"
            "../1_rigidbody/rigidbody_2.pdb\t-\t1\t-90.0\t0.0\t1.0\t0.0\t0.0\t1.0\t-\n"
            "../1_rigidbody/rigidbody_1.pdb\t-\t2\t-50.0\t3.1\t0.4\t7.2\t5.0\t0.3\t-\n")


@pytest.fixture
def finished(setup):
    """Return the set up with a fake finished run: two rigidbody models (one compressed) and a caprieval step."""
    import gzip

    setup.export_inputs()
    rigid, capri = setup.run_dir / "1_rigidbody", setup.run_dir / "2_caprieval"
    rigid.mkdir(parents=True)
    capri.mkdir()
    atoms = setup.input_paths["reference"].read_text()
    (rigid / "rigidbody_1.pdb").write_text(REMARKS.format(name="rigidbody_1.pdb", score=-50.0, vdw=-10.0, bsa=900.0)
                                           + atoms)
    with gzip.open(rigid / "rigidbody_2.pdb.gz", "wt") as handle:
        handle.write(REMARKS.format(name="rigidbody_2.pdb", score=-90.0, vdw=-30.0, bsa=1300.0) + atoms)
    (capri / "capri_ss.tsv").write_text(CAPRI_SS)
    (capri / "capri_clt.tsv").write_text("# header\n#\ncluster_rank\tcluster_id\tn\tscore\n1\t1\t2\t-70.0\n")
    (capri / "params.cfg").write_text('reference_fname = ""\nreceptor_chain = "A"\n')
    return setup


def test_step_lookup(finished):
    assert list(finished.step_paths()) == ["1_rigidbody", "2_caprieval"]
    paths = {finished._step_path(step) for step in ("rigidbody", 1, "1_rigidbody")}
    assert paths == {finished.run_dir / "1_rigidbody"}
    with pytest.raises(KeyError):
        finished._step_path("emref")


def test_scores_from_model_headers(finished, monkeypatch):
    monkeypatch.setattr(HaddockDocking, "_read_io_rows", lambda self, folder: None)
    scores = finished.score_table("rigidbody")
    assert scores["model"].tolist() == ["rigidbody_2.pdb", "rigidbody_1.pdb"]
    assert scores["rank"].tolist() == [1, 2] and scores["score"].tolist() == [-90.0, -50.0]
    assert scores["vdw"].tolist() == [-30.0, -10.0] and scores["bsa"].tolist() == [1300.0, 900.0]
    assert scores["desolv"].eq(10.8).all() and "rms-dev." not in scores
    assert all(path.is_file() for path in scores["path"]) and scores["path"].iloc[0].suffix == ".gz"


def test_capri_tables_and_docked(finished, monkeypatch):
    monkeypatch.setattr(HaddockDocking, "_read_io_rows", lambda self, folder: None)
    tables = finished.caprieval_tables()
    models = tables["models"]
    assert models["caprieval_rank"].tolist() == [1, 2] and models["reference"].eq("best model").all()
    assert [path.name for path in models["path"]] == ["rigidbody_2.pdb.gz", "rigidbody_1.pdb"]
    assert tables["clusters"]["cluster_id"].tolist() == [1]
    assert finished.caprieval_tables()["models"].equals(models)
    docked = finished.docked_model_table()
    assert docked["path"].tolist() == finished.docked_model_paths() == models["path"].tolist()
    assert docked["dockq"].tolist() == [1.0, 0.3]
    with pytest.raises(FileNotFoundError):
        finished.caprieval_tables("rigidbody")
    loaded = finished.load_docked_models(top=1)
    assert list(loaded) == ["rigidbody_2"] and sorted(loaded["rigidbody_2"].chain_ids) == ["A", "B", "C"]


def test_plot_caprieval(finished):
    import matplotlib

    matplotlib.use("Agg")
    fig, ax = finished.plot_caprieval()
    assert len(ax.collections[0].get_offsets()) == 2
    assert ax.get_xlabel() == "i-RMSD (Å)" and ax.get_ylabel() == "HADDOCK score"
    assert fig.axes[-1].get_ylabel() == "DockQ" and ax.get_title() == "Models vs best model"
    fig, ax = finished.plot_caprieval(x="fnat", y="dockq", color=None, title="t")
    assert len(fig.axes) == 1 and ax.get_title() == "t"
    with pytest.raises(ValueError, match="choose from"):
        finished.plot_caprieval(x="cluster_id")
    with pytest.raises(ValueError):
        finished.plot_caprieval(color="nope")
    matplotlib.pyplot.close("all")


def test_docked_without_caprieval(finished, monkeypatch):
    monkeypatch.setattr(HaddockDocking, "_read_io_rows", lambda self, folder: None)
    shutil.rmtree(finished.run_dir / "2_caprieval")
    (finished.run_dir / "1_rigidbody" / "io.json").write_text("{}")
    docked = finished.docked_model_table()
    assert docked["model"].tolist() == ["rigidbody_2.pdb", "rigidbody_1.pdb"]
    assert finished.docked_model_paths() == docked["path"].tolist()


# alascan (no docking)

SCAN_HEADER = "#" * 80 + "\n# `alascan` results\n#\n"
SCAN_MODEL = (SCAN_HEADER + "chain\tres\tori_resname\tend_resname\tscore\tdelta_score\tz_score\n"
              "A\t146\tARG\tALA\t-250.0\t-16.0\t-1.0\nA\t180\tARG\tALA\t-260.0\t-8.0\t1.0\n")
SCAN_CLUSTER = (SCAN_HEADER + "chain\tresid\tresname\tfull_resname\tdelta_score\tdelta_score_std\tfrac_pres\tz_score\n"
                "A\t180\tARG\tA-180-ARG\t-8.0\t0.0\t0.5\t1.0\nA\t146\tARG\tA-146-ARG\t-16.0\t0.2\t1.0\t-1.0\n")


@pytest.fixture
def alascan(sync, tmp_path):
    """Return an alascan set up of chain A against the B+C duplex."""
    return sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "scan", method="alascan", seed=7)


def test_alascan_workflow_keeps_the_complex(sync, alascan, tmp_path):
    config = _config(alascan)
    assert [name for name in config if isinstance(config[name], dict)] == [
        "topoaa", "mdref", "rmsdmatrix", "clustrmsd", "alascan", "caprieval"]
    assert "rigidbody" not in config and len(config["molecules"]) == 1
    assert config["mdref"]["iniseed"] == 7 and config["caprieval"]["ligand_chains"] == ["B"]
    assert not alascan.separated and not alascan.placements
    assert alascan._model_counts() == [1, 10, 0, 0, 0, 0]
    with pytest.raises(ValueError, match="cannot be separated"):
        sync.prepare_haddock_docking("A", "B", workdir=tmp_path, method="alascan", separate=True)
    alascan.export_inputs()
    complex_pdb = alascan.input_paths["molecules"][0]
    assert complex_pdb.name == "complex_1AAY_AB.pdb"
    parsed = PDBParser(QUIET=True).get_structure("c", complex_pdb)
    assert sorted(chain.id for chain in parsed[0]) == ["A", "B"]
    written = np.array([atom.coord for atom in parsed.get_atoms()], dtype=np.float64)
    assert np.allclose(written, _heavy_coords(sync, ["A", "B", "C"]), atol=1e-3)


def test_alascan_warns_without_interface(sync, tmp_path):
    apart = sync.subset_chains(["B", "C"])
    for atom in apart.structure.get_atoms():
        atom.coord = atom.coord + np.array([500.0, 0.0, 0.0], dtype=atom.coord.dtype)
    with pytest.warns(UserWarning, match="no interface"):
        sync.prepare_haddock_docking("A", (apart, ["B", "C"]), workdir=tmp_path, method="alascan")


def test_scan_residues(alascan, setup):
    alascan.add_scan_residues("A", [24, 46]).add_scan_residues("A", [146], numbering="author")
    assert _config(alascan)["alascan"]["resdic_A"] == [124, 146]
    with pytest.raises(ValueError, match="alascan"):
        setup.add_scan_residues("A", [24])


def test_alascan_tables(alascan):
    with pytest.raises(FileNotFoundError):
        alascan.alascan_tables()
    folder = alascan.run_dir / "4_alascan"
    folder.mkdir(parents=True)
    (folder / "scan_mdref_1.tsv").write_text(SCAN_MODEL)
    (folder / "scan_clt_1.tsv").write_text(SCAN_CLUSTER)
    tables = alascan.alascan_tables()
    assert alascan.alascan_tables(4)["models"].equals(tables["models"])
    models, clusters = tables["models"], tables["clusters"]
    assert models["model"].eq("mdref_1").all() and models["haddock_resid"].tolist() == [146, 180]
    assert models["reference_index"].tolist() == [46, 80] and models["chain_key"].eq("A").all()
    assert clusters["cluster_id"].eq(1).all() and clusters["delta_score"].tolist() == [-16.0, -8.0]
    assert clusters["reference_index"].tolist() == [46, 80] and "full_resname" not in clusters
    assert list(models.columns[:3]) == ["model", "segid", "haddock_resid"]


# haddock3 runs (low throughput)

needs_haddock = pytest.mark.skipif(shutil.which("haddock3") is None, reason="haddock3 is not installed")


@pytest.mark.haddock
@needs_haddock
def test_haddock_validates_every_method(sync, tmp_path):
    for method in HADDOCK_PRESETS:
        setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / method, method=method, ncores=2)
        if method == "targeted":
            setup.add_patch_pair_restraints([(("A", [20, 23, 26]), ("B", [3, 4, 5]))])
        setup.validate_config()
        assert not (setup.workdir / "setup_check").exists()
    bad = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "bad", module_params={"rigidbody": {"samplin": 4}})
    with pytest.raises(RuntimeError):
        bad.validate_config()


@pytest.mark.haddock
@needs_haddock
def test_haddock_oracle_run(sync, tmp_path, caplog):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="oracle", seed=7, ncores=2,
                                 module_params=LOW_THROUGHPUT)
    setup.add_contact_restraints(kind="unambiguous", scheme="sidechain", cutoff=0.45)
    with caplog.at_level("INFO"):
        setup.run_haddock()
    messages = [record.getMessage() for record in caplog.records]
    assert any(message.startswith("HADDOCK3 [2/6] rigidbody: 8/8 models") for message in messages)
    assert sum(" started | run " in message for message in messages) == 6
    assert any(message.startswith("HADDOCK3 [6/6] caprieval finished in") for message in messages)
    clusters = setup.score_table("clustfcc")["cluster_id"].nunique()
    assert setup._top_clusters == min(clusters, 5)
    selected, emref = setup.score_table("seletopclusts"), setup.score_table("emref")
    assert len(emref) == len(selected) <= 4 * setup._top_clusters
    assert any(message.startswith(f"HADDOCK3 [5/6] emref: {len(emref)}/{len(emref)} models") for message in messages)
    models = setup.caprieval_tables()["models"]
    assert len(models) == len(emref) and {"irmsd", "fnat", "dockq"} <= set(models)
    assert models["reference"].eq("reference.pdb").all()
    assert all(path.is_file() for path in setup.docked_model_paths())
    assert len(setup.score_table("rigidbody")) == 8 and emref["origin"].str.startswith("cluster_").all()
    assert setup.score_table()["model"].tolist() == [path.name.removesuffix(".gz") for path in setup.docked_model_paths()]
    docked_table = setup.docked_model_table()
    assert docked_table["path"].tolist() == setup.docked_model_paths() and {"seed", "origin", "dockq"} <= set(docked_table)
    assert len(setup.load_docked_models(top=2)) == 2
    analysed = setup.run_caprieval(["rigidbody"])["1_rigidbody"]["models"]
    assert len(analysed) == 8 and analysed["path"].map(Path.is_file).all()
    assert setup.caprieval_tables("rigidbody")["models"]["reference"].eq("reference.pdb").all()
    docked = setup.load_docked_model(0)
    assert sorted(docked.chain_ids) == ["A", "B", "C"]
    assert docked.chain_table()["mapping_fraction"].eq(1.0).all()
    setup.run_haddock(overwrite=True)


@pytest.mark.haddock
@needs_haddock
def test_haddock_abinitio_run(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="ab-initio", seed=3, ncores=2,
                                 module_params=LOW_THROUGHPUT)
    setup.run_haddock()
    tables = setup.caprieval_tables()
    assert {"models", "clusters"} <= set(tables) and tables["models"]["cluster_id"].notna().all()
    assert setup.score_table("seletopclusts")["cluster_id"].notna().all()
    assert sorted(setup.load_docked_model(0).chain_ids) == ["A", "B", "C"]


@pytest.mark.haddock
@needs_haddock
def test_haddock_run_without_reference_still_scores(sync, tmp_path):
    duplex = sync.subset_chains(["B", "C"])
    setup = sync.prepare_haddock_docking("A", (duplex, ["B", "C"]), workdir=tmp_path, method="oracle", seed=5, ncores=2,
                                 module_params=LOW_THROUGHPUT)
    setup.add_unambiguous_restraint("A", [20, 23], "B", [3, 4])
    setup.run_haddock()
    docked = setup.docked_model_table()
    assert docked["reference"].eq("best model").all() and docked["irmsd"].iloc[0] == 0


@pytest.mark.haddock
@needs_haddock
def test_haddock_run_across_structures(apo, tmp_path):
    barnase, barstar, holo = apo
    setup = HaddockDocking([(barnase, "A"), (barstar, "A")], tmp_path, method="oracle", seed=11, ncores=2,
                           module_params=LOW_THROUGHPUT)
    setup.set_reference(holo, chain_map=BARNASE_BARSTAR)
    setup.run_haddock()
    models = setup.caprieval_tables()["models"]
    assert models["reference"].eq("reference.pdb").all() and models["irmsd"].notna().all()
    fig, _ = setup.plot_caprieval()
    assert fig.axes[-1].get_ylabel() == "DockQ"
    docked = setup.load_docked_model(0)
    assert sorted(docked.chain_ids) == ["A", "AB"]                  # both partners had a chain 'A'
    assert docked.chain_table()["mapping_fraction"].eq(1.0).all()


@pytest.mark.haddock
@needs_haddock
def test_haddock_alascan_run(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="alascan", seed=7, ncores=2,
                                         module_params={"mdref": {"sampling_factor": 2}, "clustrmsd": {"n_clusters": 1}})
    setup.add_scan_residues("A", [24, 46, 70, 80])           # arginines bound to the DNA
    setup.run_haddock()
    assert [name.split("_", 1)[1] for name in setup.step_paths()] == [
        "topoaa", "mdref", "rmsdmatrix", "clustrmsd", "alascan", "caprieval"]
    tables = setup.alascan_tables()
    models = tables["models"]
    assert len(models) and models["chain_key"].eq("A").all() and set(models["reference_index"]) <= {24, 46, 70, 80}
    assert set(models["model"]) == {"mdref_1", "mdref_2"} and len(tables["clusters"])
    docked = setup.docked_model_table()
    assert len(docked) == 2 and docked["reference"].eq("reference.pdb").all() and docked["irmsd"].lt(3).all()
    assert sorted(setup.load_docked_model(0).chain_ids) == ["A", "B", "C"]
