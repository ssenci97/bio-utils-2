"""HADDOCK3 set up on 1AAY: protein chain A docked to the DNA duplex B+C.

Offline tests need no haddock3. Tests marked 'haddock' call the haddock3 executable with low-throughput settings
(run with: pytest -m haddock).
"""

import shutil
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBParser

from synced_structure import HaddockDocking, SyncedStructure
from synced_structure.constants import HADDOCK_PRESETS, HADDOCK_SEEDED_MODULES
from synced_structure.docking import cluster_contact_patches

CIF = Path(__file__).parent / "data" / "1AAY.cif"
LOW_THROUGHPUT = {"rigidbody": {"sampling": 4}, "seletop": {"select": 2}}


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
    frame = setup.residue_table
    ranges = frame.groupby(["segid", "chain_key"])["haddock_resid"].agg(["min", "max", "count"])
    assert ranges.loc[("A", "A")].tolist() == [103, 187, 85]
    assert ranges.loc[("B", "B")].tolist() == [1, 11, 11]
    assert ranges.loc[("B", "C")].tolist() == [51, 61, 11]
    assert (frame["haddock_resid"] == frame["pdb_resnum"]).all()
    assert frame.loc[frame["chain_key"] == "A", "reference_index"].iloc[0] == 3


def test_partner_renumbered_when_author_numbers_collide(sync, tmp_path):
    sync.renumber_chain("B", in_place=True)
    sync.renumber_chain("C", in_place=True)
    frame = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path).residue_table
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
    apart = sync.subset_chains(["B", "C"])
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
    assert pair.restraint_table[["segid_a", "segid_b"]].iloc[0].tolist() == ["A", "B"]
    duplex = sync.subset_chains(["B", "C"])
    assert duplex["B"].chain_id == "A"
    keyed = sync.prepare_haddock_docking("A", (duplex, ["B", "C"]), workdir=tmp_path / "keyed")
    assert keyed._locate_chain("A") == (0, "A")


def test_patch_pairs(setup):
    setup.add_patch_pair_restraints([(("A", [20, 23]), ("B", [3, 4])), (("A", [50]), ("C", [2]))])
    table = setup.restraint_table
    assert table["kind"].eq("ambiguous").all() and table["origin"].eq("patch_pair").all()
    assert table["residues_a"].tolist() == [[120, 123], [150]]
    assert table["residues_b"].tolist() == [[3, 4], [52]]
    setup.add_patch_pair_restraints([(("A", [20]), ("B", [3]))], kind="unambiguous")
    assert setup.restraint_table["kind"].iloc[-1] == "unambiguous"
    with pytest.raises(ValueError):
        setup.add_patch_pair_restraints([])


def test_cluster_patches():
    points = [(1, 1), (2, 3), (10, 10), (11, 14), (30, 1)]
    assert cluster_contact_patches(points, gap=4) == [[0, 1], [2, 3], [4]]
    assert cluster_contact_patches(points, gap=0) == [[0], [1], [2], [3], [4]]


def test_contact_restraints(setup, sync):
    contacts = sync.contact_table(contact_type="inter", chain_ids=["A", "B", "C"])
    protein = contacts[(contacts["auth_chain_1"] == "A") | (contacts["auth_chain_2"] == "A")]
    expected = set()
    for row in protein.itertuples():
        side = 1 if row.auth_chain_1 == "A" else 2
        expected.add(getattr(row, f"resnum_{side}"))
    setup.add_contact_restraints(kind="ambiguous")
    table = setup.restraint_table
    assert set().union(*table["residues_a"]) == expected
    setup.clear_restraints().add_contact_restraints(kind="unambiguous", scheme="sidechain", cutoff=0.45)
    patches = setup.contact_patch_table
    assert len(patches) == len(setup.restraint_table) and patches["n_pairs"].sum() > len(patches)
    assert set(patches["chain_b"]) == {"B", "C"}


def test_contacts_transfer_through_reference_positions(sync, tmp_path):
    direct = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "direct")
    direct.add_contact_restraints(kind="ambiguous")
    apo = SyncedStructure.from_path(CIF)
    apo.renumber_chain("A", in_place=True)
    moved = apo.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "moved")
    with pytest.raises(ValueError):
        HaddockDocking([(apo, "A"), (sync, ["B", "C"])], tmp_path / "mixed").add_contact_restraints()
    moved.add_contact_restraints(kind="ambiguous", source=sync, chain_map={"A": "A", "B": "B", "C": "C"})
    shift = lambda residues: [resid - 102 for resid in residues]
    assert moved.restraint_table["residues_a"].tolist() == [shift(r) for r in direct.restraint_table["residues_a"]]
    assert moved.restraint_table["residues_b"].tolist() == direct.restraint_table["residues_b"].tolist()


def _validate(text):
    librestraints = pytest.importorskip("haddock.libs.librestraints")
    librestraints.validate_tbldata(text)


# configuration

def test_preset_modules_seed_and_restraints(sync, tmp_path):
    for method in HADDOCK_PRESETS:
        setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / method, method=method, seed=42)
        if method == "eval":
            setup.add_patch_pair_restraints([(("A", [20, 23]), ("B", [3, 4]))])
        config = _config(setup)
        modules = [name for name, _ in HADDOCK_PRESETS[method]] + ["caprieval"]
        assert [key for key in config if isinstance(config[key], dict)] == modules
        for name in modules:
            if name in HADDOCK_SEEDED_MODULES:
                assert config[name]["iniseed"] == 42
        assert config["caprieval"]["receptor_chain"] == "A" and config["caprieval"]["ligand_chains"] == ["B"]
        assert len(config["molecules"]) == 2 and config["run_dir"].endswith("run")
    abinitio = _config(sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "a"))
    assert abinitio["rigidbody"] == {"sampling": 10000, "cmrest": True, "ranair": True, "iniseed": 917}
    assert abinitio["seletop"]["select"] == 200 and abinitio["emref"]["contactairs"] is True
    assert abinitio["seletopclusts"]["top_clusters"] == 10 and "clustfcc" in abinitio
    oracle = _config(sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "o", method="oracle"))
    assert oracle["rigidbody"]["sampling"] == 200 and oracle["seletop"]["select"] == 10
    assert oracle["rigidbody"]["unambig_fname"] == oracle["emref"]["unambig_fname"]
    assert "ambig_fname" not in oracle["rigidbody"]


def test_eval_uses_patch_pairs_as_ambiguous(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="eval")
    with pytest.raises(ValueError):
        setup.config_toml()
    setup.add_patch_pair_restraints([(("A", [20, 23]), ("B", [3, 4]))])
    config = _config(setup)
    assert config["rigidbody"]["ambig_fname"].endswith("ambig.tbl") and "unambig_fname" not in config["rigidbody"]
    no_reference = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="eval", reference=None)
    no_reference.add_patch_pair_restraints([(("A", [20]), ("B", [3]))])
    with pytest.raises(ValueError):
        no_reference.config_toml()


def test_overrides_and_guards(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="ab-initio", ncores=2,
                                 module_params=LOW_THROUGHPUT, global_params={"debug": True}, caprieval=False)
    setup.add_ambiguous_restraints("A", [20], "B", [3])
    config = _config(setup)
    assert config["rigidbody"]["sampling"] == 4 and config["seletop"]["select"] == 2
    assert config["rigidbody"]["ranair"] is False and config["rigidbody"]["cmrest"] is True
    assert config["emref"]["ambig_fname"] == config["rigidbody"]["ambig_fname"]
    assert config["debug"] is True and config["ncores"] == 2 and "caprieval" not in config
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


def test_docked_without_caprieval(finished, monkeypatch):
    monkeypatch.setattr(HaddockDocking, "_read_io_rows", lambda self, folder: None)
    shutil.rmtree(finished.run_dir / "2_caprieval")
    (finished.run_dir / "1_rigidbody" / "io.json").write_text("{}")
    docked = finished.docked_model_table()
    assert docked["model"].tolist() == ["rigidbody_2.pdb", "rigidbody_1.pdb"]
    assert finished.docked_model_paths() == docked["path"].tolist()


# haddock3 runs (low throughput)

needs_haddock = pytest.mark.skipif(shutil.which("haddock3") is None, reason="haddock3 is not installed")


@pytest.mark.haddock
@needs_haddock
def test_haddock_validates_every_method(sync, tmp_path):
    for method in HADDOCK_PRESETS:
        setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / method, method=method, ncores=2)
        if method == "eval":
            setup.add_patch_pair_restraints([(("A", [20, 23, 26]), ("B", [3, 4, 5]))])
        setup.validate_config()
        assert not (setup.workdir / "setup_check").exists()
    bad = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path / "bad", module_params={"rigidbody": {"samplin": 4}})
    with pytest.raises(RuntimeError):
        bad.validate_config()


@pytest.mark.haddock
@needs_haddock
def test_haddock_oracle_run(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="oracle", seed=7, ncores=2,
                                 module_params=LOW_THROUGHPUT)
    setup.add_contact_restraints(kind="unambiguous", scheme="sidechain", cutoff=0.45)
    setup.run_haddock()
    models = setup.caprieval_tables()["models"]
    assert len(models) == 2 and {"irmsd", "fnat", "dockq"} <= set(models)
    assert models["reference"].eq("reference.pdb").all()
    assert all(path.is_file() for path in setup.docked_model_paths())
    rigid, emref = setup.score_table("rigidbody"), setup.score_table("emref")
    assert len(rigid) == 4 and len(emref) == 4 and emref["origin"].str.startswith("rigidbody_").all()
    assert setup.score_table()["model"].tolist() == [path.name.removesuffix(".gz") for path in setup.docked_model_paths()]
    docked_table = setup.docked_model_table()
    assert docked_table["path"].tolist() == setup.docked_model_paths() and {"seed", "origin", "dockq"} <= set(docked_table)
    assert len(setup.load_docked_models(top=2)) == 2
    analysed = setup.run_caprieval(["rigidbody"])["1_rigidbody"]["models"]
    assert len(analysed) == 4 and analysed["path"].map(Path.is_file).all()
    assert setup.caprieval_tables("rigidbody")["models"]["reference"].eq("reference.pdb").all()
    docked = setup.load_docked_model(0)
    assert sorted(docked.chain_ids) == ["A", "B", "C"]
    assert docked.chain_table()["mapping_fraction"].eq(1.0).all()
    setup.run_haddock(overwrite=True)


@pytest.mark.haddock
@needs_haddock
def test_haddock_abinitio_run(sync, tmp_path):
    setup = sync.prepare_haddock_docking("A", ["B", "C"], workdir=tmp_path, method="ab-initio", seed=3, ncores=2,
                                 module_params={"rigidbody": {"sampling": 8}, "seletop": {"select": 4},
                                                "clustfcc": {"min_population": 1}})
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
