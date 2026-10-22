"""Regression and coverage tests for behaviour the smoke tests do not reach, all offline."""

import copy
import logging
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import requests
from Bio.PDB import PDBIO, MMCIFParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

matplotlib.use("Agg")

from synced_structure import SyncedStructure, rcsb, side_by_side  # noqa: E402
from synced_structure import utils  # noqa: E402
from synced_structure.residues import residue_class  # noqa: E402
from synced_structure.constants import MOIETY_NUCLEIC, MOIETY_OTHER  # noqa: E402

CIF = Path(__file__).parent / "data" / "1AAY.cif"


@pytest.fixture
def sync():
    """Return a fresh SyncedStructure of 1AAY keyed by author chain IDs."""
    return SyncedStructure.from_path(CIF)


class _Records(logging.Handler):
    """Collect log records from the package logger."""

    def __init__(self):
        super().__init__(level=logging.INFO)
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


@pytest.fixture
def log_records():
    handler = _Records()
    logger = logging.getLogger(utils.LOGGER_NAME)
    logger.addHandler(handler)
    yield handler
    logger.removeHandler(handler)


def _single_atom_residue(resname):
    residue = Residue((" ", 1, " "), resname, " ")
    residue.add(Atom("U", np.zeros(3), 0.0, 1.0, " ", "U", 1, element="U"))
    return residue


# residue classification

def test_single_atom_uranium_is_ion_not_nucleic():
    assert residue_class(_single_atom_residue("U")) == MOIETY_OTHER


def test_multi_atom_uracil_is_still_nucleic():
    residue = Residue((" ", 1, " "), "U", " ")
    for serial, name in enumerate(("P", "C1'", "C2"), start=1):
        residue.add(Atom(name, np.zeros(3), 0.0, 1.0, " ", name, serial, element=name[0]))
    assert residue_class(residue) == MOIETY_NUCLEIC


# interface metrics

def test_energy_per_partner_is_nan_when_contact_cutoff_exceeds_score_cutoff(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("synced_structure:\n  interface_contact_cutoff_a: 9.0\n  interface_score_cutoff_a: 4.0\n")
    loaded = SyncedStructure.from_path(CIF, config_path=config)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        residues = loaded.interface_residues("A", "B")
    no_score = residues["energy_per_partner"].isna()
    assert no_score.any()
    assert (residues.loc[no_score, "n_partners"] > 0).all()


# contact plots

def test_sidechain_scheme_maps_base_pairs_of_dna_duplex(sync):
    # 1AAY strands B and C pair B(i) with C(12 - i); C4' backbone points are farther apart than 0.4 nm
    assert np.isnan(sync.contact_matrix("B", "C", 0.4, scheme="ca")["matrix"]).all()
    rows, cols = np.nonzero(~np.isnan(sync.contact_matrix("B", "C", 0.4, scheme="sidechain")["matrix"]))
    assert len(rows) > 0
    assert set((rows + cols).tolist()) <= {10, 11, 12}


def test_multi_plot_refuses_more_pairs_than_colors(sync):
    for index, new_id in enumerate(["X", "Y", "W"]):
        clone = copy.deepcopy(sync["A"].chain)
        clone.id = new_id
        sync.add_chain(clone, reference_sequence="ACDEFG", entity_id=str(10 + index))
    with pytest.raises(ValueError, match="pair_colors needs 15 colors"):
        sync.plot_contacts_multi(("A", "B", "C", "X", "Y", "W"))


# chain mutation and renaming

def test_add_replace_remove_chain(sync):
    clone = copy.deepcopy(sync["A"].chain)
    clone.id = "Z"
    sync.add_chain(clone, reference_sequence="ACDEFG", entity_id="7", metadata={"note": 1})
    assert "Z" in sync.chain_ids
    assert sync.get_reference_sequence("Z") == "ACDEFG"
    assert sync.chains["Z"].metadata == {"note": 1}
    with pytest.raises(ValueError, match="already exists"):
        sync.add_chain(copy.deepcopy(clone))
    sync.replace_chain("Z", clone, entity_id="8")
    assert sync.chains["Z"].entity_id == "8"
    sync.remove_chain("Z")
    assert "Z" not in sync.chain_ids


def test_rename_chains_keeps_author_map_and_sequence(sync):
    sequence = sync.get_sequence("A")
    sync.rename_chains({"A": "X"})
    assert "X" in sync.chain_ids and "A" not in sync.chain_ids
    assert sync.get_sequence("X") == sequence
    with pytest.raises(ValueError, match="duplicate"):
        sync.rename_chains({"B": "C"})


def test_rename_by_sequence_reports_unmatched(sync):
    report = sync.rename_by_sequence({"P": sync.get_sequence("B"), "Q": "WWWWWWWWWWWWWWW"})
    statuses = dict(zip(report["new_chain_id"].fillna("-"), report["status"]))
    assert statuses["P"] == "renamed"
    assert statuses["Q"] == "unmatched_target"
    assert "P" in sync.chain_ids
    assert (report["status"] == "unmatched_chain").sum() == 2


# loading

def test_label_keys_and_pdb_input_warning(tmp_path, log_records):
    labelled = SyncedStructure.from_path(CIF, key_by="label")
    assert labelled.key_by == "label"
    structure = MMCIFParser(QUIET=True).get_structure("1AAY", str(CIF))
    pdb = tmp_path / "1AAY.pdb"
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(pdb))
    from_pdb = SyncedStructure.from_path(pdb, key_by="label")
    assert sorted(from_pdb.chain_ids) == sorted(SyncedStructure.from_path(CIF).chain_ids)
    assert any("key_by='label' is not available for PDB" in m for m in log_records.messages)


def test_remap_false_logs_missing_mapping(sync, log_records):
    sync.set_reference_sequence("B", "GCGTGGGCGTCC", remap=False)
    assert sync.chains["B"].mapped_residues == 0
    assert any("has no residue mapping until remapped" in m for m in log_records.messages)


def test_log_completion_without_resource_module(monkeypatch, log_records):
    monkeypatch.setattr(utils, "resource", None)
    utils.log_completion(0.0, "done")
    assert any("max memory: n/a" in m for m in log_records.messages)


# deprecated names and dict access

def test_deprecated_method_names_forward(sync):
    with pytest.warns(DeprecationWarning, match="block_contact_matrix"):
        block = sync.complete_contacts("A", "B")
    assert block["intra"].shape == sync.block_contact_matrix("A", "B")["intra"].shape
    with pytest.warns(DeprecationWarning, match="interface_partner_labels"):
        labels = sync.partners("A", "B")
    assert labels == sync.interface_partner_labels("A", "B")
    with pytest.warns(DeprecationWarning, match="show_structure_3d"):
        view = sync.jupy_viz(title="old")
    assert "<div" in side_by_side(view).data


def test_deprecated_from_rcsb_id_forwards(monkeypatch):
    calls = []
    monkeypatch.setattr(SyncedStructure, "from_rcsb", classmethod(lambda cls, *a, **k: calls.append((a, k)) or "ok"))
    with pytest.warns(DeprecationWarning, match="from_rcsb"):
        assert SyncedStructure.from_rcsb_id("1abc", key_by="label") == "ok"
    assert calls == [(("1abc",), {"key_by": "label"})]


def test_dict_style_chain_access_is_deprecated(sync):
    with pytest.warns(DeprecationWarning, match="use .residues"):
        assert sync["A"]["residues"] is sync["A"].residues


# RCSB search and publications, with requests mocked

class _Response:
    def __init__(self, status_code=200, payload=None):
        self.status_code = status_code
        self._payload = payload or {}
        self.sent = None

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")


def test_search_rcsb_builds_query_and_reads_hits(monkeypatch):
    seen = {}

    def fake_post(url, json, timeout):
        seen["url"], seen["json"] = url, json
        return _Response(payload={"result_set": [{"identifier": "1AAY"}, {"identifier": "2ABC"}]})

    monkeypatch.setattr(rcsb.requests, "post", fake_post)
    hits = rcsb.search_rcsb(organisms=["Homo sapiens"], resolution_max=2.0, max_rows=5)
    assert hits == ["1AAY", "2ABC"]
    assert seen["url"] == rcsb.RCSB_SEARCH_URL
    query = seen["json"]["query"]
    assert query["type"] == "group" and query["logical_operator"] == "and"
    assert len(query["nodes"]) == 2
    assert seen["json"]["request_options"]["paginate"] == {"start": 0, "rows": 5}


def test_search_rcsb_no_hits_and_validation(monkeypatch):
    monkeypatch.setattr(rcsb.requests, "post", lambda *a, **k: _Response(status_code=204))
    assert rcsb.search_rcsb(keywords=["kinase"]) == []
    with pytest.raises(ValueError, match="At least one search criterion"):
        rcsb.search_rcsb()
    with pytest.raises(ValueError, match="return_type"):
        rcsb.search_rcsb(keywords=["x"], return_type="nope")


def test_search_rcsb_pubs_keeps_partial_data_and_logs(monkeypatch, log_records):
    payload = {
        "errors": [{"message": "one entry unknown"}],
        "data": {"entries": [{
            "rcsb_id": "1AAY",
            "citation": [{"id": "primary", "rcsb_is_primary": "Y", "title": "T", "rcsb_authors": ["A. Author"],
                          "journal_abbrev": "J", "journal_volume": "1", "page_first": "1", "page_last": "2",
                          "year": 2001, "pdbx_database_id_DOI": "10/x", "pdbx_database_id_PubMed": 123}],
        }]},
    }
    monkeypatch.setattr(rcsb.requests, "post", lambda *a, **k: _Response(payload=payload))
    df = pd.DataFrame({"pdb_id": ["1aay-1"]})
    merged = rcsb.search_rcsb_pubs(df)
    assert merged.loc[0, "title"] == "T"
    assert merged.loc[0, "pubmed_id"] == 123
    assert any("partial errors" in m for m in log_records.messages)


def test_search_rcsb_pubs_missing_column():
    with pytest.raises(KeyError):
        rcsb.search_rcsb_pubs(pd.DataFrame({"x": [1]}))


# shared categorical palette

def test_categorical_colors_come_from_palette():
    from synced_structure.constants import PALETTE
    from synced_structure.plot_style import CONTACT_MAP_COLOR, MULTI_PAIR_COLORS
    from synced_structure.viz3d import CHAIN_PALETTE

    assert len(PALETTE) == 10
    assert tuple(PALETTE) == tuple(CHAIN_PALETTE)
    assert sorted(MULTI_PAIR_COLORS) == sorted(PALETTE)
    assert CONTACT_MAP_COLOR in PALETTE
