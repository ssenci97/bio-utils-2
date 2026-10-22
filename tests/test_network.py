"""Live tests against RCSB services, skipped by default; run with: pytest -m network."""

import pandas as pd
import pytest

from synced_structure import SyncedStructure, search_rcsb, fetch_rcsb_citations

pytestmark = pytest.mark.network


@pytest.fixture(scope="module")
def hemoglobin():
    """Return deoxy human hemoglobin 4HHB: alpha chains A and C, beta chains B and D."""
    return SyncedStructure.from_rcsb("4HHB")


def test_from_rcsb_hemoglobin_chains(hemoglobin):
    assert hemoglobin.chain_ids == ["A", "B", "C", "D"]
    summary = hemoglobin.chain_table()
    descriptions = dict(zip(summary["chain_id"], summary["description"]))
    assert "alpha" in descriptions["A"].lower() and "beta" in descriptions["B"].lower()
    assert summary["mapping_fraction"].eq(1.0).all()


def test_hemoglobin_interfaces(hemoglobin):
    complexes = hemoglobin.interface_tables()["complex"].set_index(["chain_a", "chain_b"])
    assert complexes.loc[("A", "B"), "n_contacts"] > 0
    assert complexes.loc[("B", "D"), "n_contacts"] == 0
    assert str(complexes.loc[("B", "D"), "dASA_A2"]) == "0.0"


def test_search_rcsb_returns_entry_ids():
    hits = search_rcsb(keywords=["hemoglobin"], max_rows=5)
    assert 0 < len(hits) <= 5
    assert all(len(hit) == 4 for hit in hits)


def test_fetch_rcsb_citations_primary_citation():
    pubs = fetch_rcsb_citations(pd.DataFrame({"pdb_id": ["4HHB"]}))
    assert len(pubs) == 1
    assert isinstance(pubs.loc[0, "title"], str) and pubs.loc[0, "title"]
