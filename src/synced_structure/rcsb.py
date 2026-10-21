"""RCSB PDB queries: attribute and full-text searches via the Search API v2, and citation lookups via the Data API."""

###################################### IMPORTS AND CONSTANTS

from typing import Any, Dict, List, Optional

import pandas as pd
import requests

from .constants import PUB_COLUMNS, RCSB_GRAPHQL_URL, RCSB_RETURN_TYPES, RCSB_SEARCH_URL

PUB_QUERY = """query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    citation {
      id rcsb_is_primary title rcsb_authors journal_abbrev journal_volume
      page_first page_last year pdbx_database_id_DOI pdbx_database_id_PubMed
    }
  }
}"""

###################################### QUERY NODES


def _text(attr, op, val, neg=False):
    """Return a text-service terminal node."""
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attr, "operator": op, "negation": neg, "value": val}}


def _full_text(val):
    """Return a full_text-service terminal node."""
    return {"type": "terminal", "service": "full_text", "parameters": {"value": val}}


def _group(nodes, op="or"):
    """Return a single node as is, or wrap several nodes in a logical group."""
    return nodes[0] if len(nodes) == 1 else {"type": "group", "logical_operator": op, "nodes": nodes}


###################################### SEARCH


def search_rcsb(
    organisms: Optional[List[str]] = None,
    keywords: Optional[List[str]] = None,
    title_contains: Optional[List[str]] = None,
    keywords_field: Optional[List[str]] = None,
    experimental_method: Optional[List[str]] = None,
    resolution_max: Optional[float] = None,
    resolution_min: Optional[float] = None,
    release_date_from: Optional[str] = None,
    release_date_to: Optional[str] = None,
    polymer_entity_count_min: Optional[int] = None,
    polymer_entity_count_max: Optional[int] = None,
    uniprot_ids: Optional[List[str]] = None,
    enzyme_classification: Optional[List[str]] = None,
    return_type: str = "entry",
    max_rows: int = 10000,
    return_all_hits: bool = False,
    include_computational: bool = False,
    sort_by: str = "score",
    sort_direction: str = "desc",
    timeout: int = 30,
) -> List[str]:
    """Search RCSB PDB using the Search API v2; values within one criterion are OR-ed, criteria are AND-ed.

    Args:
        organisms: Scientific names matched with contains_phrase on
            rcsb_entity_source_organism.ncbi_scientific_name (OR).
        keywords: Free-text terms searched via full_text service (OR).
        title_contains: Words that must appear in struct.title (OR).
        keywords_field: Words that must appear in struct_keywords.pdbx_keywords (OR).
        experimental_method: Exact matches on exptl.method (OR).
        resolution_max: Upper resolution limit in Å.
        resolution_min: Lower resolution limit in Å.
        release_date_from: Earliest release date (YYYY-MM-DD).
        release_date_to: Latest release date (YYYY-MM-DD).
        polymer_entity_count_min: Minimum number of polymer entities.
        polymer_entity_count_max: Maximum number of polymer entities.
        uniprot_ids: UniProt accessions (OR).
        enzyme_classification: EC numbers (OR).
        return_type: One of entry, assembly, polymer_entity, non_polymer_entity,
            polymer_instance, mol_definition.
        max_rows: Maximum results to return (ignored if return_all_hits=True).
        return_all_hits: Return every matching identifier.
        include_computational: Also include Computed Structure Models.
        sort_by: Sort field (score, size, code, date, ...).
        sort_direction: asc or desc.
        timeout: Request timeout in seconds.

    Returns:
        List of matching identifiers.

    Raises:
        ValueError: If no search criteria are provided or return_type is unknown.
        requests.HTTPError: On non-2xx API response.
    """
    if return_type not in RCSB_RETURN_TYPES:
        raise ValueError(f"return_type must be one of {RCSB_RETURN_TYPES}, got {return_type!r}")
    nodes = []
    if organisms:
        nodes.append(_group([_text("rcsb_entity_source_organism.ncbi_scientific_name", "contains_phrase", o) for o in organisms]))
    if keywords:
        nodes.append(_group([_full_text(k) for k in keywords]))
    if title_contains:
        nodes.append(_group([_text("struct.title", "contains_words", x) for x in title_contains]))
    if keywords_field:
        nodes.append(_group([_text("struct_keywords.pdbx_keywords", "contains_words", x) for x in keywords_field]))
    if experimental_method:
        nodes.append(_group([_text("exptl.method", "exact_match", m) for m in experimental_method]))
    if resolution_min is not None or resolution_max is not None:
        nodes.append(_text("rcsb_entry_info.resolution_combined", "range", {
            "from": resolution_min if resolution_min is not None else 0.0,
            "to": resolution_max if resolution_max is not None else 100.0,
            "include_lower": True, "include_upper": True}))
    if release_date_from or release_date_to:
        dv = {"include_lower": True, "include_upper": True}
        if release_date_from:
            dv["from"] = release_date_from
        if release_date_to:
            dv["to"] = release_date_to
        nodes.append(_text("rcsb_accession_info.initial_release_date", "range", dv))
    if polymer_entity_count_min is not None or polymer_entity_count_max is not None:
        nodes.append(_text("rcsb_entry_info.polymer_entity_count", "range", {
            "from": polymer_entity_count_min if polymer_entity_count_min is not None else 0,
            "to": polymer_entity_count_max if polymer_entity_count_max is not None else 9999,
            "include_lower": True, "include_upper": True}))
    if uniprot_ids:
        nodes.append(_group([
            _text("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "exact_match", u)
            for u in uniprot_ids
        ]))
    if enzyme_classification:
        nodes.append(_group([_text("rcsb_polymer_entity.rcsb_enzyme_class_combined.ec", "exact_match", e) for e in enzyme_classification]))
    if not nodes:
        raise ValueError("At least one search criterion must be provided")

    content = ["experimental"] + (["computational"] if include_computational else [])
    opts: Dict[str, Any] = {
        "results_content_type": content,
        "sort": [{"sort_by": sort_by, "direction": sort_direction}],
        "scoring_strategy": "combined",
    }
    if return_all_hits:
        opts["return_all_hits"] = True
    else:
        opts["paginate"] = {"start": 0, "rows": max_rows}

    payload = {"query": _group(nodes, "and"), "return_type": return_type, "request_options": opts}
    response = requests.post(RCSB_SEARCH_URL, json=payload, timeout=timeout)
    if response.status_code == 204:  # valid query with no hits
        return []
    response.raise_for_status()
    return [hit["identifier"] for hit in response.json().get("result_set", [])]


###################################### PUBLICATIONS


def search_rcsb_pubs(
    df: pd.DataFrame,
    id_col: str = "pdb_id",
    primary_only: bool = True,
    merge: bool = True,
    batch_size: int = 200,
    timeout: int = 30,
) -> pd.DataFrame:
    """Fetch the citations of the PDB entries listed in a DataFrame column from the RCSB Data API.

    IDs are stripped and upper-cased; identifiers such as '1ABC-1' or '1ABC_A' are reduced to the
    four-character entry ID. Entries unknown to RCSB get no citation rows.

    Args:
        df: DataFrame holding PDB IDs.
        id_col: Name of the column with PDB IDs.
        primary_only: Keep only the primary citation of each entry (one row per entry).
        merge: Left-merge the citations onto df (on id_col); if False, return only the citation table.
        batch_size: Number of entry IDs per GraphQL request.
        timeout: Request timeout in seconds.

    Returns:
        With merge=False, a citation table with PUB_COLUMNS. With merge=True, df plus every PUB_COLUMNS
        column except 'pdb_id', NaN where an entry has no citation; rows repeat once per citation when
        primary_only=False. Authors are joined by '; '.

    Raises:
        KeyError: If id_col is not in df.
        requests.HTTPError: On non-2xx API response.
    """
    if id_col not in df.columns:
        raise KeyError(f"Column {id_col!r} not found in DataFrame")
    entry_ids = df[id_col].dropna().astype(str).str.strip().str.upper().str[:4]
    unique_ids = list(dict.fromkeys(i for i in entry_ids if len(i) == 4))

    rows = []
    for start in range(0, len(unique_ids), batch_size):
        batch = unique_ids[start:start + batch_size]
        response = requests.post(RCSB_GRAPHQL_URL, json={"query": PUB_QUERY, "variables": {"ids": batch}}, timeout=timeout)
        response.raise_for_status()
        data = response.json()
        if data.get("errors") and not data.get("data"):
            raise requests.HTTPError(f"RCSB GraphQL error: {data['errors']}", response=response)
        for entry in (data.get("data") or {}).get("entries") or []:
            for cit in entry.get("citation") or []:
                is_primary = cit.get("rcsb_is_primary") == "Y"
                if primary_only and not is_primary:
                    continue
                rows.append([
                    entry["rcsb_id"], cit.get("id"), is_primary, cit.get("title"),
                    "; ".join(cit.get("rcsb_authors") or []), cit.get("journal_abbrev"), cit.get("journal_volume"),
                    cit.get("page_first"), cit.get("page_last"), cit.get("year"),
                    cit.get("pdbx_database_id_DOI"), cit.get("pdbx_database_id_PubMed"),
                ])
    pubs = pd.DataFrame(rows, columns=PUB_COLUMNS)
    pubs["year"] = pubs["year"].astype("Int64")
    pubs["pubmed_id"] = pubs["pubmed_id"].astype("Int64")
    if not merge:
        return pubs
    keyed = df.assign(_entry_id=entry_ids.reindex(df.index))
    pubs = pubs.rename(columns={"pdb_id": "_entry_id"})
    return keyed.merge(pubs, on="_entry_id", how="left").drop(columns="_entry_id")
