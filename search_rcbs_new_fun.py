import requests
from typing import List, Optional, Dict, Any

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
    """Search RCSB PDB using the Search API v2.

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
        ValueError: If no search criteria are provided.
        requests.HTTPError: On non-2xx API response.
    """
    def t(attr, op, val, neg=False):
        return {"type": "terminal", "service": "text", "parameters": {"attribute": attr, "operator": op, "negation": neg, "value": val}}
    def ft(val):
        return {"type": "terminal", "service": "full_text", "parameters": {"value": val}}
    def g(nodes, op="or"):
        return nodes[0] if len(nodes) == 1 else {"type": "group", "logical_operator": op, "nodes": nodes}

    nodes = []
    if organisms:
        nodes.append(g([t("rcsb_entity_source_organism.ncbi_scientific_name", "contains_phrase", o) for o in organisms]))
    if keywords:
        nodes.append(g([ft(k) for k in keywords]))
    if title_contains:
        nodes.append(g([t("struct.title", "contains_words", x) for x in title_contains]))
    if keywords_field:
        nodes.append(g([t("struct_keywords.pdbx_keywords", "contains_words", x) for x in keywords_field]))
    if experimental_method:
        nodes.append(g([t("exptl.method", "exact_match", m) for m in experimental_method]))
    if resolution_min is not None or resolution_max is not None:
        nodes.append(t("rcsb_entry_info.resolution_combined", "range", {
            "from": resolution_min if resolution_min is not None else 0.0,
            "to": resolution_max if resolution_max is not None else 100.0,
            "include_lower": True, "include_upper": True}))
    if release_date_from or release_date_to:
        dv = {"include_lower": True, "include_upper": True}
        if release_date_from: dv["from"] = release_date_from
        if release_date_to: dv["to"] = release_date_to
        nodes.append(t("rcsb_accession_info.initial_release_date", "range", dv))
    if polymer_entity_count_min is not None or polymer_entity_count_max is not None:
        nodes.append(t("rcsb_entry_info.polymer_entity_count", "range", {
            "from": polymer_entity_count_min if polymer_entity_count_min is not None else 0,
            "to": polymer_entity_count_max if polymer_entity_count_max is not None else 9999,
            "include_lower": True, "include_upper": True}))
    if uniprot_ids:
        nodes.append(g([t("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "exact_match", u) for u in uniprot_ids]))
    if enzyme_classification:
        nodes.append(g([t("rcsb_polymer_entity.rcsb_enzyme_class_combined.ec", "exact_match", e) for e in enzyme_classification]))
    if not nodes:
        raise ValueError("At least one search criterion must be provided")

    query = g(nodes, "and")
    content = ["experimental"] + (["computational"] if include_computational else [])
    opts: Dict[str, Any] = {"results_content_type": content, "sort": [{"sort_by": sort_by, "direction": sort_direction}], "scoring_strategy": "combined"}
    if return_all_hits:
        opts["return_all_hits"] = True
    else:
        opts["paginate"] = {"start": 0, "rows": max_rows}

    r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json={"query": query, "return_type": return_type, "request_options": opts}, timeout=timeout)
    r.raise_for_status()
    return [h["identifier"] for h in r.json().get("result_set", [])]