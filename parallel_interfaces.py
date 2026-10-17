#!/usr/bin/env python
"""Parallel extraction of RCSB interfaces for UniProtKB entries using SyncedStructure."""
###################################### dependencies
import argparse, datetime, os, pickle, re, resource, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import pandas as pd, requests, yaml
from tqdm import tqdm
from SyncedStructure import SyncedStructure
###################################### constants
SCRIPT = Path(__file__).stem
ROOT = Path(__file__).resolve().parent
YAML_PATH = ROOT / "configs" / "param_configs.yaml"
UNIPROT_RE = r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_CIF = "https://files.rcsb.org/download/{}.cif"
###################################### configuration
P = {"outdir": str(ROOT / "data" / "interfaces"), "rows": 250, "workers": len(os.sched_getaffinity(0))}
###################################### helpers
def read_ids(paths):
    ids = [l.strip().upper() for p in paths for l in Path(p).read_text().splitlines() if l.strip()]
    valid = list(dict.fromkeys(i for i in ids if re.match(UNIPROT_RE, i)))
    for i in dict.fromkeys(ids):
        if i not in valid: print(f"[PYTHON-INFO] invalid UniProtKB id skipped: {i}")
    return valid

def fetch_fasta(uid, out):
    r = requests.get(UNIPROT_FASTA.format(uid), timeout=60)
    r.raise_for_status()
    (out / "uniprot" / f"{uid}.fasta").write_text(r.text)
    return "".join(r.text.splitlines()[1:])

def search_pdb(uid):
    q = {"query": {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "in", "value": [uid]}}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": P["rows"]}}}
    r = requests.post(RCSB_SEARCH, json=q, timeout=60)
    return [x["identifier"] for x in r.json().get("result_set", [])] if r.status_code == 200 else []

def consensus(seqs):
    return "".join("1" if "1" in c else "X" if not any(x in "01" for x in c) else "0" for c in zip(*seqs))

def process(task):
    uid, pdb_id, seq, out = task
    raw = requests.get(RCSB_CIF.format(pdb_id), timeout=120).text
    sync = SyncedStructure.from_rcsb_id(pdb_id)
    summ = sync.summary()
    hit = summ[summ["external_ref"] == uid]
    if hit.empty: return None
    focus = hit["auth_chain_id"].values[0]
    partners = sync._chain_contact_partners()
    if not partners.get(focus): return None
    keep = set(partners[focus]) | {focus}
    for chain in set(summ["auth_chain_id"]) - keep: sync.remove_chain(chain)
    sync.set_reference_sequence(chain_id=focus, sequence=seq, remap=True)
    sync.remap_chain(focus)
    cts = sync.compute_contacts()
    contacts = cts[(cts["auth_chain_1"] == focus) | (cts["auth_chain_2"] == focus)]
    interface = sync.get_reference_contact_sequence(focus)
    base = f"{uid}_{pdb_id}_{focus}"
    sync.summary().assign(pdb_id=pdb_id).to_csv(out / "interface_models" / f"{uid}_{pdb_id}_chains.tsv", sep="\t", index=False)
    (out / "interface_models" / f"{base}_model.cif").write_text(raw)
    with open(out / "interface_models" / f"{base}_model.pkl", "wb") as f: pickle.dump(sync, f)
    contacts.to_csv(out / "contacts" / f"{base}_contacts.tsv", sep="\t", index=False)
    return {"uniprot_id": uid, "pdb_id": pdb_id, "chain_id": focus, "neighbours": ",".join(sorted(keep)), "n_contacts": len(contacts), "interface_sequence": interface}
###################################### main
def main():
    start = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="*")
    ap.add_argument("--outdir", default=P["outdir"])
    ap.add_argument("--workers", type=int, default=P["workers"])
    ap.add_argument("--yaml-in", action="store_true")
    ap.add_argument("--yaml-out", action="store_true")
    a = ap.parse_args()
    if a.yaml_in: a.__dict__.update(yaml.safe_load(YAML_PATH.read_text())[SCRIPT])
    if a.yaml_out:
        cfg = (yaml.safe_load(YAML_PATH.read_text()) if YAML_PATH.exists() else None) or {}
        cfg[SCRIPT] = {"inputs": a.inputs, "outdir": a.outdir, "workers": a.workers}
        YAML_PATH.parent.mkdir(parents=True, exist_ok=True)
        YAML_PATH.write_text(yaml.safe_dump(cfg))
    out = Path(a.outdir)
    for d in ("interface_models", "contacts", "uniprot"): (out / d).mkdir(parents=True, exist_ok=True)
    seqs, tasks = {}, []
    for uid in read_ids(a.inputs):
        try: seqs[uid] = fetch_fasta(uid, out)
        except Exception as e:
            print(f"[PYTHON-INFO] {uid} fasta failed: {e}")
            continue
        tasks += [(uid, pdb_id, seqs[uid], out) for pdb_id in search_pdb(uid)]
    rows = {uid: [] for uid in seqs}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(process, t): t[:2] for t in tasks}
        for f in tqdm(as_completed(futs), total=len(futs), desc="[PYTHON-INFO] structures"):
            uid, pdb_id = futs[f]
            try:
                r = f.result()
                if r: rows[uid].append(r)
            except Exception as e: print(f"[PYTHON-INFO] {uid} {pdb_id} skipped: {e}")
    for uid, rs in rows.items():
        if not rs:
            print(f"[PYTHON-INFO] {uid}: no interfaces found")
            continue
        rs.append({"uniprot_id": uid, "pdb_id": "CONSENSUS", "interface_sequence": consensus([r["interface_sequence"] for r in rs])})
        pd.DataFrame(rs).to_csv(out / f"{uid}_interfaces.tsv", sep="\t", index=False)
        print(f"[PYTHON-INFO] {uid}: {len(rs) - 1} interfaces saved")
    mem = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024
    print(f"[PYTHON-INFO] {datetime.datetime.now():%Y-%m-%d %H:%M:%S} | execution {time.time() - start:.1f}s | max memory {mem:.1f} MB | workers {a.workers}")

if __name__ == "__main__":
    main()


