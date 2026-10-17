#!/usr/bin/env python
"""Pipelined RCSB interface extraction with serial FASTA and PDB search, threaded mmCIF downloads, incremental per-ID collection, double-guarded file handling, resource monitoring and a final pf_chart.png."""
###################################### dependencies
import argparse, datetime, os, pickle, re, resource, sys, threading, time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, psutil, requests, yaml
from tqdm import tqdm
from synced_structure import SyncedStructure
###################################### constants
SCRIPT = Path(__file__).stem
ROOT = Path(__file__).resolve().parent
YAML_PATH = ROOT / "configs" / "param_configs.yaml"
UNIPROT_RE = r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_CIF = "https://files.rcsb.org/download/{}.cif"
TAG = "[PYTHON-INFO]"
CHART = "pf_chart.png"
###################################### configuration
NPROC = len(os.sched_getaffinity(0))
TOTAL_GB = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30
P = {"outdir": "data/interfaces", "rows": 250, "max_pdbs": 20, "timeout": 120, "mem_per_worker_gb": 2.5, "workers": max(1, min(NPROC, int(TOTAL_GB // 2.5))), "max_tasks_per_child": 4, "io_threads": 4, "sample_s": 1.0}
DIRS = {"fasta": "uniprotkb_fasta", "pdb_ids": "pdb_ids", "cif": "raw_cifs", "pkl": "model_pkl", "chains": "model_tsv", "contacts": "contacts_tsv", "interfaces": "interface_tsv"}
BAR = dict(file=sys.stderr, leave=True, dynamic_ncols=False, ncols=100, mininterval=0.5)
###################################### helpers
def ok_file(p):
    return p.exists() and p.stat().st_size > 0

def atomic_write(p, text):
    tmp = p.with_name(p.name + ".part")
    tmp.write_text(text)
    tmp.replace(p)

def read_ids(paths):
    ids = [l.strip().upper() for p in paths for l in Path(p).read_text().splitlines() if l.strip()]
    valid = list(dict.fromkeys(i for i in ids if re.match(UNIPROT_RE, i)))
    for i in dict.fromkeys(i for i in ids if i not in valid): tqdm.write(f"{TAG} invalid UniProtKB id skipped: {i}", file=sys.stderr)
    return valid

def fetch_fasta(uid, d, s):
    f = d["fasta"] / f"{uid}.fasta"
    if not ok_file(f):
        r = s.get(UNIPROT_FASTA.format(uid), timeout=60)
        r.raise_for_status()
        atomic_write(f, r.text)
    seq = "".join(f.read_text().splitlines()[1:])
    if not seq: raise ValueError("empty FASTA sequence")
    return seq

def search_pdb(uid, d, s):
    f = d["pdb_ids"] / f"{uid}.txt"
    if not ok_file(f):
        q = {"query": {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "in", "value": [uid]}}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": P["rows"]}}}
        r = s.post(RCSB_SEARCH, json=q, timeout=60)
        pdbs = [x["identifier"] for x in r.json().get("result_set", [])] if r.status_code == 200 else []
        atomic_write(f, "\n".join(pdbs[:P["max_pdbs"]]))
    return f.read_text().split()[:P["max_pdbs"]]

def fetch_cif(uid, pdb_id, d, s):
    f = d["cif"] / f"{uid}_{pdb_id}.cif"
    if not ok_file(f):
        r = s.get(RCSB_CIF.format(pdb_id), timeout=P["timeout"])
        r.raise_for_status()
        atomic_write(f, r.text)
    return f

def is_done(uid, d):
    return ok_file(d["interfaces"] / f"{uid}_interfaces.tsv")

def consensus(seqs):
    return "".join("1" if "1" in c else "X" if not any(x in "01" for x in c) else "0" for c in zip(*seqs))

def process(task):
    uid, pdb_id, seq, d = task
    cif = d["cif"] / f"{uid}_{pdb_id}.cif"
    if not ok_file(cif): raise FileNotFoundError(cif.name)
    sync = SyncedStructure.from_path(cif)
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
    atomic_write(d["chains"] / f"{uid}_{pdb_id}_chains.tsv", sync.summary().assign(pdb_id=pdb_id).to_csv(sep="\t", index=False))
    atomic_write(d["contacts"] / f"{base}_contacts.tsv", contacts.to_csv(sep="\t", index=False))
    pkl = d["pkl"] / f"{base}_model.pkl"
    tmp = pkl.with_name(pkl.name + ".part")
    with open(tmp, "wb") as f: pickle.dump(sync, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(pkl)
    return {"uniprot_id": uid, "pdb_id": pdb_id, "chain_id": focus, "neighbours": ",".join(sorted(keep)), "n_contacts": len(contacts), "interface_sequence": interface}

def write_uid(uid, rows, d):
    out = d["interfaces"] / f"{uid}_interfaces.tsv"
    rows.append({"uniprot_id": uid, "pdb_id": "CONSENSUS", "interface_sequence": consensus([r["interface_sequence"] for r in rows])})
    atomic_write(out, pd.DataFrame(rows).to_csv(sep="\t", index=False))

def finish(uid, rows, d):
    if rows:
        write_uid(uid, rows, d)
        tqdm.write(f"{TAG} {uid}: {len(rows)} interfaces saved", file=sys.stderr)
    else: tqdm.write(f"{TAG} {uid}: no interfaces found", file=sys.stderr)

def dispatch(ex, waiting, futs, remaining, seqs, d, bar_u, block=False):
    for uid in [u for u, (cf, _) in waiting.items() if block or all(f.done() for f in cf)]:
        cf, pdbs = waiting.pop(uid)
        wait(cf)
        if any(f.exception() for f in cf):
            tqdm.write(f"{TAG} {uid} mmCIF download failed, skipped", file=sys.stderr)
            bar_u.update(1)
            continue
        good = [p for p in pdbs if ok_file(d["cif"] / f"{uid}_{p}.cif")]
        if len(good) < len(pdbs): tqdm.write(f"{TAG} {uid}: {len(pdbs) - len(good)} mmCIF missing on disk, skipped", file=sys.stderr)
        if not good:
            bar_u.update(1)
            continue
        for p in good:
            remaining[uid] += 1
            futs[ex.submit(process, (uid, p, seqs[uid], d))] = (uid, p)

def collect(futs, rows, remaining, d, bar_u, block=False):
    if block: wait(list(futs))
    for f in [f for f in list(futs) if f.done()]:
        uid, p = futs.pop(f)
        try:
            r = f.result()
            if r: rows[uid].append(r)
        except Exception as e: tqdm.write(f"{TAG} {uid} {p} failed: {e}", file=sys.stderr)
        remaining[uid] -= 1
        if remaining[uid] == 0:
            finish(uid, rows.pop(uid, []), d)
            bar_u.update(1)

def monitor(stop, log, interval):
    me = psutil.Process()
    psutil.cpu_percent(None)
    while not stop.wait(interval):
        try:
            kids = me.children(recursive=True)
            rss = me.memory_info().rss
            for k in kids:
                try: rss += k.memory_info().rss
                except psutil.Error: pass
            log.append((time.time(), psutil.cpu_percent(None), rss / 2**30, len(kids)))
        except psutil.Error: pass

def plot_monitor(log, start, path):
    if not log: return
    t = [x[0] - start for x in log]
    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(9, 8))
    ax[0].plot(t, [x[1] for x in log], color="tab:blue"); ax[0].set_ylabel("CPU %"); ax[0].set_ylim(0, 100)
    ax[1].plot(t, [x[2] for x in log], color="tab:orange"); ax[1].set_ylabel("RSS (GB)")
    ax[2].step(t, [x[3] for x in log], where="post", color="tab:green"); ax[2].set_ylabel("workers"); ax[2].set_xlabel("elapsed (s)")
    for a_ in ax: a_.grid(alpha=0.3)
    fig.suptitle("Pipeline resource monitor")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)
###################################### main
def main():
    start = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", nargs="+", default=[])
    ap.add_argument("--outdir", default=P["outdir"])
    ap.add_argument("--workers", type=int, default=P["workers"])
    ap.add_argument("--io-threads", type=int, default=P["io_threads"])
    ap.add_argument("--yaml-in", nargs="?", const=str(YAML_PATH), default=None)
    ap.add_argument("--yaml-out", nargs="?", const=str(YAML_PATH), default=None)
    a = ap.parse_args()
    if a.yaml_in: a.__dict__.update(yaml.safe_load(Path(a.yaml_in).read_text())[SCRIPT])
    if a.yaml_out:
        yp = Path(a.yaml_out)
        cfg = (yaml.safe_load(yp.read_text()) if yp.exists() else None) or {}
        cfg[SCRIPT] = {"input": a.input, "outdir": a.outdir, "workers": a.workers, "io_threads": a.io_threads}
        yp.parent.mkdir(parents=True, exist_ok=True)
        yp.write_text(yaml.safe_dump(cfg))
    out = Path(a.outdir)
    d = {k: out / v for k, v in DIRS.items()}
    for p in d.values(): p.mkdir(parents=True, exist_ok=True)
    W = max(1, min(a.workers, NPROC))
    tqdm.write(f"{TAG} nproc={NPROC} ram={TOTAL_GB:.1f}GB workers={W} io_threads={a.io_threads} max_pdbs={P['max_pdbs']} outdir={out.resolve()}", file=sys.stderr)
    ids = read_ids(a.input)
    s = requests.Session()
    seqs, waiting, futs, remaining, rows = {}, {}, {}, defaultdict(int), defaultdict(list)
    stop, log = threading.Event(), []
    mon = threading.Thread(target=monitor, args=(stop, log, P["sample_s"]), daemon=True)
    mon.start()
    bar_u = tqdm(total=len(ids), desc="UniProtKB entries", position=0, unit="id", **BAR)
    bar_f = tqdm(total=len(ids), desc="FASTA downloaded ", position=1, unit="fasta", **BAR)
    bar_c = tqdm(total=0, desc="mmCIF downloaded  ", position=2, unit="cif", **BAR)
    with ProcessPoolExecutor(max_workers=W, max_tasks_per_child=P["max_tasks_per_child"]) as ex, ThreadPoolExecutor(max_workers=a.io_threads) as io:
        for uid in ids:
            if is_done(uid, d):
                bar_f.update(1); bar_u.update(1); continue
            try:
                seqs[uid] = fetch_fasta(uid, d, s)
                bar_f.update(1)
            except Exception as e:
                tqdm.write(f"{TAG} {uid} FASTA failed: {e}", file=sys.stderr); bar_u.update(1); continue
            try:
                pdbs = search_pdb(uid, d, s)
            except Exception as e:
                tqdm.write(f"{TAG} {uid} RCSB search failed: {e}", file=sys.stderr); bar_u.update(1); continue
            if not pdbs:
                tqdm.write(f"{TAG} {uid}: no PDB entries", file=sys.stderr); bar_u.update(1); continue
            bar_c.total += len(pdbs); bar_c.refresh()
            cf = []
            for p in pdbs:
                f = io.submit(fetch_cif, uid, p, d, s)
                f.add_done_callback(lambda _: bar_c.update(1))
                cf.append(f)
            waiting[uid] = (cf, pdbs)
            dispatch(ex, waiting, futs, remaining, seqs, d, bar_u)
            collect(futs, rows, remaining, d, bar_u)
        dispatch(ex, waiting, futs, remaining, seqs, d, bar_u, block=True)
        collect(futs, rows, remaining, d, bar_u, block=True)
    bar_u.close(); bar_f.close(); bar_c.close()
    stop.set(); mon.join()
    plot_monitor(log, start, out / CHART)
    tqdm.write(f"{TAG} resource chart saved to {(out / CHART).resolve()}", file=sys.stderr)
    mem = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024
    tqdm.write(f"{TAG} {datetime.datetime.now():%Y-%m-%d %H:%M:%S} | execution {time.time() - start:.1f}s | max memory {mem:.1f} MB | workers {W}", file=sys.stderr)

if __name__ == "__main__":
    main()


