#!/usr/bin/env python
# coding: utf-8

# # SYNCEDSTRUCTURE
# ### **A Class to handle, crossreference and annotate Structural Proteomics Data**
# 
# Credits:
# LICENCES:

# In[1]:


### <DEPENDENCIES and SyncedStructure Methods display> 
import re, urllib.request, requests, io, copy, os, pickle
import pandas as pd
import seaborn as sns
from synced_structure import SyncedStructure
from viz3d import side_by_side

def get_uniprot(uid):
    u = str(uid).strip().upper()
    if not re.match(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$", u):
        raise ValueError()
    buf = io.StringIO(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{u}.fasta").read().decode())
    return u, "".join(buf.read().splitlines()[1:])


pd.set_option("display.width", 1000)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_rows", 100)

synced_structure_options = dir(SyncedStructure)
d = [x for x in synced_structure_options if x.startswith("__")]
h = [x for x in synced_structure_options
    if x.startswith("_") and not x.startswith("__")]
e = [x for x in synced_structure_options if not x.startswith("_")]
for g in (d, h, e):
  if g:
    [print(g[i : i + 4]) for i in range(0, len(g), 4)]
    print()
### <\DEPENDENCIES and SyncedStructure Methods display>


# In[2]:


### SYMBOLIC ONE LINE CHANGE


# # Human Timeless Homolog Interfaces

# In[3]:


import requests

upid_list = ["Q9UNS1"]
payload = {
    "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "operator": "in",
            "value": upid_list
        }
    },
    "return_type": "entry",
    "request_options": {"paginate": {"start": 0, "rows": 250}}
}
### r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=payload)
### print([x["identifier"] for x in r.json().get("result_set", [])])


# In[4]:


import os, pickle
import pandas as pd
pdb_ids = ['4XHT', '4XHU', '4XHW', '5MQI', '6T9Q', '6TAZ', '7PFO', '7PLO', '8B9D', '9O34']
TIM_INTERFACES_PATH = "tutorial-data/TIM/TIM-interfaces-raw.pkl"

try:
    if os.path.exists(TIM_INTERFACES_PATH):
        with open(TIM_INTERFACES_PATH, "rb") as f:
            synced_list = pickle.load(f)
    else:
        raise FileNotFoundError
except (EOFError, FileNotFoundError):
    synced_list = [SyncedStructure.from_rcsb_id(pid) for pid in pdb_ids]
    os.makedirs("tutorial-data/TIM", exist_ok=True)
    with open(TIM_INTERFACES_PATH, "wb") as f:
        pickle.dump(synced_list, f)

sync_summaries = []
for pdb_id, sync in zip(pdb_ids, synced_list):
    summary=sync.summary()
    summary['pdb_id']=pdb_id
    sync_summaries.append(summary)
sync_summaries = [sync.summary() for sync in synced_list]

sync_df = pd.concat(sync_summaries, ignore_index=True)


# In[5]:


print(sync_df)


# In[6]:


# 1. Populate tim_neighbours with contact partners
tim_neighbours = {}
for pdb_id, sync in zip(pdb_ids, synced_list):
    try:
        summ = sync.summary()
        if summ is None or (hasattr(summ, 'empty') and summ.empty):
            raise ValueError("Empty summary output")

        mask = summ['external_ref'] == upid_list[0]
        matched = summ[mask]
        if matched.empty:
            raise ValueError("Empty output from mask filter")

        chain_id = matched['auth_chain_id'].values[0]
        if chain_id is None or (hasattr(chain_id, '__len__') and len(chain_id) == 0):
            raise ValueError("Empty chain_id")

        partners = sync._chain_contact_partners()
        if not partners or chain_id not in partners:
            raise KeyError("chain_id not found in contact partners")            

        neighbours = list(partners[chain_id])
        if not neighbours:
            # Fixed: added missing 'f' prefix here
            raise ValueError(f"Empty neighbours output for {pdb_id} {chain_id}")

        if chain_id not in neighbours:
            neighbours.append(chain_id)

        tim_neighbours[pdb_id] = neighbours

    except Exception as e:
        print(f"[WARNING] {e}")
        continue


# 2. Remove non-partner chains
for pdb_id, sync in zip(pdb_ids, synced_list):
    if pdb_id not in tim_neighbours:
        continue
    try:
        keep = set(tim_neighbours[pdb_id])
        chain_ids = sync.summary()["auth_chain_id"].tolist()
        to_remove = [chain for chain in chain_ids if chain not in keep]
        if to_remove:
            print(f"{pdb_id}: removing non-partner chains {to_remove}")
            for chain in to_remove:
                sync.remove_chain(chain)

    except Exception as e:
        print(f"[WARNING] Error processing {pdb_id}: {e}")
        continue


# In[7]:


upseq = get_uniprot(upid_list[0])
print(upseq)

viz_list = []
for pdb_id, sync in zip(pdb_ids, synced_list):
    if pdb_id not in tim_neighbours:
        continue
    keep = tim_neighbours[pdb_id]
    focus = keep[-1]
    sync.set_reference_sequence(chain_id=focus, sequence=upseq[1], remap=True)
    sync.remap_chain(focus)
    title = f"Timeless ({focus}) in {pdb_id}"
    PALETTE = dict(zip(keep, ["#ced4da"] * len(keep)))
    PALETTE[focus] = '#a4161a'
    plot = sync.jupy_viz(title=title, colors=PALETTE, size=320)
    viz_list.append(plot)


# In[8]:


side_by_side(*viz_list)


# ## Extracting Interfaces for Chains of interest
# 
# Here, we extract the `contact_sequences` for the chain of interest in each protein complex.

# In[9]:


### <EXAMPLE CONTACTS>
cts = sync.compute_contacts()
def _filter_chain(chain_id, df):
    df_filtered = df[(df['auth_chain_1']==chain_id)|(df['auth_chain_2']==chain_id)]
    return df_filtered

print(f"All: {len(cts)}")
print(f"Filtered: {len(_filter_chain("K",cts))}")
### <EXAMPLE CONTACTS>


# In[12]:


interface_list = []
for pdb_id, sync in zip(pdb_ids, synced_list):
    if pdb_id not in tim_neighbours:
        continue
    keep = tim_neighbours[pdb_id]
    focus = keep[-1]
    interface = sync.get_reference_contact_sequence(focus)
    interface_list.append(interface)


# In[23]:


print(f"-350 ...{interface_list[0][-350:]}")
print()
print(f"OVERALL INTERFACE for Timeless in {pdb_ids} mapped to Q9UNS1")
print()
print("".join('1' if '1' in c else 'X' if not any(x in '01' for x in c) else '0' for c in zip(*interface_list)))

