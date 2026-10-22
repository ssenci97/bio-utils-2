# 1. Getting started

**Goal:** install the package, load a structure, read its chains, find its contacts and draw a contact map. About 5 minutes.

**You need:** Python 3.10 or newer and an internet connection (to download the structure).

## Step 1: install

```bash
pip install "git+https://github.com/ssenci97/bio-utils-2.git"
```

If `pip` fails while building `mdtraj`, install it from conda-forge first (`conda install -c conda-forge mdtraj`) and repeat the command.

## Step 2: load a structure

We use [1AAY](https://www.rcsb.org/structure/1AAY): the Zif268 zinc-finger protein bound to a DNA duplex. It is small, so everything runs in seconds.

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")
print(s.chain_ids)
```

```text
['B', 'C', 'A']
```

These are the **chain keys**: the author chain IDs, in file order. A local file works the same way: `SyncedStructure.from_path("my_structure.cif")` (or `.pdb`).

## Step 3: look at the chains

```python
summary = s.chain_summary()
print(summary[["chain_id", "auth_chain_id", "moiety", "description", "sequence_length", "reference_length"]])
```

```text
  chain_id auth_chain_id   moiety                                      description  sequence_length  reference_length
0        A             B  nucleic  DNA (5'-D(*AP*GP*CP*GP*TP*GP*GP*GP*CP*GP*T)-3')               11                11
1        B             C  nucleic  DNA (5'-D(*TP*AP*CP*GP*CP*CP*CP*AP*CP*GP*C)-3')               11                11
2        C             A  protein             PROTEIN (ZIF268 ZINC FINGER PEPTIDE)               85                90
```

Read this table before anything else:

- The protein has **author ID `A`** but **label ID `C`**. Methods take the key (`"A"` here), and tables show both IDs.
- The protein has 85 residues in the file (`sequence_length`) out of 90 in its full sequence (`reference_length`). The other 5 were not modelled.

See [Concepts](../concepts.md#one-chain-two-ids) for why a chain has two IDs.

## Step 4: find contacts

```python
contacts = s.compute_contacts()          # residue pairs from different chains closer than 0.8 nm
print(len(contacts))
print(contacts[["auth_chain_1", "resnum_1", "resname_1", "auth_chain_2", "resnum_2", "resname_2", "min_distance_nm"]].head(3))
```

```text
8
  auth_chain_1  resnum_1 resname_1 auth_chain_2  resnum_2 resname_2  min_distance_nm
0            B         3        DC            A       153       HIS           0.7058
1            B         3        DC            A       156       THR           0.6842
2            B         4        DG            A       153       HIS           0.7950
```

Contact distances are in **nm**. The default scheme measures one point per residue (Cα for protein, C4′ for DNA); [guide 3](03-contacts.md) shows the side-chain scheme, which finds more protein–DNA contacts.

## Step 5: draw a contact map

```python
fig, ax = s.plot_contact_map("A", "B", crop=True)   # protein A against DNA strand B
fig.savefig("contacts_A_B.pdf")
```

`crop=True` zooms on the residues in contact; without it, the map shows the whole protein. In Jupyter the figure is shown automatically.

## Step 6 (Jupyter): look at it in 3D

```python
view = s.show_structure_3d(style="cartoon")
view   # in a notebook, this displays an interactive 3D viewer
```

The viewer needs a browser with internet access, because it loads 3Dmol.js from the web.

## Optional: a configuration file

Every loader accepts `config_path=` pointing to a YAML file. Only the keys you set change, and the file is never written:

```yaml
synced_structure:
  save_outdir: data/saved          # default folder for exports
  interface_contact_cutoff_a: 4.5  # Å, heavy-atom contact cutoff for interface metrics
  interface_score_cutoff_a: 12.0   # Å, distance range of the interface contact score
```

## Common mistakes

- **Using a table's `chain_id` as a key.** `chain_id` is the *label* ID. With the default `key_by="auth"`, pass `auth_chain_id` values (or the keys in `s.chain_ids`).
- **Thinking in Å for contacts.** `cutoff=4.5` means 4.5 nm here, which is huge; write `cutoff=0.45`.
- **Comparing residues by `sequence_index`.** Two copies of a protein can miss different residues; compare `reference_index` instead ([guide 2](02-chains-and-numbering.md)).

## Next

- [2. Chains and numbering](02-chains-and-numbering.md): IDs, reference sequences, renaming.
- [3. Contacts](03-contacts.md): contact tables, strings, matrices and maps in detail.
