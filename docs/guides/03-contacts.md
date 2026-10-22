# 3. Contacts

**Goal:** list residue contacts, turn them into per-residue strings, matrices and statistics, and draw contact maps.

**You need:** [guide 1](01-getting-started.md). Distances in this guide are in **nm** (see [Units](../concepts.md#units)).

## Step 1: choose a scheme and a cutoff

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")

coarse = s.compute_contacts()                                  # scheme="ca", cutoff=0.8 nm
detailed = s.compute_contacts(scheme="sidechain", cutoff=0.45)
print(len(coarse), len(detailed))
```

```text
8 69
```

| Scheme | Measures | Typical cutoff |
|---|---|---|
| `"ca"` (default) | One point per residue: Cα (protein), C4′ (nucleic acids) | 0.8 nm |
| `"sidechain"` | Closest heavy atom of side chains (protein) or bases (nucleic acids) | 0.4–0.5 nm |

`contact_type` selects `"inter"` (between chains, the default), `"intra"` (within a chain) or `"all"`.

## Step 2: restrict to the chains you care about

`inter` means *any* two different chains, so in 1AAY it also counts the base pairs between the two DNA strands:

```python
print(detailed.groupby(["auth_chain_1", "auth_chain_2"]).size())

protein_dna = s.compute_contacts(scheme="sidechain", cutoff=0.45, chain_ids=["A", "B"])
print(len(protein_dna))
```

```text
auth_chain_1  auth_chain_2
B             A               24
              C               27
C             A               18
dtype: int64
24
```

## Step 3: contact strings

A contact string has one character per residue of a chain: `1` in contact, `0` not in contact, `X` no contact point.

```python
print(s.contact_string("A", scheme="sidechain", cutoff=0.45))
print(s.reference_contact_string("A", scheme="sidechain", cutoff=0.45))
```

```text
0000000000000001111001000000000000000000000101100100000000000000000000010110010000000
XX0000000000000001111001000000000000000000000101100100000000000000000000010110010000000XXX
```

The second string is laid out on the full reference sequence: the `X`s are the residues that were never modelled. Use it to compare chains or structures position by position. `s.export_contact_strings("contacts.tsv")` writes the strings of every chain to a table.

## Step 4: matrices and statistics

```python
pair = s.contact_matrix("A", "B", scheme="sidechain", cutoff=0.45)
print(pair["matrix"].shape)          # protein residues x DNA residues, distance in nm or NaN

stats = s.contact_statistics(scheme="sidechain", cutoff=0.45)
print(list(stats))
print(stats["pairs"].query("contact_type == 'inter'")[["chain_a", "chain_b", "contact_pairs"]])
```

```text
(85, 11)
['pairs', 'profile', 'separation', 'moiety_pairs']
  chain_a chain_b  contact_pairs
1       A       B             27
2       A       C             24
4       B       C             18
```

The statistics tables name chains by **label ID**: here `C` is the protein, `A` and `B` the DNA strands. `stats["profile"]` counts contacts per residue, which is the quickest way to find interface hot spots.

## Step 5: contact maps

```python
fig, ax = s.plot_contact_map("A", "B", scheme="sidechain", cutoff=0.45, crop=True)
fig.savefig("map_A_B.pdf")

fig, ax = s.plot_multi_contact_map(["A", "B", "C"], scheme="sidechain", cutoff=0.45)
fig.savefig("map_all.pdf")
```

- `plot_contact_map` draws one chain pair. `crop=True` zooms on the contacting residues. Pairs drawn in gray are within 1.5 times the cutoff (set `gray_cutoff=` to change it), which shows near misses.
- `plot_multi_contact_map` puts several chains on one block map: each pair in its own color, intra-chain contacts in gray, chain boundaries dashed.
- Axis ticks are author residue numbers.

## Common mistakes

- **Writing the cutoff in Å.** `cutoff=4.5` is 4.5 nm. Use `0.45`.
- **Forgetting that `inter` includes every chain pair.** Pass `chain_ids=` to look at one interface.
- **Reading statistics chains as keys.** `contact_statistics` reports label IDs; map them back with `{record.chain_id: key for key, record in s.chains.items()}`.

## Next

- [4. Interfaces](04-interfaces.md): buried area and interface metrics in Å.
- [5. Export](05-export.md): save maps and tables in one go.
