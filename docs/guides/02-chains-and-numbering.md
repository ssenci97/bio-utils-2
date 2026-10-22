# 2. Chains and numbering

**Goal:** know which chain you are talking to, read a residue's three numbers, and keep them correct when you change a reference sequence, renumber, rename or subset chains.

**You need:** [guide 1](01-getting-started.md) done. The terms used here are defined in [Concepts](../concepts.md).

## Step 1: keys, label IDs and author IDs

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")          # keyed by author IDs (the default)
protein = s["A"]                               # a ChainRecord
print(protein.chain_id, protein.auth_chain_id)
```

```text
C A
```

The protein is stored under key `"A"` (its author ID); its label ID is `"C"`. Any method that takes a chain accepts the key, the label ID or the author ID, and **keys are matched first**. So with the default keys, `"A"` is the protein, even though the DNA strand `B` has label ID `"A"`.

To key chains by label IDs instead, load with `key_by="label"`:

```python
by_label = SyncedStructure.from_rcsb("1AAY", key_by="label")
print(by_label.chain_ids)
```

```text
['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
```

With label keys, ions and water get chains of their own (`D` to `I` here), because the PDB gives them separate label IDs.

## Step 2: the three numbers of a residue

```python
mapping = s.get_residue_mapping("A")
print(mapping.head(3))
```

```text
  chain_id auth_chain_id  pdb_resnum  sequence_index  reference_index
0        C             A         103               1                3
1        C             A         104               2                4
2        C             A         105               3                5
```

The first modelled residue has author number 103, is position 1 of what was observed and position 3 of the full 90-residue sequence. `s.residue_table` is the same information for every chain at once, with residue names and moieties.

## Step 3: use your own reference sequence

The mmCIF entity sequence is the default reference. To number residues against another sequence, for example a UniProt isoform with an extra leading residue, set it; the mapping is rebuilt immediately:

```python
uniprot_like = "MS" + s.get_reference_sequence("A")[1:]   # one extra residue at the start
s.set_reference_sequence("A", uniprot_like)
print(s.get_residue_mapping("A").head(3))
```

```text
  chain_id auth_chain_id  pdb_resnum  sequence_index  reference_index
0        C             A         103               1                4
1        C             A         104               2                5
2        C             A         105               3                6
```

Every table and figure that reports reference positions now follows the new numbering. Record where the sequence came from with `s.set_refseq_id("A", "P08046-2")`.

## Step 4: renumber residues

```python
numbered = SyncedStructure.from_rcsb("1AAY")
numbered.apply_reference_numbering(["A"])     # author numbers := reference positions, in place
print(numbered.get_residue_mapping("A")["pdb_resnum"].head(3).tolist())
```

```text
[3, 4, 5]
```

`apply_reference_numbering` also drops residues that have no reference position. To renumber a chain 1..N without touching the structure, use `renumbered = s.renumber_chain("A")` (it returns a Biopython chain).

## Step 5: rename and subset

```python
renamed = SyncedStructure.from_rcsb("1AAY")
renamed.rename_chains({"A": "P"})              # key A becomes P
print(renamed.chain_ids)

report = SyncedStructure.from_rcsb("1AAY").rename_chains_by_sequence({"ZF": s.get_sequence("A")})
print(report)

pair = s.subset_chains(["A", "B"])             # a new structure with only these chains
print(pair.structure_id, pair.chain_ids)
```

```text
['B', 'C', 'P']
  new_chain_id old_chain_id old_auth_chain_id  identity_percent           status
0           ZF            C                 A             100.0          renamed
1          NaN            A                 B               NaN  unmatched_chain
2          NaN            B                 C               NaN  unmatched_chain
['B', 'A']
```

- With `key_by="auth"` the key *is* the author ID, so renaming changes the author ID; the label ID stays.
- `rename_chains_by_sequence` needs at least 90 % identity over at least half of the target sequence (`min_identity`, `min_coverage`). Chains it cannot place are listed as `unmatched_chain`.
- `subset_chains` keeps reference sequences and metadata; the original structure is not changed.

To find which chains carry a sequence without renaming anything, use `s.find_chains_by_sequence("RPYACPVESCDRRFSRSDEL")`.

## Common mistakes

- **Mixing ID systems.** Tables report `chain_id` (label) and `auth_chain_id` (author) side by side. When you take a chain from a table to pass it back to a method, use the column that matches `key_by`.
- **Comparing chains by `sequence_index`.** Copies of the same protein often miss different residues, so `sequence_index` 50 is not the same residue in two chains. Use `reference_index`.
- **Editing the Biopython structure directly.** If you change `s.structure` yourself, call `s.remap_all_chains()` afterwards, or use `add_chain`, `replace_chain` and `remove_chain`, which keep the mapping in sync.

## Next

- [3. Contacts](03-contacts.md)
- [Concepts](../concepts.md) for the definitions used here.
