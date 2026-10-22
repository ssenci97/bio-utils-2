# 4. Interfaces

**Goal:** measure how two chains pack against each other: buried surface, contacting atoms, interface residues and their partners.

**You need:** [guide 1](01-getting-started.md). Interface metrics use heavy atoms and are in **Å** and **Å²**, unlike contacts (nm). See [Units](../concepts.md#units).

## Step 1: one number at a time

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")

print(s.interface_buried_area("A", "B"))         # Å², surface buried when protein A meets DNA strand B
print(s.interface_buried_fraction("A", "B"))     # buried area / surface of the two isolated chains
print(s.interface_atom_contact_count("A", "B"))  # heavy-atom pairs within 4.5 Å
print(s.interface_min_distance_a("A", "B"))      # closest heavy atoms, Å
```

```text
1933.63
0.2016
324
2.548
```

`interface_initial_sasa` and `interface_contact_score` complete the set. The 4.5 Å contact cutoff and the 12 Å score range come from the configuration (`interface_contact_cutoff_a`, `interface_score_cutoff_a`; see [guide 1](01-getting-started.md#optional-a-configuration-file)).

## Step 2: the whole picture for one pair

```python
complex_row = s.interface_complex("A", "B")       # dict with every metric above, plus residue counts
residues = s.interface_residues("A", "B")         # one row per interface residue
print(complex_row["iface_res_a"], complex_row["iface_res_b"])
print(residues.sort_values("dSASA_A2", ascending=False)[["auth_chain", "resid", "resname", "dSASA_A2", "n_partners"]].head(3))
```

```text
25 11
   auth_chain  resid resname  dSASA_A2  n_partners
31          B      7      DG    132.75           6
29          B      5      DT    123.79           6
27          B      3      DC    111.57           4
```

25 protein residues and all 11 nucleotides of strand B are at the interface. `residues["partners"]` lists, for each residue, the residues it touches on the other chain (as `label ID + number + name`, for example `C146 ARG`).

## Step 3: every partner of a chain

```python
table = s.chain_partner_table("A")                # one row per chain in contact with A
print(table[["auth_chain_a", "auth_chain_b", "dASA_A2", "n_contacts"]])

pairs = s.chain_partner_structures("A")           # one two-chain SyncedStructure per partner
print(list(pairs))
```

```text
  auth_chain_a auth_chain_b  dASA_A2  n_contacts
0            A            B  1933.63         324
1            A            C  1045.77          70
['1AAY_A_B', '1AAY_A_C']
```

`chain_partner_table("A", level="residues")` gives the residue rows of all those interfaces. For every chain pair of the structure at once, use `s.interface_metrics_all()`, which returns `{"complex": table, "residues": table}`.

## Common mistakes

- **Comparing interface distances with contact distances.** `interface_min_distance_a` is in Å and uses every heavy atom (2.548 Å here); `compute_contacts` reports nm and measures the scheme's points (the closest side-chain contact of this pair is 0.2676 nm).
- **Expecting ligand residues in the tables.** Interface residues and partners are the polymer residues (protein and nucleic acids) of the two chains; ligands, ions and water are not listed.
- **Reading partner labels as keys.** `C146 ARG` uses the label ID `C`, which is the protein with author ID `A`.

## Next

- [5. Export](05-export.md): write these tables to TSV.
- [6. Docking](06-docking.md): dock chains with HADDOCK3.
