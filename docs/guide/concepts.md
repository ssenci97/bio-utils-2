# Concepts

A deposited structure is an incomplete copy of a biological sequence. Some residues are not modelled, the depositors choose the numbering, and every chain has two IDs. `synced_structure` keeps all of these, and every table reports them side by side.

## Three numbers per residue

| Number | Column | Meaning | Example (5D39, chain A, first modelled residue) |
|---|---|---|---|
| Author number | `pdb_resnum` | The residue number written in the file | 125 |
| Sequence index | `sequence_index` | Position in the **observed** sequence, the residues present in the file | 1 |
| Reference index | `reference_index` | Position in the **reference** sequence, the full sequence of the chain's entity | 6 |

Only `reference_index` means the same residue in every chain, structure and model. It is what makes contacts from two structures, or a docking restraint and a docked model, comparable. All indices in returned tables are **1-based**.

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")
s.chain_residue_table("A").head(3)     # pdb_resnum 103, 104, 105 → reference_index 3, 4, 5
```

## Chain IDs

Every chain has two IDs in mmCIF, plus the name `synced_structure` stores it under:

| Name | Column | Meaning |
|---|---|---|
| Label ID | `chain_id` | ID assigned by the PDB (`label_asym_id`). |
| Author ID | `auth_chain_id` | ID chosen by the depositors (`auth_asym_id`); the one most viewers show. |
| Chain key | `chain_key` | The name the chain is stored under: its author ID by default (`key_by="auth"`), its label ID with `key_by="label"`, or the new name after `rename_chains`. |

The same letter can name **different** chains in the two systems. In 1AAY, label chain `A` is author chain `B`. Every method accepts a key, a label ID or an author ID, and keys are matched first. Every table that lists chains has a `chain_key` column (`chain_key_a`/`chain_key_b` for pairs), so a value read from a table can always be passed back to a method:

```python
pairs = s.contact_stat_tables()["pairs"].query("contact_type == 'inter'")
key_a, key_b = pairs[["chain_key_a", "chain_key_b"]].iloc[0]
fig, ax = s.plot_contact_map(key_a, key_b)
```

## Glossary

observed sequence
:   One-letter sequence of the polymer residues present in the file.

reference sequence
:   Full sequence of the chain's entity, read from the mmCIF file or set with `set_reference_sequence`.

mapping
:   The observed → reference correspondence of a chain, built by global alignment. It is rebuilt when a chain's observed or reference sequence changes, and kept otherwise, so a mapping set by `apply_reference_numbering` survives later renames.

moiety
:   Type of a residue: `protein`, `nucleic` or `other` (ligands, ions, water). Polymer residues are `protein` and `nucleic`.

contact point
:   One point per residue: Cα (protein), C4′ (nucleic) or the heavy-atom centroid (other).

contact scheme
:   `ca`: distance between contact points. `sidechain`: closest side-chain (protein) or base (nucleic) heavy atoms.

contact type
:   `inter` (between chains), `intra` (within a chain) or `all`.

contact string
:   One character per residue: `1` contact, `0` no contact, `X` no contact point.

interface
:   The residues of two chains that bury surface against each other, measured with **all heavy atoms** at the `interface_contact_cutoff_a` distance (Å). This is not the same as a contact, which uses the contact scheme and a cutoff in nm.

partner
:   One rigid body of a HADDOCK3 docking; it can hold several chains.
