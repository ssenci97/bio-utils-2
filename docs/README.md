# Documentation

Start with **Concepts** if you are new: it defines the words the rest of the pages use. Then follow the guides in order; each one takes a few minutes and runs on the same small example structure, [1AAY](https://www.rcsb.org/structure/1AAY).

| Page | Read it to… |
|---|---|
| [Concepts](concepts.md) | understand chain keys, label vs author IDs, the three residue numbers, moieties and units (with diagrams) |
| [1. Getting started](guides/01-getting-started.md) | install, load a structure, see its chains, find contacts, draw a first map |
| [2. Chains and numbering](guides/02-chains-and-numbering.md) | address chains safely, set reference sequences, renumber, rename and subset |
| [3. Contacts](guides/03-contacts.md) | compute contacts, contact strings, matrices, statistics and contact maps |
| [4. Interfaces](guides/04-interfaces.md) | measure buried surface, interface residues and partners |
| [5. Export](guides/05-export.md) | save tables, figures, panel grids and 3D pages |
| [6. Docking with HADDOCK3](guides/06-docking.md) | set up, run and read back a HADDOCK3 docking |
| [API index](api_index.md) | look up any method, its arguments and the configuration keys |
| [Development](development.md) | conventions, tests and checks for contributors |

Every guide is executed by the test suite (`tests/test_docs.py`), so the code you copy from it works.
