# Conventions

## Names say what they return

A name starts with a verb only when the verb fixes the return type or a side effect. Otherwise the name is a noun whose last word is the type. Names stay as short as that allows.

=== "Verbs"

    | Verb | Returns |
    |---|---|
    | `set_`, `add_`, `remove_`, `replace_`, `rename_`, `renumber_`, `remap_`, `apply_`, `clear_`, `validate_` | the object itself, so calls chain |
    | `subset_`, `from_`, `load_` | a new `SyncedStructure` (or a dict of them) |
    | `prepare_` | a `HaddockDocking` set up, nothing run |
    | `find_` | a DataFrame of matches |
    | `fetch_` | a DataFrame from a web service |
    | `search_` | a list of IDs |
    | `plot_`, `build_` | a matplotlib figure: `(Figure, Axes)` / `Figure` |
    | `show_` | IPython HTML |
    | `export_<name>` | the path(s) written for `<name>` |
    | `run_` | runs a haddock3 subprocess |
    | `get_` | a stored value, paired with `set_` |

=== "Noun endings"

    | Ending | Type |
    |---|---|
    | `_table` | `pandas.DataFrame`; always a method, called with `()` |
    | `_tables` | `dict[str, DataFrame]` |
    | `_matrix` | dict holding a matrix and its labels |
    | `_metrics`, `_map` | `dict` |
    | `_ids` | `list[str]` |
    | `_paths` | file or folder paths |
    | `_string`, `_sequence`, `_text`, `_toml` | `str` |
    | `_structures` | `dict[str, SyncedStructure]` |
    | `_trajectory` | `mdtraj.Trajectory` |

Renames replace the old name outright, without deprecated aliases; the [changelog](../changelog.md) lists every rename.

## Tables first

Analyses return DataFrames, or dicts of them. Figures are drawn from those tables, and `export_` methods write them. Columns are shared across tables:

- `structure_id` in every table that can mix structures.
- `chain_key`, `chain_id`, `auth_chain_id` for one chain; `chain_key_a`, `chain_a`, `auth_chain_a` (and `_b`) for a pair.
- `sequence_index`, `reference_index` (1-based) and `pdb_resnum` for residues.

## Units

- **nm** for contacts: `contact_table`, contact matrices, contact maps and every contact `cutoff=`. Columns end in `_nm`.
- **Å** and **Å²** for interface metrics and docking. Columns and settings end in `_A`, `_A2` or `_a`.

## Arguments

- Chain arguments accept a key, label ID or author ID. Methods that need chains have no default chain.
- `moieties` defaults to the polymer moieties, `("protein", "nucleic")`, everywhere.
- `outdir` names the output folder of every `export_` method; `None` uses the `save_outdir` setting.
