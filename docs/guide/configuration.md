# Configuration

Every loader accepts `config_path=`, a YAML file whose `synced_structure` section overrides these defaults. Keys left out keep their defaults, and the file is never written.

| Key | Default | Meaning |
|---|---|---|
| `save_outdir` | `'data/saved'` | Default output folder of every `export_` method. |
| `aln_match_score` | `1` | Local alignment match score (alignment tables, sequence matching). |
| `aln_mismatch_score` | `-1` | Local alignment mismatch score. |
| `aln_open_gap_score` | `-2` | Local alignment gap opening score. |
| `aln_extend_gap_score` | `-0.5` | Local alignment gap extension score. |
| `interface_contact_cutoff_a` | `4.5` | Heavy-atom contact cutoff of interface metrics, Å. |
| `interface_score_cutoff_a` | `12.0` | Distance range of the interface contact score, Å. |
| `interface_block_size` | `1000` | Rows per block when computing interface atom distances; lower it to save memory. |

```yaml
synced_structure:
  save_outdir: results
  interface_contact_cutoff_a: 5.0
```

```python
from synced_structure import SyncedStructure

s = SyncedStructure.from_rcsb("1AAY")               # defaults
print(s.config["interface_contact_cutoff_a"])        # 4.5
```

<!-- skip-test -->
```python
s = SyncedStructure.from_rcsb("1AAY", config_path="config.yaml")
```

## Logging

Progress lines go to the `synced_structure` logger and are printed with a `[PYTHON-INFO]` prefix until your application configures logging. To silence them:

<!-- skip-test -->
```python
import logging

logging.getLogger("synced_structure").handlers.clear()
```
