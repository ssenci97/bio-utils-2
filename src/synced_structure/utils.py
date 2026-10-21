"""Configuration loading and timestamped logging."""

###################################### IMPORTS AND CONSTANTS

import datetime
import resource
import sys
import time
from pathlib import Path

import yaml

###################################### CONFIG AND LOGGING


def resolve_config(yaml_path, defaults, section):
    """Return a section of a YAML config with missing keys filled from defaults; never writes the file.

    Args:
        yaml_path: Path to the YAML file, or None to use defaults only.
        defaults: Dict of default values for the section.
        section: Top-level YAML key that holds this section.
    """
    if yaml_path is None:
        return dict(defaults)
    path = Path(yaml_path).expanduser()
    data = (yaml.safe_load(path.read_text()) if path.is_file() else None) or {}
    return {**defaults, **(data.get(section) or {})}


def log_info(message):
    """Print a [PYTHON-INFO] log line."""
    print(f"[PYTHON-INFO] {message}")


def finish_log(start_time, message):
    """Log a timestamped line with elapsed time and peak memory.

    Args:
        start_time: Value of time.perf_counter() taken at the start of the task.
        message: Description of the completed task.
    """
    elapsed = time.perf_counter() - start_time
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    mem_mb = usage / 1024**2 if sys.platform == "darwin" else usage / 1024
    timestamp = datetime.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S")
    log_info(f"{timestamp} | {message} | execution time: {elapsed:.2f}s | max memory: {mem_mb:.1f} MB")
