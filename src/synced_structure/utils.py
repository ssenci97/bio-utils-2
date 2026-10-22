"""Configuration loading and timestamped logging.

Log lines go to the 'synced_structure' logger. Until the application configures it, lines are printed to stdout
with the [PYTHON-INFO] prefix; clear its handlers to silence them.
"""

###################################### IMPORTS AND CONSTANTS

from __future__ import annotations

import datetime
import importlib
import importlib.util
import logging
import sys
import time
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

import yaml

from ._types import PathLike

# Windows has no resource module
resource: ModuleType | None = importlib.import_module("resource") if importlib.util.find_spec("resource") else None

LOGGER_NAME = "synced_structure"
LOG_PREFIX = "[PYTHON-INFO] "

logger = logging.getLogger(LOGGER_NAME)
if not logger.handlers:
    _handler = logging.StreamHandler(sys.stdout)
    _handler.setFormatter(logging.Formatter(LOG_PREFIX + "%(message)s"))
    logger.addHandler(_handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False

###################################### CONFIG AND LOGGING


def resolve_config(yaml_path: PathLike | None, defaults: Mapping[str, Any], section: str) -> dict[str, Any]:
    """Return a section of a YAML config with missing keys filled from defaults; never writes the file.

    Args:
        yaml_path: Path to the YAML file, or None to use defaults only.
        defaults: Default values for the section.
        section: Top-level YAML key that holds this section.

    Returns:
        The section's values over the defaults.
    """
    if yaml_path is None:
        return dict(defaults)
    path = Path(yaml_path).expanduser()
    data = (yaml.safe_load(path.read_text()) if path.is_file() else None) or {}
    return {**defaults, **(data.get(section) or {})}


def log_info(message: str) -> None:
    """Log an informational message on the package logger.

    Args:
        message: Text to log.
    """
    logger.info(message)


def _peak_memory_mb() -> float | None:
    """Return the peak resident memory of this process.

    Returns:
        Peak memory in MB, or None where the resource module is unavailable.
    """
    if resource is None:
        return None
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return usage / 1024**2 if sys.platform == "darwin" else usage / 1024


def log_completion(start_time: float, message: str) -> None:
    """Log a timestamped line with elapsed time and peak memory.

    Args:
        start_time: Value of time.perf_counter() taken at the start of the task.
        message: Description of the completed task.
    """
    elapsed = time.perf_counter() - start_time
    mem_mb = _peak_memory_mb()
    memory = "n/a" if mem_mb is None else f"{mem_mb:.1f} MB"
    timestamp = datetime.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S")
    log_info(f"{timestamp} | {message} | execution time: {elapsed:.2f}s | max memory: {memory}")
