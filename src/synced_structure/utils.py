"""Configuration loading and timestamped logging.

Log lines go to the 'synced_structure' logger. Until the application configures it, lines are printed to stdout
with the [PYTHON-INFO] prefix; clear its handlers to silence them.
"""

###################################### IMPORTS AND CONSTANTS

import datetime
import logging
import sys
import time
import warnings
from pathlib import Path

import yaml

try:
    import resource
except ImportError:  # Windows has no resource module
    resource = None

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


def deprecated_method(old_name, new_name):
    """Return a method that warns that old_name is deprecated and forwards the call to new_name.

    Args:
        old_name: Name of the deprecated method.
        new_name: Name of the method that replaces it, on the same object.
    """
    def method(self, *args, **kwargs):
        warnings.warn(f"{old_name}() is deprecated; use {new_name}().", DeprecationWarning, stacklevel=2)
        return getattr(self, new_name)(*args, **kwargs)

    method.__name__ = old_name
    method.__doc__ = f"Deprecated alias of {new_name}()."
    return method


def log_info(message):
    """Log an informational message on the package logger."""
    logger.info(message)


def _peak_memory_mb():
    """Return the peak resident memory of this process in MB, or None where the resource module is unavailable."""
    if resource is None:
        return None
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return usage / 1024**2 if sys.platform == "darwin" else usage / 1024


def log_completion(start_time, message):
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
