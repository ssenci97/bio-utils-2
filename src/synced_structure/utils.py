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
import math
import os
import sys
import time
from collections.abc import Iterator, Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
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

###################################### DISTANCES


def blocked_distances(xa: np.ndarray, xb: np.ndarray, block: int) -> Iterator[tuple[int, np.ndarray]]:
    """Walk the pairwise distances of two point sets a block of xa rows at a time, bounding memory.

    Args:
        xa: (n, 3) coordinates.
        xb: (m, 3) coordinates, in the same unit.
        block: Number of xa rows per block.

    Yields:
        (first xa row of the block, distances of shape (rows in block, m)).
    """
    for start in range(0, len(xa), block):
        yield start, np.linalg.norm(xa[start:start + block, None] - xb[None], axis=-1)

###################################### CPU CORES

CGROUP_ROOT = Path("/sys/fs/cgroup")
# batch schedulers export the cores granted to the job; Slurm sets SLURM_CPUS_PER_TASK only with --cpus-per-task
SCHEDULER_CORE_VARS = ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "PBS_NCPUS", "NCPUS", "LSB_DJOB_NUMPROC", "NSLOTS")


def _cgroup_quota_cores(proc_cgroup: Path = Path("/proc/self/cgroup"), root: Path = CGROUP_ROOT) -> int | None:
    """Return the CPU quota of this process's cgroups, rounded up to whole cores.

    Reads cpu.max (cgroup v2) or cpu.cfs_quota_us and cpu.cfs_period_us (cgroup v1) of the process's cgroup and of
    each parent, and keeps the tightest. Containers, Kubernetes, systemd slices and Slurm set these quotas.

    Args:
        proc_cgroup: The /proc/<pid>/cgroup file naming the process's cgroups.
        root: Mount point of the cgroup hierarchy.

    Returns:
        The quota in cores, or None without a quota or on systems without cgroups.
    """
    try:
        lines = proc_cgroup.read_text().splitlines()
    except OSError:
        return None
    quotas: list[float] = []
    for line in lines:
        _, controllers, path = line.split(":", 2)
        if controllers == "":                                   # cgroup v2: 0::/path
            base, files = root, ["cpu.max"]
        elif "cpu" in controllers.split(","):                   # cgroup v1: n:cpu,cpuacct:/path
            base = next((root / name for name in ("cpu,cpuacct", "cpu") if (root / name).is_dir()), root / "cpu")
            files = ["cpu.cfs_quota_us", "cpu.cfs_period_us"]
        else:
            continue
        folder = base / path.lstrip("/")
        # inside a cgroup namespace the process path may not exist below the mount; the mount itself is the cgroup
        for directory in (folder, *folder.parents):
            if not directory.is_relative_to(base):
                break
            try:
                # cpu.max holds 'quota period'; v1 keeps them in two files
                quota, period = " ".join((directory / name).read_text() for name in files).split()[:2]
            except (OSError, ValueError):
                continue
            if quota not in ("max", "-1") and float(period) > 0:
                quotas.append(float(quota) / float(period))
    return max(1, math.ceil(min(quotas))) if quotas else None


def core_count() -> int:
    """Return the number of CPU cores this process may use.

    os.cpu_count() counts every core of the node, also when a job or container is granted only some of them. This
    takes the smallest of the CPU affinity mask (taskset, Slurm, numactl), the cgroup CPU quota (containers,
    Kubernetes, systemd, Slurm), and the cores a batch scheduler grants (Slurm, PBS, LSF, SGE environment variables).

    Returns:
        The number of usable cores, at least 1.
    """
    limits = [len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)]
    quota = _cgroup_quota_cores()
    if quota is not None:
        limits.append(quota)
    for name in SCHEDULER_CORE_VARS:
        value = os.environ.get(name, "")
        # Slurm may write SLURM_CPUS_ON_NODE as '8(x2)'; the count comes first
        digits = value.split("(")[0].strip()
        if digits.isdigit() and int(digits) > 0:
            limits.append(int(digits))
            break
    return max(1, min(limits))

###################################### CONFIG AND LOGGING


def load_config(yaml_path: PathLike | None, defaults: Mapping[str, Any], section: str) -> dict[str, Any]:
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
