#!/usr/bin/env python
"""Parallel interface analysis over a directory of structures with ensemble-wide residue frequencies."""

###################################### IMPORTS AND CONSTANTS

import multiprocessing
import os
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from constants import CONFIG_SECTION, DEFAULT_CONFIG, DEFAULT_CONFIG_PATH, INTERFACE_COMPLEX_COLUMNS, INTERFACE_RESIDUE_COLUMNS, SUPPORTED_SUFFIXES
from helpers import finish_log, resolve_config

###################################### ENSEMBLE FREQUENCIES


class EnsembleMixin:
    @classmethod
    def add_interface_frequencies(cls, frame):
        """Add single_freq and pair_freq columns counting the structures containing each residue or pair.

        Args:
            frame: Residue-level interface DataFrame with a structure column.
        """
        frame = frame.copy()
        keys = (frame["chain"].astype(str) + frame["resid"].astype(str) + " " + frame["resname"].astype(str)).tolist()
        partner_lists = [cls._split_partners(p) for p in frame["partners"]]
        singles, pairs = defaultdict(set), defaultdict(set)
        for name, key, plist in zip(frame["structure"], keys, partner_lists):
            singles[key].add(name)
            for partner in plist:
                pairs["|".join(sorted((key, partner)))].add(name)
        frame["single_freq"] = [len(singles[k]) for k in keys]
        frame["pair_freq"] = [
            max((len(pairs["|".join(sorted((k, p)))]) for p in plist), default=0)
            for k, plist in zip(keys, partner_lists)
        ]
        return frame

###################################### ENSEMBLE RUNNER

    @classmethod
    def interface_ensemble(cls, input_dir, chain_a="A", chain_b="B", workers=None,
                           config_path=DEFAULT_CONFIG_PATH, outdir=None):
        """Compute interfaces for every structure in a directory in parallel.

        Args:
            input_dir: Directory containing structure files.
            chain_a: First chain identifier.
            chain_b: Second chain identifier.
            workers: Number of worker processes, or None for all CPUs.
            config_path: Path to the YAML config file.
            outdir: Optional directory for the TSV score tables.

        Returns:
            Dict with 'complex' and 'residues' DataFrames and a 'skipped' list.
        """
        start = time.perf_counter()
        input_dir = Path(input_dir).expanduser()
        files = sorted(p for p in input_dir.iterdir() if p.suffix.lower() in SUPPORTED_SUFFIXES)
        if not files:
            raise FileNotFoundError(f"No structure files found in {input_dir}")
        resolve_config(config_path, DEFAULT_CONFIG, CONFIG_SECTION)
        jobs = [(str(p), chain_a, chain_b, str(config_path), cls) for p in files]
        complexes, frames, skipped = [], [], []
        with multiprocessing.Pool(workers or os.cpu_count()) as pool:
            for name, complex_row, frame, error in tqdm(
                pool.imap_unordered(_interface_job, jobs), total=len(jobs), desc="interfaces",
            ):
                if error:
                    skipped.append(name)
                    print(f"[PYTHON-INFO] {name}: skipped, {error}")
                    continue
                complexes.append(complex_row)
                frames.append(frame)
        complex_frame = pd.DataFrame(complexes, columns=INTERFACE_COMPLEX_COLUMNS).sort_values("structure", ignore_index=True)
        residue_frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=INTERFACE_RESIDUE_COLUMNS)
        residue_frame = cls.add_interface_frequencies(residue_frame).sort_values(["structure", "chain", "resid"], ignore_index=True)
        if outdir:
            outdir = cls._prepare_dir(outdir)
            tag = input_dir.resolve().name
            complex_frame.to_csv(outdir / f"{tag}_complex_scores.tsv", sep="\t", index=False)
            residue_frame.to_csv(outdir / f"{tag}_residue_scores.tsv", sep="\t", index=False)
        finish_log(start, f"Interfaces {chain_a}-{chain_b} for {len(complexes)} structures, skipped {len(skipped)}")
        return {"complex": complex_frame, "residues": residue_frame, "skipped": skipped}

###################################### WORKERS


def _interface_job(job):
    """Compute one structure's interface in a worker process.

    Args:
        job: Tuple of (path, chain_a, chain_b, config_path, owner_class).

    Returns:
        Tuple of (name, complex_row, residue_frame, error_message_or_None).
    """
    path, chain_a, chain_b, config_path, owner = job
    name = Path(path).stem
    try:
        result = owner.from_Path(path, config_path=config_path).interface(chain_a, chain_b)
        return name, result["complex"], result["residues"], None
    except Exception as ex:
        return name, None, None, str(ex).replace("\n", " ")
