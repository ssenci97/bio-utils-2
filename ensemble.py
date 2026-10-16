#!/usr/bin/env python
"""Parallel all-pairs interface analysis over a directory of structures with ensemble-wide residue frequencies, runnable from the command line with YAML parameters."""

###################################### IMPORTS AND CONSTANTS

import argparse
import multiprocessing
import os
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml
from tqdm import tqdm

from constants import CONFIG_SECTION, DEFAULT_CONFIG, DEFAULT_CONFIG_PATH, INTERFACE_COMPLEX_COLUMNS, INTERFACE_RESIDUE_COLUMNS, SUPPORTED_SUFFIXES
from helpers import finish_log, resolve_config
from loaders import KEY_MODES

SCRIPT_NAME = "ensemble"

###################################### ENSEMBLE FREQUENCIES


class EnsembleMixin:
    @classmethod
    def add_interface_frequencies(cls, frame):
        """Add single_freq and pair_freq columns counting the structures containing each residue or residue pair.

        Args:
            frame: Residue-level interface DataFrame with structure, chain, resid, resname, partners columns.
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
    def interface_ensemble(cls, input_dir, workers=None, config_path=DEFAULT_CONFIG_PATH,
                           key_by="auth", outdir=None):
        """Compute all chain-pair interfaces for every structure in a directory in parallel.

        Args:
            input_dir: Directory containing structure files.
            workers: Number of worker processes, or None for all CPUs.
            config_path: Path to the YAML config file.
            key_by: 'auth' (default) keys chains by author IDs, 'label' by label IDs.
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
        jobs = [(str(p), str(config_path), key_by, cls) for p in files]
        complexes, residues, skipped = [], [], []
        with multiprocessing.Pool(workers or os.cpu_count()) as pool:
            for name, complex_frame, residue_frame, error in tqdm(
                pool.imap_unordered(_interface_job, jobs), total=len(jobs), desc="structures", unit="file",
            ):
                if error:
                    skipped.append(name)
                    print(f"[PYTHON-INFO] {name}: skipped, {error}")
                    continue
                complexes.append(complex_frame)
                residues.append(residue_frame)
        complex_all = (pd.concat(complexes, ignore_index=True) if complexes
                       else pd.DataFrame(columns=INTERFACE_COMPLEX_COLUMNS))
        complex_all = complex_all.sort_values(["structure", "chain_a", "chain_b"], ignore_index=True)
        residue_all = (pd.concat(residues, ignore_index=True) if residues
                       else pd.DataFrame(columns=INTERFACE_RESIDUE_COLUMNS))
        residue_all = cls.add_interface_frequencies(residue_all).sort_values(
            ["structure", "pair", "chain", "resid"], ignore_index=True,
        )
        if outdir:
            outdir = cls._prepare_dir(outdir)
            tag = input_dir.resolve().name
            complex_all.to_csv(outdir / f"{tag}_complex_scores.tsv", sep="\t", index=False)
            residue_all.to_csv(outdir / f"{tag}_residue_scores.tsv", sep="\t", index=False)
        finish_log(start, f"All-pair interfaces for {len(complexes)} structures, skipped {len(skipped)}")
        return {"complex": complex_all, "residues": residue_all, "skipped": skipped}

###################################### WORKERS


def _interface_job(job):
    """Compute all chain-pair interfaces of one structure in a worker process.

    Args:
        job: Tuple of (path, config_path, key_by, owner_class).

    Returns:
        Tuple of (name, complex_frame, residue_frame, error_message_or_None).
    """
    path, config_path, key_by, owner = job
    name = Path(path).stem
    try:
        result = owner.from_Path(path, config_path=config_path, key_by=key_by).interface_all()
        return name, result["complex"], result["residues"], None
    except Exception as ex:
        return name, None, None, str(ex).replace("\n", " ")

###################################### COMMAND LINE


def main():
    """Parse arguments, optionally read them from YAML, optionally write them to YAML, and run the ensemble."""
    parser = argparse.ArgumentParser(description="All-pair interface analysis over a directory of structures.")
    parser.add_argument("--input-dir", default=None)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--config-path", default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--key-by", default="auth", choices=KEY_MODES)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--yaml-in", default=None)
    parser.add_argument("--yaml-out", default=None)
    args = parser.parse_args()
    if args.yaml_in:
        saved = (yaml.safe_load(Path(args.yaml_in).read_text()) or {}).get(SCRIPT_NAME) or {}
        parser.set_defaults(**saved)
        args = parser.parse_args()
    if args.yaml_out:
        path = Path(args.yaml_out)
        data = (yaml.safe_load(path.read_text()) if path.exists() else None) or {}
        data[SCRIPT_NAME] = {k: v for k, v in vars(args).items() if k not in ("yaml_in", "yaml_out")}
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(data, sort_keys=False))
    if not args.input_dir:
        parser.error("--input-dir is required, either as a flag or in the --yaml-in file.")
    from synced_structure import SyncedStructure
    SyncedStructure.interface_ensemble(
        args.input_dir, workers=args.workers, config_path=args.config_path,
        key_by=args.key_by, outdir=args.outdir,
    )


if __name__ == "__main__":
    main()