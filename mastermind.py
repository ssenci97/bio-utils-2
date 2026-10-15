#!/usr/bin/env python
"""Mastermind entry point that runs the SyncedStructure example with YAML-driven parameters."""

###################################### IMPORTS AND CONSTANTS

import argparse
import time
from pathlib import Path

import numpy as np
import yaml

from constants import DEFAULT_CONTACT_TYPE, DEFAULT_CUTOFF, SEED
from helpers import finish_log
from synced_structure import SyncedStructure

SCRIPT_NAME = Path(__file__).stem
DEFAULT_YAML = "configs/param_configs.yaml"
PARAMS = {"pdb_id": "1UBQ", "cutoff": DEFAULT_CUTOFF, "contact_type": DEFAULT_CONTACT_TYPE}

###################################### ARGUMENTS AND YAML


def parse_args():
    """Parse --yaml-in and --yaml-out command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--yaml-in", nargs="?", const=DEFAULT_YAML, default=None)
    parser.add_argument("--yaml-out", nargs="?", const=DEFAULT_YAML, default=None)
    return parser.parse_args()


def load_params(path, defaults):
    """Return the script's parameters from a YAML file, merged over defaults.

    Args:
        path: YAML file path.
        defaults: Default parameter dict.
    """
    data = yaml.safe_load(Path(path).read_text()) or {}
    return {**defaults, **(data.get(SCRIPT_NAME) or {})}


def save_params(path, params):
    """Write the script's parameters under its name in a YAML file, keeping other sections.

    Args:
        path: YAML file path.
        params: Parameter dict to store.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (yaml.safe_load(path.read_text()) if path.exists() else None) or {}
    data[SCRIPT_NAME] = params
    path.write_text(yaml.safe_dump(data, sort_keys=False))

###################################### MAIN


def main():
    """Run the SyncedStructure example workflow with parameters from the command line and YAML."""
    start = time.perf_counter()
    args = parse_args()
    params = load_params(args.yaml_in, PARAMS) if args.yaml_in else dict(PARAMS)
    if args.yaml_out:
        save_params(args.yaml_out, params)
    np.random.seed(SEED)
    structure = SyncedStructure.from_rcsb_id(params["pdb_id"])
    print(structure)
    print(structure.summary())
    print(structure.alignment_stats)
    print(structure.residue_map.head())
    if structure.chain_ids:
        chain_id = structure.chain_ids[0]
        contacts = structure.compute_contacts(cutoff=params["cutoff"], contact_type=params["contact_type"])
        print(f"[PYTHON-INFO] Chain: {chain_id}")
        print(f"[PYTHON-INFO] Observed sequence: {structure.get_sequence(chain_id)}")
        print(f"[PYTHON-INFO] Contact sequence:  {structure.get_contact_sequence(chain_id, contacts=contacts)}")
        print(structure.get_contact_summary(contacts=contacts))
    structure.save_all(cutoff=params["cutoff"])
    finish_log(start, "SyncedStructure mastermind completed")


if __name__ == "__main__":
    main()
