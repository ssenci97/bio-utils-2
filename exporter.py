#!/usr/bin/env python
"""Writing PDB, FASTA, alignment, contact, and interface tables with a YAML manifest, using standard and author chain IDs."""

###################################### IMPORTS AND CONSTANTS

import copy
import datetime
import io
import itertools
import time
from pathlib import Path

import pandas as pd
import yaml
from Bio.PDB import PDBIO
from tqdm import tqdm

from constants import DEFAULT_CUTOFF, NUMBERINGS
from helpers import finish_log

###################################### PATHS


class ExportMixin:
    @staticmethod
    def _prepare_path(path):
        """Expand a file path, create its parent directory, and return it as a Path."""
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def _prepare_dir(path):
        """Expand a directory path, create it, and return it as a Path."""
        path = Path(path).expanduser()
        path.mkdir(parents=True, exist_ok=True)
        return path

###################################### STRUCTURE AND SEQUENCE EXPORT

    def to_pdb_string(self, chain_id=None, numbering="original"):
        """Return the structure or one chain as PDB text with author chain IDs and the chosen numbering.

        Args:
            chain_id: Chain to export, or None for all chains.
            numbering: 'original', 'sequential', or 'reference'.
        """
        if numbering not in NUMBERINGS:
            raise ValueError(f"numbering must be one of {NUMBERINGS}.")
        selected = {self._require_chain(chain_id)} if chain_id is not None else set(self.chains)
        structure = copy.deepcopy(self.structure)
        model = next(structure.get_models())
        for chain in list(model.get_chains()):
            if str(chain.id) not in selected:
                model.detach_child(chain.id)
        for chain in list(model.get_chains()):
            data = self.chains[str(chain.id)]
            if numbering == "sequential":
                self._renumber(chain)
            elif numbering == "reference" and data["sequence_index_to_ref_index"]:
                self._keep_mapped(chain, data["sequence_index_to_ref_index"])
            chain.id = data["auth_chain_id"]
        writer = PDBIO()
        writer.set_structure(structure)
        output = io.StringIO()
        writer.save(output)
        return output.getvalue()

    def save_pdb(self, output_path, chain_id=None, numbering="original"):
        """Write the structure or one chain as a PDB file and return its path.

        Args:
            output_path: Destination file path.
            chain_id: Chain to export, or None for all chains.
            numbering: 'original', 'sequential', or 'reference'.
        """
        output_path = self._prepare_path(output_path)
        output_path.write_text(self.to_pdb_string(chain_id=chain_id, numbering=numbering))
        return output_path

    def to_fasta(self, output_path, reference=False):
        """Write observed or reference sequences as FASTA with standard and author chain IDs in headers.

        Args:
            output_path: Destination file path.
            reference: Whether to write reference sequences, skipping chains without one.
        """
        lines = []
        for data in self.chains.values():
            if reference and not data["reference_sequence"]:
                continue
            sequence = data["reference_sequence"] if reference else data["sequence"]
            suffix = " reference" if reference else ""
            lines.extend([f">{self.structure_id}_{data['chain_id']} auth_{data['auth_chain_id']}{suffix}", sequence])
        output_path = self._prepare_path(output_path)
        output_path.write_text("\n".join(lines) + "\n")
        return output_path

    def to_biopython(self):
        """Return a deep copy of the underlying Biopython Structure."""
        return copy.deepcopy(self.structure)

    def export_alignments(self, outdir=None):
        """Write alignment stats and positions as TSV files and return the output directory.

        Args:
            outdir: Output directory, or None to use the configured alignment directory.
        """
        outdir = self._prepare_dir(outdir or self.config["aln_outdir"])
        self.alignment_stats.to_csv(outdir / f"{self.structure_id}_stats.tsv", sep="\t", index=False)
        self.alignment_positions.to_csv(outdir / f"{self.structure_id}_positions.tsv", sep="\t", index=False)
        return outdir

###################################### SYSTEMATIC SAVING

    def save_all(self, outdir=None, numbering=NUMBERINGS, chain_pairs=None, cutoff=DEFAULT_CUTOFF):
        """Save all structure, mapping, alignment, contact, sequence, and interface tables with a manifest.

        Args:
            outdir: Output directory, or None to use the configured save directory.
            numbering: Numbering modes to write as PDB files.
            chain_pairs: Chain pairs for interface tables, or None for all pairs.
            cutoff: Contact cutoff in nm.

        Returns:
            List of written file paths, including the manifest.
        """
        start = time.perf_counter()
        outdir = self._prepare_dir(outdir or self.config["save_outdir"])
        sid = self.structure_id
        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all")
        tasks = [(f"{sid}_{mode}.pdb", lambda p, m=mode: self.save_pdb(p, numbering=m)) for mode in numbering]
        tasks += [
            (f"{sid}_chains.tsv", lambda p: self.summary().to_csv(p, sep="\t", index=False)),
            (f"{sid}_residue_map.tsv", lambda p: self.residue_map.to_csv(p, sep="\t", index=False)),
            (f"{sid}_alignment_stats.tsv", lambda p: self.alignment_stats.to_csv(p, sep="\t", index=False)),
            (f"{sid}_alignment_positions.tsv", lambda p: self.alignment_positions.to_csv(p, sep="\t", index=False)),
            (f"{sid}_contacts.tsv", lambda p: contacts.to_csv(p, sep="\t", index=False)),
            (f"{sid}_contact_summary.tsv", lambda p: self.get_contact_summary(contacts=contacts).to_csv(p, sep="\t", index=False)),
            (f"{sid}_contact_sequences.csv", lambda p: self.export_contact_sequences(p, cutoff=cutoff, contact_type="all")),
            (f"{sid}_contact_sequences_reference.csv", lambda p: self.export_contact_sequences(p, cutoff=cutoff, contact_type="all", reference=True)),
            (f"{sid}_sequences.fasta", lambda p: self.to_fasta(p)),
            (f"{sid}_reference_sequences.fasta", lambda p: self.to_fasta(p, reference=True)),
            (f"{sid}_moiety_summary.tsv", lambda p: self.moiety_summary.to_csv(p, sep="\t", index=False)),
        ]
        for chain_a, chain_b in (chain_pairs or itertools.combinations(self.chain_ids, 2)):
            chain_a, chain_b = self._require_chain(chain_a), self._require_chain(chain_b)
            pair = f"{self.chains[chain_a]['chain_id']}{self.chains[chain_b]['chain_id']}"
            try:
                result = self.interface(chain_a, chain_b)
            except ValueError as ex:
                print(f"[PYTHON-INFO] {sid}: interface {pair} skipped, {ex}")
                continue
            tasks.append((f"{sid}_{pair}_interface_complex.tsv",
                          lambda p, r=result: pd.DataFrame([r["complex"]]).to_csv(p, sep="\t", index=False)))
            tasks.append((f"{sid}_{pair}_interface_residues.tsv",
                          lambda p, r=result: r["residues"].to_csv(p, sep="\t", index=False)))
        written = []
        for name, writer in tqdm(tasks, desc="saving", unit="file"):
            path = outdir / name
            writer(path)
            written.append(str(path))
        manifest_path = outdir / f"{sid}_manifest.yaml"
        manifest = {
            "structure_id": sid,
            "file_format": self.file_format,
            "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "config": self.config,
            "files": written,
        }
        manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
        written.append(str(manifest_path))
        finish_log(start, f"Saved {len(written)} files for {sid} to {outdir}")
        return written

    def export_alignments_and_contacts(self, outdir=None, cutoff=DEFAULT_CUTOFF):
        """Alias for save_all with original numbering only, kept for compatibility."""
        return self.save_all(outdir=outdir, numbering=("original",), cutoff=cutoff)