#!/usr/bin/env python
"""Writing structure, alignment, contact, and interface tables per object, grouped by entity type, with a YAML manifest and cross-object concatenation."""

###################################### IMPORTS AND CONSTANTS

import copy
import datetime
import io
import itertools
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio.PDB import PDBIO
from tqdm import tqdm

from constants import DEFAULT_CUTOFF, NUMBERINGS
from helpers import finish_log

ENTITY_TABLES = ("chain_summary", "residue_summary", "contacts", "pairs", "pair_summary")
ALL_PAIRS = "all"
MERGED_PAIR_KEYS = ("chain_a", "chain_b")
RESIDUE_PAIR_KEYS = ("a_chain", "a_resnum", "b_chain", "b_resnum")
RESIDUE_PAIR_IDENTITY = ("a_auth_chain", "a_resname", "b_auth_chain", "b_resname")
RESIDUE_PAIR_COLUMNS = [
    "structure_id", "a_chain", "a_auth_chain", "a_resnum", "a_resname",
    "b_chain", "b_auth_chain", "b_resnum", "b_resname",
    "ca_min_distance_nm", "iface_min_dist_A", "iface_pair_energy",
    "in_ca_contact", "in_interface",
]
SIDE_1 = {"chain": "chain_1", "auth_chain": "auth_chain_1", "resnum": "resnum_1", "resname": "resname_1"}
SIDE_2 = {"chain": "chain_2", "auth_chain": "auth_chain_2", "resnum": "resnum_2", "resname": "resname_2"}

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

    @staticmethod
    def _normalize_numbering(numbering):
        """Return numbering modes as a tuple, accepting one mode string or an iterable of modes.

        Raises:
            ValueError: If any mode is not one of NUMBERINGS.
        """
        modes = (numbering,) if isinstance(numbering, str) else tuple(numbering)
        invalid = [mode for mode in modes if mode not in NUMBERINGS]
        if invalid:
            raise ValueError(f"numbering must be one of {NUMBERINGS}, got {invalid}.")
        return modes

    def _normalize_chain_pairs(self, chain_pairs):
        """Return a list of (chain_a, chain_b) tuples from None, 'all', or an iterable of pairs.

        Args:
            chain_pairs: None or 'all' for every chain pair, or an iterable of two-chain tuples.

        Raises:
            ValueError: If an entry is not a pair of chain identifiers.
        """
        if chain_pairs is None or (isinstance(chain_pairs, str) and chain_pairs.lower() == ALL_PAIRS):
            return list(itertools.combinations(self.chain_ids, 2))
        if isinstance(chain_pairs, str):
            raise ValueError(f"chain_pairs must be None, '{ALL_PAIRS}', or an iterable of (chain_a, chain_b) tuples.")
        pairs = []
        for entry in chain_pairs:
            if isinstance(entry, str) or len(entry) != 2:
                raise ValueError(f"Each chain pair must be a (chain_a, chain_b) tuple, got {entry!r}.")
            pairs.append((entry[0], entry[1]))
        return pairs

    def _tag(self, frame):
        """Return a copy of a table with a structure_id column first, if it has none."""
        if "structure_id" in frame.columns:
            return frame
        tagged = frame.copy()
        tagged.insert(0, "structure_id", self.structure_id)
        return tagged

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

###################################### ENTITY TABLES: CHAINS

    def _chain_summary_table(self, contacts, interface_complex):
        """Return one row per chain with metadata, inter-chain CA contact count, buried area, and partner count.

        Buried area is the sum of dASA over every interface the chain takes part in, so a residue
        buried against two partners contributes to both pairs.

        Args:
            contacts: Contacts DataFrame from compute_contacts.
            interface_complex: Complex-level interface DataFrame from interface_all.
        """
        base = self._tag(self.summary())
        labels = base["chain_id"]
        inter = contacts[contacts["contact_type"] == "inter"]
        ca_counts = pd.concat([inter["chain_1"], inter["chain_2"]]).value_counts()
        sides = pd.concat([
            interface_complex[["chain_a", "dASA_A2"]].set_axis(["chain", "dASA_A2"], axis=1),
            interface_complex[["chain_b", "dASA_A2"]].set_axis(["chain", "dASA_A2"], axis=1),
        ])
        buried = sides.groupby("chain")["dASA_A2"].sum()
        links = pd.concat([
            interface_complex[["chain_a", "chain_b"]].set_axis(["chain", "partner"], axis=1),
            interface_complex[["chain_b", "chain_a"]].set_axis(["chain", "partner"], axis=1),
        ])
        partner_counts = links.groupby("chain")["partner"].nunique()
        base["ca_inter_contacts"] = labels.map(ca_counts).fillna(0).astype(int)
        base["iface_dASA_sum_A2"] = labels.map(buried).fillna(0.0).round(2)
        base["iface_partner_chains"] = labels.map(partner_counts).fillna(0).astype(int)
        return base

###################################### ENTITY TABLES: RESIDUES

    def _residue_summary_table(self, contacts, interface_residues):
        """Return one row per observed protein residue with mapping, CA contact, and interface data.

        Interface sums are taken over every interface row for the residue, so a residue
        buried against two partners reports both contributions.

        Args:
            contacts: Contacts DataFrame from compute_contacts.
            interface_residues: Residue-level interface DataFrame from interface_all.
        """
        base = self.residue_map
        ends = pd.concat([
            contacts[["chain_1", "resnum_1"]].set_axis(["chain_id", "pdb_resnum"], axis=1),
            contacts[["chain_2", "resnum_2"]].set_axis(["chain_id", "pdb_resnum"], axis=1),
        ])
        ca = ends.groupby(["chain_id", "pdb_resnum"]).size().rename("ca_n_contacts").reset_index()
        iface = (
            interface_residues.groupby(["chain", "resid"])
            .agg(
                iface_n_pairs=("pair", "size"),
                iface_n_partner_contacts=("n_partners", "sum"),
                iface_dSASA_sum_A2=("dSASA_A2", "sum"),
                iface_contact_energy_sum=("contact_energy", "sum"),
                iface_min_dist_A=("min_dist_A", "min"),
            )
            .reset_index()
            .rename(columns={"chain": "chain_id", "resid": "pdb_resnum"})
        )
        table = base.merge(ca, on=["chain_id", "pdb_resnum"], how="left")
        table = table.merge(iface, on=["chain_id", "pdb_resnum"], how="left")
        counts = ["ca_n_contacts", "iface_n_pairs", "iface_n_partner_contacts"]
        table[counts] = table[counts].fillna(0).astype(int)
        table["iface_dSASA_sum_A2"] = table["iface_dSASA_sum_A2"].round(2)
        table["iface_contact_energy_sum"] = table["iface_contact_energy_sum"].round(4)
        return table

###################################### ENTITY TABLES: CONTACTS

    def _contacts_table(self, contacts):
        """Return all residue-pair contacts within the cutoff, intra- and inter-chain.

        Args:
            contacts: Contacts DataFrame from compute_contacts.
        """
        return self._tag(contacts)

###################################### ENTITY TABLES: PAIRS

    @staticmethod
    def _order_residue_pairs(frame, side_1, side_2, keep=()):
        """Return residue pairs with a_ and b_ fields, ordered by (chain, resnum) so that side a sorts first.

        Args:
            frame: DataFrame with one residue on each side, named by side_1 and side_2.
            side_1: Column names for the first residue, keyed by chain, auth_chain, resnum, resname.
            side_2: Column names for the second residue, keyed the same way.
            keep: Extra columns copied through unchanged.
        """
        c1 = frame[side_1["chain"]].astype(str).to_numpy()
        c2 = frame[side_2["chain"]].astype(str).to_numpy()
        r1 = frame[side_1["resnum"]].astype(int).to_numpy()
        r2 = frame[side_2["resnum"]].astype(int).to_numpy()
        swap = (c1 > c2) | ((c1 == c2) & (r1 > r2))
        out = {}
        for field in ("chain", "auth_chain", "resnum", "resname"):
            v1 = frame[side_1[field]].to_numpy()
            v2 = frame[side_2[field]].to_numpy()
            out[f"a_{field}"] = np.where(swap, v2, v1)
            out[f"b_{field}"] = np.where(swap, v1, v2)
        for column in keep:
            out[column] = frame[column].to_numpy()
        return pd.DataFrame(out, index=frame.index)

    def _residue_pair_interface(self, chain_pairs):
        """Return one row per heavy-atom residue pair across chains within the interface contact cutoff.

        Distances are minimum heavy-atom distances in A. Pair energy uses the same hydropathy
        expression as the residue-level table, applied to this single pair.

        Args:
            chain_pairs: List of (chain_a, chain_b) tuples.
        """
        cutoff = self.config["interface_contact_cutoff_a"]
        score_cutoff = self.config["interface_score_cutoff_a"]
        block = self.config["interface_block_size"]
        rows = []
        for chain_a, chain_b in chain_pairs:
            chain_a, chain_b = self._require_chain(chain_a), self._require_chain(chain_b)
            ga, gb = self._chain_geometry(chain_a), self._chain_geometry(chain_b)
            if not len(ga["xyz"]) or not len(gb["xyz"]):
                print(f"[PYTHON-INFO] {self.structure_id}: residue pairs {chain_a}-{chain_b} skipped, no heavy atoms")
                continue
            dmin = self._residue_min_distances(ga["xyz"], ga["owner"], gb["xyz"], gb["owner"],
                                               len(ga["residues"]), len(gb["residues"]), block)
            auth_a = self.chains[chain_a]["auth_chain_id"]
            auth_b = self.chains[chain_b]["auth_chain_id"]
            for i, j in zip(*np.nonzero(dmin <= cutoff)):
                distance = float(dmin[i, j])
                energy = (-0.5 * (ga["hydro"][i] + gb["hydro"][j]) * (1 - distance / score_cutoff)
                          if distance <= score_cutoff else 0.0)
                res_a, res_b = ga["residues"][i], gb["residues"][j]
                rows.append({
                    "chain_1": ga["chain_id"], "auth_chain_1": auth_a,
                    "resnum_1": res_a.id[1], "resname_1": res_a.get_resname(),
                    "chain_2": gb["chain_id"], "auth_chain_2": auth_b,
                    "resnum_2": res_b.id[1], "resname_2": res_b.get_resname(),
                    "iface_min_dist_A": round(distance, 3),
                    "iface_pair_energy": round(float(energy), 4),
                })
        return pd.DataFrame(rows, columns=[*SIDE_1.values(), *SIDE_2.values(),
                                           "iface_min_dist_A", "iface_pair_energy"])

    def _pairs_table(self, contacts, chain_pairs):
        """Return one row per residue pair joining inter-chain CA contacts with heavy-atom interface pairs.

        Only inter-chain contacts are merged, because the interface is defined between chains.
        Distances keep their units: ca_min_distance_nm in nm, iface_min_dist_A in A.
        The in_ca_contact and in_interface flags show which source each pair came from.

        Args:
            contacts: Contacts DataFrame from compute_contacts.
            chain_pairs: List of (chain_a, chain_b) tuples.
        """
        inter = contacts[contacts["contact_type"] == "inter"]
        ca_side = self._order_residue_pairs(inter, SIDE_1, SIDE_2, keep=("min_distance_nm",))
        ca_side = ca_side.rename(columns={"min_distance_nm": "ca_min_distance_nm"})

        interface = self._residue_pair_interface(chain_pairs)
        iface_side = self._order_residue_pairs(interface, SIDE_1, SIDE_2,
                                               keep=("iface_min_dist_A", "iface_pair_energy"))

        merged = ca_side.merge(iface_side, on=list(RESIDUE_PAIR_KEYS), how="outer",
                               indicator=True, suffixes=("_ca", "_if"))
        for column in RESIDUE_PAIR_IDENTITY:
            merged[column] = merged.pop(f"{column}_ca").combine_first(merged.pop(f"{column}_if"))
        merged["in_ca_contact"] = merged["_merge"] != "right_only"
        merged["in_interface"] = merged["_merge"] != "left_only"
        merged = merged.drop(columns="_merge")
        merged.insert(0, "structure_id", self.structure_id)
        return merged[RESIDUE_PAIR_COLUMNS].sort_values(list(RESIDUE_PAIR_KEYS), ignore_index=True)

###################################### ENTITY TABLES: PAIR SUMMARY

    @staticmethod
    def _ordered_pairs(frame, first, second):
        """Return a copy with chain_a and chain_b columns, ordered so that chain_a sorts before or equal to chain_b.

        Args:
            frame: DataFrame holding two chain-ID columns.
            first: Name of the first chain-ID column.
            second: Name of the second chain-ID column.
        """
        a = frame[first].astype(str)
        b = frame[second].astype(str)
        ordered = frame.copy()
        ordered["chain_a"] = a.where(a <= b, b)
        ordered["chain_b"] = b.where(a <= b, a)
        return ordered

    def _pair_summary_table(self, contacts, interface_complex):
        """Return one row per chain pair with inter-chain contact counts and complex-level interface metrics.

        Contact distances stay in nm (ca_ prefix). Interface distances stay in A (iface_ prefix).
        Pairs present in only one source are kept, with zero contact counts where a source has none.

        Args:
            contacts: Contacts DataFrame from compute_contacts.
            interface_complex: Complex-level interface DataFrame from interface_all.
        """
        auth_of = {data["chain_id"]: data["auth_chain_id"] for data in self.chains.values()}

        inter = contacts[contacts["contact_type"] == "inter"]
        contact_side = self._ordered_pairs(inter, "chain_1", "chain_2")
        contact_side = (
            contact_side.groupby(list(MERGED_PAIR_KEYS), as_index=False)
            .agg(ca_n_contacts=("min_distance_nm", "size"), ca_min_distance_nm=("min_distance_nm", "min"))
        )

        interface_side = self._ordered_pairs(interface_complex, "chain_a", "chain_b")
        interface_side = interface_side.drop(columns=["structure", "auth_chain_a", "auth_chain_b"])
        interface_side = interface_side.rename(columns={
            column: f"iface_{column}" for column in interface_side.columns if column not in MERGED_PAIR_KEYS
        })

        merged = contact_side.merge(interface_side, on=list(MERGED_PAIR_KEYS), how="outer")
        merged[["ca_n_contacts", "iface_n_contacts"]] = merged[["ca_n_contacts", "iface_n_contacts"]].fillna(0).astype(int)
        merged.insert(0, "structure_id", self.structure_id)
        merged.insert(2, "auth_chain_a", merged["chain_a"].map(auth_of))
        merged.insert(4, "auth_chain_b", merged["chain_b"].map(auth_of))
        return merged.sort_values(list(MERGED_PAIR_KEYS), ignore_index=True)

###################################### SYSTEMATIC SAVING

    def save_all(self, outdir=None, numbering=NUMBERINGS, chain_pairs=None, cutoff=DEFAULT_CUTOFF):
        """Save one structure into <outdir>/<structure_id>/, grouped by kind: structure, tables, and auxiliary.

        The folder layout is:
            <outdir>/<sid>/structure/   PDB and FASTA files
            <outdir>/<sid>/tables/      chain_summary, residue_summary, contacts, pairs, pair_summary
            <outdir>/<sid>/auxiliary/   alignment and contact-string files
            <outdir>/<sid>/<sid>_manifest.yaml

        Args:
            outdir: Parent directory holding one folder per structure, or None to use save_outdir.
            numbering: A numbering mode ('original', 'sequential', 'reference') or an iterable of modes.
            chain_pairs: None or 'all' for every chain pair, or an iterable of (chain_a, chain_b) tuples.
            cutoff: Contact cutoff in nm.

        Returns:
            List of written file paths, including the manifest.

        Raises:
            ValueError: If numbering or chain_pairs is invalid. Raised before any files are written.
        """
        numbering = self._normalize_numbering(numbering)
        pairs = self._normalize_chain_pairs(chain_pairs)
        start = time.perf_counter()
        sid = self.structure_id
        root = self._prepare_dir(Path(outdir or self.config["save_outdir"]) / sid)
        structure_dir = self._prepare_dir(root / "structure")
        tables_dir = self._prepare_dir(root / "tables")
        aux_dir = self._prepare_dir(root / "auxiliary")

        contacts = self.compute_contacts(cutoff=cutoff, contact_type="all")
        interface = self.interface_all(pairs)
        tables = {
            "chain_summary": self._chain_summary_table(contacts, interface["complex"]),
            "residue_summary": self._residue_summary_table(contacts, interface["residues"]),
            "contacts": self._contacts_table(contacts),
            "pairs": self._pairs_table(contacts, pairs),
            "pair_summary": self._pair_summary_table(contacts, interface["complex"]),
        }

        structure_tasks = [(structure_dir / f"{sid}_{mode}.pdb", lambda p, m=mode: self.save_pdb(p, numbering=m))
                           for mode in numbering]
        structure_tasks += [
            (structure_dir / f"{sid}_sequences.fasta", lambda p: self.to_fasta(p)),
            (structure_dir / f"{sid}_reference_sequences.fasta", lambda p: self.to_fasta(p, reference=True)),
        ]
        table_tasks = [(tables_dir / f"{name}.tsv", lambda p, df=frame: df.to_csv(p, sep="\t", index=False))
                       for name, frame in tables.items()]
        aux_tasks = [
            (aux_dir / "alignment_stats.tsv", lambda p: self.alignment_stats.to_csv(p, sep="\t", index=False)),
            (aux_dir / "alignment_positions.tsv", lambda p: self.alignment_positions.to_csv(p, sep="\t", index=False)),
            (aux_dir / "contact_sequences.csv",
             lambda p: self.export_contact_sequences(p, cutoff=cutoff, contact_type="all")),
            (aux_dir / "contact_sequences_reference.csv",
             lambda p: self.export_contact_sequences(p, cutoff=cutoff, contact_type="all", reference=True)),
        ]

        written = []
        for path, writer in tqdm(structure_tasks + table_tasks + aux_tasks, desc="saving", unit="file"):
            writer(path)
            written.append(str(path))
        manifest_path = root / f"{sid}_manifest.yaml"
        manifest = {
            "structure_id": sid,
            "file_format": self.file_format,
            "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "config": self.config,
            "chain_pairs": [list(pair) for pair in pairs],
            "tables": list(ENTITY_TABLES),
            "files": written,
        }
        manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
        written.append(str(manifest_path))
        finish_log(start, f"Saved {len(written)} files for {sid} to {root}")
        return written

    def export_alignments_and_contacts(self, outdir=None, cutoff=DEFAULT_CUTOFF):
        """Alias for save_all with original numbering only, kept for compatibility."""
        return self.save_all(outdir=outdir, numbering=("original",), cutoff=cutoff)

###################################### CONCATENATION

    @classmethod
    def concatenate_reports(cls, root_dir, outdir=None):
        """Concatenate the entity tables of every object folder into one file per entity type.

        Each object folder under root_dir must contain a tables/ subfolder, as written by save_all
        with outdir=root_dir. Output files go to root_dir/all/ (or outdir), one row per entity
        across all objects.

        Args:
            root_dir: Parent directory holding one folder per structure.
            outdir: Output directory, or None to write into root_dir/all.

        Returns:
            Dict of {entity_table: output_path}.
        """
        start = time.perf_counter()
        root = Path(root_dir).expanduser()
        outdir = cls._prepare_dir(outdir or root / "all")
        frames = defaultdict(list)
        objects = 0
        for obj in sorted(p for p in root.iterdir() if (p / "tables").is_dir()):
            objects += 1
            for entity in ENTITY_TABLES:
                path = obj / "tables" / f"{entity}.tsv"
                if path.exists():
                    frames[entity].append(pd.read_csv(path, sep="\t"))
        written = {}
        for entity in ENTITY_TABLES:
            if not frames[entity]:
                continue
            combined = pd.concat(frames[entity], ignore_index=True)
            output_path = outdir / f"{entity}.tsv"
            combined.to_csv(output_path, sep="\t", index=False)
            written[entity] = output_path
        finish_log(start, f"Concatenated {len(written)} entity tables from {objects} objects in {root}")
        return written