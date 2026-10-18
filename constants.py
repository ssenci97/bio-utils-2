#!/usr/bin/env python
"""Constants, default configuration, and column schemas for synchronized structure analysis."""

###################################### IMPORTS AND CONSTANTS

from Bio.PDB.Polypeptide import protein_letters_3to1_extended

SEED = 42
RCSB_URL = "https://files.rcsb.org/download/{}.cif"
DEFAULT_CUTOFF = 0.8
DEFAULT_CONTACT_TYPE = "inter"
DEFAULT_CONFIG_PATH = "configs/param_configs.yaml"
CONFIG_SECTION = "synced_structure"
DEFAULT_CONFIG = {
    "chain_seqs_outdir": "data/chains_seqs",
    "chain_refseqs_outdir": "data/chains_refseqs",
    "aln_outdir": "data/aln",
    "save_outdir": "data/saved",
    "aln_match_score": 1,
    "aln_mismatch_score": -1,
    "aln_open_gap_score": -2,
    "aln_extend_gap_score": -0.5,
    "interface_contact_cutoff_a": 4.5,
    "interface_score_cutoff_a": 12.0,
    "interface_block_size": 1000,
}
SUPPORTED_SUFFIXES = {".pdb", ".cif", ".mmcif"}
NUMBERINGS = ("original", "sequential", "reference")
EXPORT_SUFFIX = ".csv"

MOIETY_PROTEIN = "protein"
MOIETY_NUCLEIC = "nucleic"
MOIETY_LIGAND = "ligand"
MOIETY_ION = "ion"
MOIETY_WATER = "water"
MOIETIES = (MOIETY_PROTEIN, MOIETY_NUCLEIC, MOIETY_LIGAND, MOIETY_ION, MOIETY_WATER)
DEFAULT_CONTACT_MOIETIES = (MOIETY_PROTEIN, MOIETY_NUCLEIC)
NUCLEIC_RESNAMES = frozenset({"A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DU", "DI"})
WATER_RESNAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})
REPRESENTATIVE_ATOMS = {MOIETY_PROTEIN: ("CA",), MOIETY_NUCLEIC: ("C4'", "C1'", "P")}
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

GLOBAL_MATCH_SCORE = 2.0
GLOBAL_MISMATCH_SCORE = -1.0
GLOBAL_GAP_OPEN_SCORE = -5.0
GLOBAL_GAP_EXTEND_SCORE = -0.5
CONTACT_COLUMNS = [
    "chain_1", "auth_chain_1", "resnum_1", "resname_1", "sequence_index_1",
    "chain_2", "auth_chain_2", "resnum_2", "resname_2", "sequence_index_2",
    "min_distance_nm", "contact_type", "moiety_1", "moiety_2",
]
CONTACT_SEQUENCE_COLUMNS = [
    "structure_id", "chain_id", "auth_chain_id", "sequence", "contacts",
    "contact_type", "mdtraj_scheme", "cutoff_nm",
]
CHAIN_COLUMNS = ["structure_id", "chain_id", "auth_chain_id", "entity_id", "description", "residues"]
RESIDUE_COLUMNS = [
    "structure_id", "chain_id", "auth_chain_id", "pdb_resnum", "resname",
    "sequence_index", "reference_index", "moiety",
]
MOIETY_SUMMARY_COLUMNS = ["structure_id", "chain_id", "auth_chain_id", "moiety", "residues", "atoms"]

MAPPING_COLUMNS = ["chain_id", "auth_chain_id", "pdb_resnum", "sequence_index", "reference_index"]
INTERFACE_COMPLEX_COLUMNS = [
    "structure", "chain_a", "auth_chain_a", "chain_b", "auth_chain_b",
    "len_res_a", "len_res_b", "len_atoms_a", "len_atoms_b", "dASA_A2",
    "initial_sasa_A2", "frac_dASA", "n_contacts", "iface_res_a", "iface_res_b",
    "min_dist_A", "contact_score",
]
INTERFACE_RESIDUE_COLUMNS = [
    "structure", "pair", "chain", "auth_chain", "resid", "resname", "dSASA_A2", "isolated_sasa_A2", "rel_dSASA",
    "n_partners", "min_dist_A", "contact_energy", "energy_per_partner", "partners",
]
HYDRO = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}
TO1 = {k.upper(): v for k, v in protein_letters_3to1_extended.items()}
TIEN = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}