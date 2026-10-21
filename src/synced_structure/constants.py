"""Constants, default configuration, and column schemas for synchronized structure analysis."""

###################################### IMPORTS AND CONSTANTS

from Bio.Data.PDBData import nucleic_letters_3to1_extended
from Bio.PDB.Polypeptide import protein_letters_3to1_extended

RCSB_URL = "https://files.rcsb.org/download/{}.cif"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_RESOLUTION_MIN_DEFAULT = 0.0
RCSB_RESOLUTION_MAX_DEFAULT = 100.0
RCSB_COUNT_MAX_DEFAULT = 9999
RCSB_RETURN_TYPES = ("entry", "assembly", "polymer_entity", "non_polymer_entity", "polymer_instance", "mol_definition")
PUB_COLUMNS = [
    "pdb_id", "citation_id", "is_primary", "title", "authors", "journal", "volume",
    "page_first", "page_last", "year", "doi", "pubmed_id",
]
CONTACTS_CUTOFF_NM = 0.8
DEFAULT_CONTACT_TYPE = "inter"
CONTACT_TYPES = ("inter", "intra", "all")
# MDTraj scheme name for each contact scheme
MD_SCHEMES = {"ca": "ca", "sidechain": "closest-heavy"}
CONTACT_SCHEMES = tuple(MD_SCHEMES)
CONFIG_SECTION = "synced_structure"
DEFAULT_CONFIG = {
    "save_outdir": "data/saved",
    "aln_match_score": 1,
    "aln_mismatch_score": -1,
    "aln_open_gap_score": -2,
    "aln_extend_gap_score": -0.5,
    "interface_contact_cutoff_a": 4.5,
    "interface_score_cutoff_a": 12.0,
    "interface_block_size": 1000,
}
# Categorical palette for chains, chain pairs, and contacts across the package
PALETTE = ["#ECDE90", "#93DEB6", "#F1B584", "#76BBE2", "#DA919B",
           "#899E64", "#8E81B9", "#A66C58", "#30767A", "#7A5271"]
SUPPORTED_SUFFIXES = {".pdb", ".cif", ".mmcif"}
EXPORT_SUFFIX = ".tsv"
EXPORT_SEP = "\t"
KEY_MODES = ("auth", "label")

MOIETY_PROTEIN = "protein"
MOIETY_NUCLEIC = "nucleic"
MOIETY_LIGAND = "ligand"
MOIETY_ION = "ion"
MOIETY_WATER = "water"
MOIETIES = (MOIETY_PROTEIN, MOIETY_NUCLEIC, MOIETY_LIGAND, MOIETY_ION, MOIETY_WATER)
POLYMER_MOIETIES = (MOIETY_PROTEIN, MOIETY_NUCLEIC)
DEFAULT_CONTACT_MOIETIES = POLYMER_MOIETIES

### token dictionaries
PROTEIN_3TO1 = {k.strip().upper(): v for k, v in protein_letters_3to1_extended.items()}
NUCLEIC_TO1 = {k.strip().upper(): v for k, v in nucleic_letters_3to1_extended.items()}
PROTEIN_RESNAMES = frozenset(PROTEIN_3TO1)
NUCLEIC_RESNAMES = frozenset(NUCLEIC_TO1)
WATER_RESNAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})
ION_RESNAMES = frozenset({
    # Alkali & Alkaline Earth Metals
    "LI", "NA", "K", "RB", "CS",
    "MG", "CA", "SR", "BA",
    # Transition & Heavy Metals (commonly in metalloproteins/cofactors)
    "ZN", "FE", "CU", "MN", "NI", "CO", "CD", "HG", "PB", "AL",
    "AG", "AU", "PT", "PD", "SN", "SB", "BI", "CR", "V", "MO",
    "W", "TI", "GA", "IN", "TL", "SC", "Y", "LA", "CE", "PR",
    "ND", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
    "LU", "U", "TH",
    # Halides & Simple Anions
    "CL", "BR", "I", "F",
    # Common Physiological / Buffer Counter-Ions
    "SO4", "PO4", "NO3", "NH4", "ACT"
})
REPRESENTATIVE_ATOMS = {MOIETY_PROTEIN: ("CA",), MOIETY_NUCLEIC: ("C4'", "C1'", "P")}
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})
# phosphate atoms of a nucleotide; sugar atoms are recognized by a prime or asterisk in the name
NUCLEIC_PHOSPHATE_ATOMS = frozenset({"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"})

ANGSTROM_PER_NM = 10.0
COORD_BLOCK_ROWS = 256
INTERFACE_ENERGY_PREFACTOR = -0.5

GLOBAL_MATCH_SCORE = 2.0
GLOBAL_MISMATCH_SCORE = -1.0
GLOBAL_GAP_OPEN_SCORE = -5.0
GLOBAL_GAP_EXTEND_SCORE = -0.5
ALN_SCORE_MISMATCH_WEIGHT = 0.5
ALN_SCORE_GAP_WEIGHT = 1.5

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
MATCH_COLUMNS = ["chain_id", "auth_chain_id", "identity_percent", "aligned_residues", "query_coverage_percent"]
RENAME_COLUMNS = ["new_chain_id", "old_chain_id", "old_auth_chain_id", "identity_percent", "status"]
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
PAIR_STAT_COLUMNS = [
    "chain_a", "chain_b", "contact_type", "contact_pairs", "contacted_a", "contacted_b",
    "fraction_a", "fraction_b", "contact_density", "min_nm", "median_nm", "mean_nm",
]
PROFILE_COLUMNS = [
    "structure_id", "chain_id", "auth_chain_id", "sequence_index", "resnum", "resname",
    "moiety", "n_contacts", "nearest_nm",
]
SEPARATION_COLUMNS = ["chain_id", "separation", "min_distance_nm"]
MOIETY_PAIR_COLUMNS = ["moiety_pair", "contact_pairs", "median_nm", "mean_nm"]
COMPLETE_MIN_SEPARATION = 2
# Kyte-Doolittle hydropathy per one-letter residue code
KYTE_DOOLITTLE_HYDROPATHY = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}
# Tien et al. maximum accessible surface area in A^2 per three-letter residue name
TIEN_MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}
