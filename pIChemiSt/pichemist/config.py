import pichemist
from pichemist.model import ACDPKaFlag

# Shared config
ROUNDING_DIGITS = 4

# FASTA definitions
KNOWN_BASIC_RESIDUES = ["K", "R", "H"]
KNOWN_ACIDIC_RESIDUES = ["D", "E", "C", "Y", "U"]
KNOWN_RESIDUES = [
    "G",
    "A",
    "S",
    "P",
    "V",
    "T",
    "C",
    "I",
    "L",
    "N",
    "D",
    "Q",
    "K",
    "E",
    "M",
    "H",
    "F",
    "R",
    "Y",
    "W",
    "X",
    "Z",
    "B",
    "U",
]

PKA_JSON_INDICES = {"primary": 0, "n-terminal": 1, "c-terminal": 2}

REFERENCE_PKA_SET = "IPC2_peptide"

PKA_SETS_NAMES = [
    "IPC2_peptide",
    "IPC_peptide",
    "ProMoST",
    "Gauci",
    "Grimsley",
    "Thurlkill",
    "Lehninger",
    "Toseland",
]

ALL_KNOWN_PKA_SETS_NAMES = [
    "ProMoST",
    "IPC_peptide",
    "IPC2_peptide",
    "Gauci",
    "Bjellqvist",
    "Rodwell",
    "Grimsley",
    "Thurlkill",
    "EMBOSS",
    "DTASelect",
    "Solomon",
    "Sillero",
    "Lehninger",
    "Toseland",
    "Nozaki",
    "Dawson",
]

# SMARTS definitions
SKIP_SMARTS_NAMES = []
PKA_LIMITS = {
    "acid_1": -5,
    "acid_2": 12,
    "base_1": 2,
    "base_2": 15,
}
AA_TYPE_KEYS = {
    "capped": "d_capped_aa_smarts",
    "nterm": "d_nterm_free_aa_smarts",
    "cterm": "d_cterm_free_aa_smarts",
}

# ACD configuration
ACD_METHOD = ACDPKaFlag.GALAS

# Plot definitions
PLOT_LINE_WIDTHS = {"w1": 4.0, "w2": 3.0, "w3": 2.0, "w4": 1.0}

# Paths
MODULE_DIR = pichemist.__path__[0]
DATA_DIR = f"{MODULE_DIR}/data"
FASTA_PKA_DATA_DIR = f"{DATA_DIR}/fasta"
FASTA_PKA_SETS_DIR = f"{FASTA_PKA_DATA_DIR}/standardised"
SMARTS_PKA_DATA_DIR = f"{DATA_DIR}/smarts"
SS_SMARTS_PKA_SET_FILEPATH = (
    f"{SMARTS_PKA_DATA_DIR}/standardised/ss_smarts_set_standardised.json"
)
AA_SMARTS_SET_FILEPATH = (
    f"{SMARTS_PKA_DATA_DIR}/standardised/aa_smarts_set_standardised.json"
)
PH_Q_FILE_PREFIX = "out_ph_q_curve"
