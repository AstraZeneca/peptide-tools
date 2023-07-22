import os
import pichemist


# FASTA definitions
KNOWN_BASIC_RESIDUES = ["K", "R", "H"]
KNOWN_ACIDIC_RESIDUES = ["D", "E", "C", "Y", "U"]
KNOWN_RESIDUES = ["G", "A", "S", "P", "V", "T", "C",
                  "I", "L", "N", "D", "Q", "K", "E",
                  "M", "H", "F", "R", "Y", "W", "X",
                  "Z", "B", "U"]

PKA_JSON_INDICES = {
     "primary": 0,
     "n-terminal": 1,
     "c-terminal": 2
}

REFERENCE_PKA_SET = "IPC2_peptide"

PKA_SETS_NAMES = ["IPC2_peptide",
                  "IPC_peptide",
                  "ProMoST",
                  "Gauci",
                  "Grimsley",
                  "Thurlkill",
                  "Lehninger",
                  "Toseland"]

ALL_KNOWN_PKA_SETS_NAMES = ["ProMoST",
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
                            "Dawson"]

# SMARTS definitions
SKIP_SMARTS_NAMES = []
PKA_LIMITS = {
    "acid_1": -5,
    "acid_2": 12,
    "base_1": 2,
    "base_2": 15,
}

# Plot definitions
PLOT_LINE_WIDTHS = {
    "w1": 4.0,
    "w2": 3.0,
    "w3": 2.0,
    "w4": 1.0
}

# Paths
SRC_DIR = pichemist.__path__[0]
ROOT_DIR = os.path.dirname(SRC_DIR)
DATA_DIR = f"{ROOT_DIR}/data"
FASTA_PKA_DATA_DIR = f"{DATA_DIR}/fasta"
FASTA_PKA_SETS_DIR = f"{FASTA_PKA_DATA_DIR}/standardised"
SMARTS_PKA_DATA_DIR = f"{DATA_DIR}/smarts"
SS_SMARTS_PKA_SET_FILEPATH = f"{SMARTS_PKA_DATA_DIR}/standardised/" \
                            "ss_smarts_set_standardised.json"
AA_SMARTS_SET_FILEPATH = f"{SMARTS_PKA_DATA_DIR}/standardised/" \
                         "aa_smarts_set_standardised.json"
TITRATION_FILE_PREFIX = "out_titration_curve"
