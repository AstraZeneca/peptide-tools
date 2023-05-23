import os
import pichemist

# Definitions
KNOWN_BASIC_RESIDUES = ["K", "R", "H"]
KNOWN_ACIDIC_RESIDUES = ["D", "E", "C", "Y", "U"]
KNOWN_RESIDUES = ["G", "A", "S", "P", "V", "T", "C",
                  "I", "L", "N", "D", "Q", "K", "E",
                  "M", "H", "F", "R", "Y", "W", "X",
                  "Z", "B", "U"]

PKA_JSON_TYPE_MATCHING = {
     "acidic": "acidic",
     "basic": "basic",
     "terminus_ionizable": "terminus_ionizable"
}

PKA_JSON_INDICES = {
     "primary": 0,
     "n-terminal": 1,
     "c-terminal": 2
}

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
