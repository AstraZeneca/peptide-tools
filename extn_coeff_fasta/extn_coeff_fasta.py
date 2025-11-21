import argparse
import csv
import json
import os
import textwrap
import warnings
from typing import Any
from typing import Dict

from Bio import SeqIO

# Global definitions
known_basic_res = ["K", "R", "H"]
known_acidic_res = ["D", "E", "C", "Y", "U"]
known_res = [
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
    "𝒞",
]

aa_table = [
    ("Ala", "A", "Alanine"),
    ("Arg", "R", "Arginine"),
    ("Asn", "N", "Asparagine"),
    ("Asp", "D", "Aspartic acid"),
    ("Cys", "C", "Cysteine"),
    ("Gln", "Q", "Glutamine"),
    ("Glu", "E", "Glutamic acid"),
    ("Gly", "G", "Glycine"),
    ("His", "H", "Histidine"),
    ("Ile", "I", "Isoleucine"),
    ("Leu", "L", "Leucine"),
    ("Lys", "K", "Lysine"),
    ("Met", "M", "Methionine"),
    ("Phe", "F", "Phenylalanine"),
    ("Pro", "P", "Proline"),
    ("Pyl", "O", "Pyrrolysine"),
    ("Ser", "S", "Serine"),
    ("Sec", "U", "Selenocysteine"),
    ("Thr", "T", "Threonine"),
    ("Trp", "W", "Tryptophan"),
    ("Tyr", "Y", "Tyrosine"),
    ("Val", "V", "Valine"),
    ("Asx", "B", "Aspartic acid or Asparagine"),
    ("Glx", "Z", "Glutamic acid or Glutamine"),
    ("Xaa", "X", "Any amino acid"),
    ("Xle", "J", "Leucine or Isoleucine"),
    ("TERM", "", "Termination codon"),
]


# Turns a dictionary into a class
class Dict2Class(object):
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])


def read_fasta_file(input_file: str) -> Dict[int, Dict[str, Any]]:
    """
    Reads a FASTA file and returns a dictionary mapping unique IDs
    to molecule information.

    Parameters
    ----------
    input_file : str
        Path to the FASTA file.

    Returns
    -------
    Dict[int, Dict[str, Any]]
        Dictionary where each key is a unique molecule ID, and
        the value is a dictionary containing 'mol_name', 'mol_obj',
        and 'fasta' sequence.
    """
    _, ext = os.path.splitext(input_file)
    if ext.lower() != ".fasta":
        warnings.warn(
            'File extension is not ".fasta". Assuming FASTA format and continuing.',
            UserWarning,
        )

    mol_supply_json: Dict[int, Dict[str, Any]] = {}
    with open(input_file, "r") as handle:
        for mol_id, record in enumerate(SeqIO.parse(handle, "fasta"), start=1):
            mol_supply_json[mol_id] = {
                "mol_name": record.id,
                "mol_obj": None,
                "fasta": str(record.seq),
            }

    return mol_supply_json


def calc_extn_coeff(options=None):
    if options is None:
        options = {}

    args = Dict2Class(options)

    # Determine molecule input
    if getattr(args, "seq", ""):
        # Single fasta input
        mol_supply_json = {
            1: {"mol_name": "unknown", "mol_obj": None, "fasta": args.seq}
        }
    elif getattr(args, "inputFile", ""):
        # Input from fasta file
        mol_supply_json = read_fasta_file(args.inputFile)
    elif getattr(args, "inputJSON", ""):
        # Input from JSON string
        mol_supply_json = json.loads(args.inputJSON)
    elif getattr(args, "inputDict", None):
        # Input from dictionary
        mol_supply_json = args.inputDict
    else:
        raise ValueError(
            "Error: Provide either seq, inputFile, inputJSON, or inputDict."
        )

    # Calculate extinction coefficients
    dict_out = {}
    for mol_id, mol_data in mol_supply_json.items():
        fasta = mol_data["fasta"]
        coeffs = calc_extn_coeff_single_sequence(fasta, args)
        coeffs["mol_name"] = mol_data.get("mol_name", f"mol_{mol_id}")
        dict_out[mol_id] = coeffs
    return dict_out


def calc_extn_coeff_single_sequence(orig_sequence, args):
    # Convert d Amino Acids
    sequence = _convert_sequence_left_handed_aa(orig_sequence)

### Temporarily do it at the peptide_tools_master level
#   # Convert Cysteines
#   if args.num_disulfide_cys_bonds:
#       num_disulfide_cys_bonds = args.num_disulfide_cys_bonds
#       num_cys = sequence.count("C")
#       # We can't have more bonds than Cysteine residues
#       if num_cys / num_disulfide_cys_bonds < 2:
#           raise ValueError(
#               f"Specified {num_disulfide_cys_bonds} `num_disulfide_cys_bonds` "
#               f"but only {num_cys} Cysteine residues are present."
#           )
#       num_c_to_replace = num_disulfide_cys_bonds * 2
#       replaced_sequence = ""
#       for c in sequence:
#           if num_c_to_replace > 0:
#               if c == "C":
#                   replaced_sequence += "𝒞"
#                   num_c_to_replace -= 1
#                   continue
#           replaced_sequence += c
#       sequence = replaced_sequence

    # Validation
    for R in sequence:
        if R not in known_res:
            raise Exception(
                f"Error: Residue {R} is not known. "
                "Please use X if this is a non-canonical residue."
            )

    # Residue counts
    residue_counts = {res: sequence.count(res) for res in known_res}
    nPepBond = len(sequence) - 1

    # Coefficients dictionaries
    e205_coeffs = {
        "W": 20400,
        "F": 8600,
        "Y": 6080,
        "H": 5200,
        "M": 1830,
        "R": 1350,
        "𝒞": 1100,
        "C": 690,
        "N": 400,
        "Q": 400,
    }

    e214_coeffs = {
        "W": 29050,
        "F": 5200,
        "Y": 5375,
        "H": 5125,
        "M": 980,
        "R": 102,
        "C": 225,
        "𝒞": 225,
        "N": 136,
        "Q": 142,
        "A": 32,
        "D": 58,
        "E": 78,
        "G": 21,
        "I": 45,
        "L": 45,
        "K": 41,
        "S": 34,
        "T": 41,
        "V": 43,
    }

    # Special case: Proline for e214
    nP_rest = sequence[1:].count("P")
    nP_first = 1 if sequence[0] == "P" else 0

    # Compute sums
    e205 = (
        sum(residue_counts.get(res, 0) * coeff for res, coeff in e205_coeffs.items())
        + 2780 * nPepBond
    )
    e214 = (
        sum(residue_counts.get(res, 0) * coeff for res, coeff in e214_coeffs.items())
        + 2675 * nP_rest
        + 30 * nP_first
        + 923 * nPepBond
    )
    e280 = (
        5500 * residue_counts.get("W", 0)
        + 1490 * residue_counts.get("Y", 0)
        + 0.5 * 125 * residue_counts.get("𝒞", 0)
    )

    return {"fasta": orig_sequence, "e205": e205, "e214": e214, "e280": e280}


def _convert_sequence_left_handed_aa(sequence: str):
    return sequence.upper()


def print_stdout(dict_in):
    separator = "=" * 120
    references = [
        'Pace, Vajdos, Fee, Grimsley, "How to measure and predict the molar absorption coefficient of a protein", Protein Science 1995, 4, 2411-2423',  # noqa
        'Kuipers, Gruppen, "Prediction of molar extinction coefficients of proteins and peptides ...", J. Agric. Food Chem. 2007, 55, 5445',  # noqa
        'Anthis, Clore, "Sequence-specific determination of protein and peptide concentrations by absorbance at 215 nm", Protein Science 2013, 22, 851',  # noqa
    ]

    for molid, dict_extn_coeff in dict_in.items():
        mol_name = dict_extn_coeff["mol_name"]
        fasta_seq = dict_extn_coeff["fasta"]

        print(separator)
        print("--- Molar absorption coefficient at different wavelengths ---")
        print(f"Mol ID: {mol_name}")
        print(f"Sequence: {fasta_seq}")
        print(f"e(205 nm) = {dict_extn_coeff['e205']} (M*cm)^-1")
        print(f"e(214 nm) = {dict_extn_coeff['e214']} (M*cm)^-1")
        print(f"e(280 nm) = {dict_extn_coeff['e280']} (M*cm)^-1")
        print()
        for ref in references:
            print(ref)
        print()  # Add a blank line between entries


if __name__ == "__main__":

    # Parse options
    usage_text = textwrap.dedent(
        """
        extn_coeff_fasta.py calculates peptide extinction coefficients
        (molar absorption coefficients) based on a FASTA sequence.

        Usage:
            python extn_coeff_fasta.py -s GGKGD
    """
    )

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            "extn_coeff_fasta.py calculates peptide extinction coefficients "
            "(molar absorption coefficients) based on a FASTA sequence."
        ),
        usage=usage_text,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-s", action="store", dest="seq", help="peptide sequence", default=""
    )
    parser.add_argument(
        "-i",
        dest="inputFile",
        help="input file with molecule structure. smi or sdf",
        default="",
    )
    parser.add_argument(
        "-o",
        dest="outputFile",
        help="output file with molecule structure. csv",
        default="",
    )
    parser.add_argument(
        "--json",
        default=False,
        action="store_true",
        dest="l_json",
        help="Print output as JSON",
    )
#    parser.add_argument(
#        "--num_disulfide_cys_bonds", default=0, dest="num_disulfide_cys_bonds", type=int
#    )
    args = parser.parse_args()

    dict_extn_coeff = calc_extn_coeff(args.__dict__)

    # Output
    if args.outputFile == "":  # output plain text
        if args.l_json:
            print(json.dumps(dict_extn_coeff, indent=2))
        else:
            print_stdout(dict_extn_coeff)

    else:  # output file

        known_out_file_types = [".csv"]
        filename, out_fext = os.path.splitext(args.outputFile)
        if out_fext not in known_out_file_types:
            raise Exception("Error! Output file extention not in supported file types")

        elif out_fext == ".csv":
            with open(args.outputFile, "w") as csv_f:
                csv_w = csv.writer(csv_f)
                count = 0
                for mi in dict_extn_coeff.keys():
                    count += 1
                    if count == 1:
                        header = ["mol_name", "fasta", "e205", "e214", "e280"]
                        csv_w.writerow(header)

                    row = []
                    row += [dict_extn_coeff[mi]["mol_name"]]
                    row += [dict_extn_coeff[mi]["fasta"]]
                    row += [dict_extn_coeff[mi]["e205"]]
                    row += [dict_extn_coeff[mi]["e214"]]
                    row += [dict_extn_coeff[mi]["e280"]]
                    csv_w.writerow(row)

        # print info
        dict_file = {
            "outputFile": args.outputFile,
            "outputInfo": "Number of molecules processed:"
            + str(len(dict_extn_coeff.keys())),
        }
        print(json.dumps(dict_file))
