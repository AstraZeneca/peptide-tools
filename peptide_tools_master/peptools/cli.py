#!/usr/bin/python
import argparse

from peptools.io import configure_chemical_parameters
from peptools.io import configure_runtime_parameters
from peptools.io import generate_input
from peptools.io import generate_output
from peptools.io import generate_parameter_set
from peptools.wrapper import run_peptide_master


__prog__ = "Peptide Tools Master"
__doc__ = """TODO"""


def arg_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Input filepath (or SMILES string or FASTA)",
        required=True,
    )

    ### pIChemiSt.py keys
    parser.add_argument(
        "--print_fragment_pkas",
        action="store_true",
        dest="print_fragment_pkas",
        help="Print the fragments with corresponding pKas used in pI calcution.",
        default=False,
    )
    parser.add_argument(
        "--print_pka_set",
        action="store_true",
        dest="print_pka_set",
        help="Print out stored pka sets.",
        default=False,
    )

    ### pI_fasta.py keys
    parser.add_argument(
        "--ionized_Cterm",
        dest="ionized_Cterm",
        action="store_true",
        help="is C-terminus ionized [COO-]?",
        default=True,
    )
    parser.add_argument(
        "--ionized_Nterm",
        dest="ionized_Nterm",
        action="store_true",
        help="is N-terminus ionized [N+]?",
        default=True,
    )
    parser.add_argument(
        "-p",
        action="store",
        dest="NPhosphateGroups",
        help="Number of phosphorilated residues. Phosphorilated residues must be denoted as X in the sequence. default = 0",
        default=0,
        type=int,
    )
    parser.add_argument(
        "-l",
        action="store",
        dest="NAlkylLysGroups",
        help="Number of monoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0",
        default=0,
        type=int,
    )
    parser.add_argument(
        "-d",
        action="store",
        dest="NDiAlkylLysGroups",
        help="Number of dinoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0",
        default=0,
        type=int,
    )

    args = parser.parse_args()
    return args


# def main(input_data):


if __name__ == "__main__":
    args = arg_parser()
    mol_supply_json, io_params = generate_input(args.input)
    params = generate_parameter_set(args, io_params)
    dict_out_peptide_tools_master = run_peptide_master(mol_supply_json, params)
    generate_output(mol_supply_json, dict_out_peptide_tools_master, params)
