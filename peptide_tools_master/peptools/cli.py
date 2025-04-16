#!/usr/bin/python
import argparse
import sys

from peptools.io import generate_input
from peptools.io import generate_output
from peptools.io import generate_parameter_set
from peptools.utils import str2bool
from peptools.wrapper import run_peptide_master


__prog__ = "Peptide Tools Master"
__doc__ = """TODO"""


def arg_parser(args):
    parser = argparse.ArgumentParser(prog=__prog__, description=__doc__)
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Input filepath (or SMILES string or FASTA).",
        required=True,
    )

    parser.add_argument(
        "--print_fragment_pkas",
        action="store_true",
        dest="print_fragment_pkas",
        help="Print the fragments with corresponding pKas used in pI calcution.",
        default=False,
    )

    parser.add_argument(
        "--generate_fragment_images",
        default=False,
        action="store_true",
        dest="generate_fragment_images",
        help="Generate 2D depiction of the frgament smiles in base64 format.",
    )

    parser.add_argument(
        "--print_pka_set",
        action="store_true",
        dest="print_pka_set",
        help="Print out stored pka sets.",
        default=False,
    )

    ### keys for fasta input
    parser.add_argument(
        "--ionizable_nterm",
        type=str2bool,
        default=True,
        dest="ionizable_nterm",
        help="Applies to FASTA input only. "
        "If set to 'false' the N-terminus is capped. "
        "If set to 'true' the N-terminus is free amine. ",
    )
    parser.add_argument(
        "--ionizable_cterm",
        type=str2bool,
        default=True,
        dest="ionizable_cterm",
        help="Applies to FASTA input only. "
        "If set to 'false' the C-terminus is capped. "
        "If set to 'true' the C-terminus is free amine. ",
    )
    if not args:
        args = ["-h"]
    return parser.parse_args()


def main():
    """Entry point for the CLI."""
    args = arg_parser(sys.argv[1:])
    mol_supply_json, io_params = generate_input(args.input)
    params = generate_parameter_set(args, io_params)
    dict_out = run_peptide_master(mol_supply_json, params)
    generate_output(mol_supply_json, dict_out, params)


if __name__ == "__main__":
    main()
