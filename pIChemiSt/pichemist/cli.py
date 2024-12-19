import argparse
import sys

from pichemist.api import pichemist_from_dict
from pichemist.config import PH_Q_FILE_PREFIX
from pichemist.io import generate_input
from pichemist.io import output_results
from pichemist.model import InputFormat
from pichemist.model import MODELS
from pichemist.model import OutputFormat
from pichemist.model import PKaMethod
from pichemist.utils import get_logger
from rdkit import RDLogger

# Configure logging
log = get_logger(__name__)
RDLogger.DisableLog("rdApp.info")

__prog__ = "pIChemiSt"
__doc__ = """Calculates the isoelectric point (pI) of a molecule by cutting
its amide bonds, retreiving the pKa values for known AAs, predicting pKas
of unknown fragments using pKaMatcher or ACDlabs, and finally calculating
the pI using the Henderson-Hasselbalch equation."""


def arg_parser(args):
    parser = argparse.ArgumentParser(prog=__prog__, description=__doc__)
    parser.add_argument(
        "-i", dest="input", help="Input filepath (or SMILES string)", default=None
    )
    parser.add_argument(
        "-if",
        dest="input_format",
        help="Format of the input",
        choices=MODELS[InputFormat],
        default=InputFormat.SMILES_FILE.value,
    )
    parser.add_argument("-o", dest="output_file", help="Output filepath", default=None)
    parser.add_argument(
        "-of",
        dest="output_format",
        help="Format of the output",
        choices=MODELS[OutputFormat],
        default=OutputFormat.CONSOLE.value,
    )
    parser.add_argument(
        "--plot_ph_q_curve",
        default=False,
        action="store_true",
        dest="plot_ph_q_curve",
        help="Generate an image of the " "pH/Q curve into a file",
    )
    parser.add_argument(
        "-pp",
        dest="ph_q_curve_file_prefix",
        help="pH/Q plot file prefix",
        default=PH_Q_FILE_PREFIX,
    )
    parser.add_argument(
        "--print_fragment_pkas",
        default=False,
        action="store_true",
        dest="print_fragment_pkas",
        help="Print the fragments with corresponding " "pKas used in pI calcution.",
    )
    parser.add_argument(
        "--generate_fragment_base64_images",
        default=False,
        action="store_true",
        dest="generate_fragment_base64_images",
        help="Generate 2D depiction of the frgament smiles in base64 format.",
    )
    parser.add_argument(
        "--ionizable_nterm",
        default=True,
        action="store_true",
        dest="ionizable_nterm",
        help="For FASTA input only. Set if the N-terminus is not capped (free amine)",
    )
    parser.add_argument(
        "--ionizable_cterm",
        default=True,
        action="store_true",
        dest="ionizable_cterm",
        help="For FASTA input only. Set if the C-terminus is not capped (free carboxylic acid)",
    )
    parser.add_argument(
        "--method",
        choices=MODELS[PKaMethod],
        default=PKaMethod.PKA_MATCHER.value,
        help="Method for the prediction of the pKa of unknown fragments",
    )
    if not args:
        args = ["-h"]
    return parser.parse_args(args)


def run_pichemist(args):
    """High-level wrapper of pIChemiSt."""
    input_dict = generate_input(args.input_format, args.input)

    if (
        args.input_format == "smiles_stdin"
        or args.input_format == "smiles_file"
        or args.input_format == "sdf"
    ):
        input_format = "structure"
    elif args.input_format == "fasta_stdin" or args.input_format == "fasta_file":
        input_format = "fasta"
    else:
        input_format = "unknown"
        raise ValueError("input_format is not known in run_pichemist")

    output_dict = pichemist_from_dict(
        input_dict,
        args.method,
        ph_q_curve_file_prefix=args.ph_q_curve_file_prefix,
        plot_ph_q_curve=args.plot_ph_q_curve,
        print_fragments=args.print_fragment_pkas,
        input_type=input_format,
        ionizable_nterm=args.ionizable_nterm,
        ionizable_cterm=args.ionizable_cterm,
        generate_fragment_base64_images=args.generate_fragment_base64_images,
    )
    output_results(
        input_dict,
        output_dict,
        args.output_file,
        args.output_format,
        args.method,
        args.print_fragment_pkas,
    )


def main():
    """Entry point for the CLI."""
    args = arg_parser(sys.argv[1:])
    run_pichemist(args)


if __name__ == "__main__":
    main()
