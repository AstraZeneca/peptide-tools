import os

from peptools.chem import get_fasta_from_mol
from peptools.io.fasta import _is_input_fasta
from peptools.io.fasta import configure_fasta_input
from peptools.io.fasta import read_fasta_file
from peptools.io.file import FileFormatException
from peptools.io.input import ACCEPTED_FILE_FORMATS
from peptools.io.input import InputFileExtension
from peptools.io.multi import is_input_multiline
from peptools.io.multi import multiline_input_to_filepath
from peptools.io.structure import _is_input_smi
from peptools.io.structure import configure_smi_input
from peptools.io.structure import read_structure_file
from rdkit import Chem


class IOException(Exception):
    pass


class IOParameters:
    def __init__(self):
        self.mol_name = "none"
        self.filepath = None
        self.filepath_prefix = None
        self.input_filepath = None
        self.input_file_extension = None
        self.output_filename = None
        self.output_file_extension = None
        self.output_dir = None
        self.delete_temp_file = False


class RuntimeParameters:
    def __init__(self):
        self.generate_plots = True
        self.print_fragment_pkas = False
        self.calc_extn_coeff = False
        self.calc_pIChemiSt = False
        self.calc_pI_fasta = False


class ChemicalParameters:
    def __init__(
        self,
        ionized_Cterm,
        ionized_Nterm,
        NPhosphateGroups,
        NAlkylLysGroups,
        NDiAlkylLysGroups,
    ):
        self.ionized_Cterm = ionized_Cterm
        self.ionized_Nterm = ionized_Nterm
        self.NPhosphateGroups = NPhosphateGroups
        self.NAlkylLysGroups = NAlkylLysGroups
        self.NDiAlkylLysGroups = NDiAlkylLysGroups


def generate_input(input_data):
    input_data = input_data.encode("utf-8").decode("unicode_escape")
    params = IOParameters()
    mol_supply_json = dict()
    input_data = _polish_input(input_data, params)

    # Validate input
    if not input_data:
        raise IOException("Input data is empty.")

    # Multiline input
    if is_input_multiline(input_data):
        input_data = multiline_input_to_filepath(input_data, params)

    # Input is a file path
    if os.path.exists(input_data):
        mol_supply_json = read_file(input_data, params)

    # Input is FASTA
    elif _is_input_fasta(input_data):
        mol_supply_json = configure_fasta_input(input_data, params)

    # Input is SMILES
    elif _is_input_smi(input_data):
        mol_supply_json = configure_smi_input(input_data, params)
    else:
        raise FileFormatException()
    return mol_supply_json, params


def _polish_input(input_data, params):
    input_data = input_data.strip()
    input_data = input_data.replace("ENDOFLINE", "\n")
    return input_data


def read_file(input_data, params):
    params.input_filepath = input_data
    params.workdir = os.path.dirname(params.input_filepath)
    params.filename = os.path.splitext(os.path.basename(params.input_filepath))[0]

    # Validation
    params.filepath_prefix, params.input_file_extension = os.path.splitext(
        params.input_filepath
    )
    if params.input_file_extension not in ACCEPTED_FILE_FORMATS:
        raise FileFormatException(
            "Extension not supported: " + params.input_file_extension
        )

    # Configure output file
    params.output_file_extension = ".csv"
    if params.input_file_extension == InputFileExtension.SDF:
        params.output_file_extension = ".sdf"
    params.output_filename = (
        f"{params.filepath_prefix}_OUTPUT{params.output_file_extension}"
    )

    # Read file
    if params.input_file_extension in [InputFileExtension.SDF, InputFileExtension.SMI]:
        mol_supply_json = read_structure_file(input_data)
    elif params.input_file_extension == InputFileExtension.FASTA:
        mol_supply_json = read_fasta_file(input_data)
    else:
        raise FileFormatException()

    # Delete if temporary file
    if params.delete_temp_file:
        os.remove(params.input_filepath)
    return mol_supply_json


def configure_runtime_parameters(args, io_params):
    params = RuntimeParameters()
    params.generate_plot = False
    params.print_fragment_pkas = bool(args.print_fragment_pkas)
    if io_params.input_file_extension in [
        None,
        InputFileExtension.SMI,
        InputFileExtension.SDF,
    ]:  # sic - None is assumed to be SMI from STDIN
        params.calc_extn_coeff = True
        params.calc_pIChemiSt = True
        params.calc_pI_fasta = False
    elif io_params.input_file_extension == InputFileExtension.FASTA:
        params.calc_extn_coeff = True
        params.calc_pI_fasta = True
        params.calc_pIChemiSt = False
    return params


def configure_chemical_parameters(args):
    return ChemicalParameters(
        args.ionized_Cterm,
        args.ionized_Nterm,
        args.NPhosphateGroups,
        args.NAlkylLysGroups,
        args.NDiAlkylLysGroups,
    )
