import os
import tempfile

from peptools.chem import get_fasta_from_mol
from peptools.chem import get_fasta_from_smiles
from rdkit import Chem


class FileFormatException(Exception):
    pass


class IOException(Exception):
    pass


class RuntimeParameters:
    def __init__(self):
        self.mol_name = "none"
        self.filepath = None
        self.filepath_prefix = None
        self.input_filepath = ""  # TODO
        self.input_file_extension = None
        self.output_filename = ""  # TODO
        self.output_file_extension = None
        self.output_dir = None
        self.delete_temp_file = False
        self.generate_plots = True
        self.calc_extn_coeff = False
        self.calc_pIChemiSt = False
        self.calc_pI_fasta = False


class FileExtension:
    SDF = ".sdf"
    SMI = ".smi"
    FASTA = ".fasta"
    SMILES = ".smiles"


ACCEPTED_FILE_FORMATS = [".sdf", ".smi", ".smiles", ".fasta"]


def generate_input(input_data):
    input_data = input_data.encode("utf-8").decode("unicode_escape")
    params = RuntimeParameters()
    mol_supply_json = dict()
    input_data = _polish_input(input_data, params)

    # Validate input
    if not input_data:
        raise IOException("Input data is empty.")

    # Multiline input
    if is_input_multiline(input_data):
        input_data = multiline_input_into_filepath(input_data, params)

    # Input is a file path
    if os.path.exists(input_data):
        mol_supply_json = read_file(input_data, params)

    elif not input_data.isalpha():
        # print("Input recognised as SMILES")
        # Assume it is smiles, if contains not only letters
        smiles_elements = input_data.split()
        input_data = smiles_elements[0]
        if len(smiles_elements) > 1:
            input_data = smiles_elements[0]
            mol_name = smiles_elements[1]

        mol_unique_ID = 1
        smi = input_data
        fasta = get_fasta_from_smiles(smi)
        params.calc_extn_coeff = True
        params.calc_pI_fasta = False
        params.calc_pIChemiSt = True
        mol = Chem.MolFromSmiles(smi)
        mol_supply_json[mol_unique_ID] = {
            "mol_name": params.mol_name,
            "mol_obj": mol,
            "fasta": get_fasta_from_mol(mol),
        }

    elif input_data.isalpha():
        # Assume it is FASTA, if contains only letters
        # print("Input is FASTA")
        mol_unique_ID = 1
        fasta = input_data

        # smi = ""
        params.calc_extn_coeff = True
        params.calc_pI_fasta = True
        params.calc_pIChemiSt = False
        mol_supply_json[mol_unique_ID] = {
            "mol_name": params.mol_name,
            "mol_obj": None,
            "fasta": fasta,
        }

    else:
        raise Exception(
            "ERROR: input not recongnized: not smiles, not fasta, not a known databaase ID. Must be a bug. Contact developer. Exit."
        )
    return mol_supply_json, params


def _polish_input(input_data, params):
    input_data = input_data.strip()
    input_data = input_data.replace("ENDOFLINE", "\n")
    return input_data


def is_input_multiline(input_data):
    return "\n" in input_data


def multiline_input_into_filepath(input_data, params):
    input_list = input_data.split("\n")
    suffix = recognize_input_suffix(input_list)
    tf = tempfile.NamedTemporaryFile(
        prefix="tmp_peptide_tools_master", suffix=suffix, delete=False
    )
    input_data = tf.name
    with open(input_data, "w") as f:
        for line in input_list:
            f.write(line + "\n")
    params.delete_temp_file = True
    return input_data


def recognize_input_suffix(input_list):
    if _is_input_fasta(input_list):
        suffix = ".fasta"
    elif _is_input_sdf(input_list):
        suffix = ".sdf"
    else:
        suffix = ".smi"
    return suffix


def _is_input_fasta(input_list):
    # https://en.wikipedia.org/wiki/FASTA_format
    if input_list[0][0] == ">" or input_list[0][0] == ";":
        return True
    return False


def _is_input_sdf(input_data):
    return "$$$$" in input_data


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
    if params.input_file_extension == ".sdf":
        params.output_file_extension = ".sdf"
    params.output_filename = (
        f"{params.filepath_prefix}_OUTPUT{params.output_file_extension}"
    )

    # Runtime parameters
    params.generate_plot = False
    params.calc_extn_coeff = True
    params.calc_pIChemiSt = True
    params.calc_pI_fasta = False
    if params.input_file_extension == ".fasta":
        params.calc_pI_fasta = True
        params.calc_pIChemiSt = False

    # Read file
    if params.input_file_extension in [".sdf", ".smi", ".smiles"]:
        mol_supply_json = read_structure_file(input_data)
    elif params.input_file_extension == ".fasta":
        mol_supply_json = read_fasta_file(input_data)

    # Delete if temporary file
    if params.delete_temp_file:
        os.remove(params.input_filepath)
    return mol_supply_json


def read_structure_file(inputFile):

    filename, ext = os.path.splitext(inputFile)

    # Initialize file reader
    if ext == FileExtension.SDF:
        suppl = Chem.SDMolSupplier(inputFile)
    elif ext == FileExtension.SMI:
        suppl = Chem.SmilesMolSupplier(inputFile, titleLine=False)
    else:
        raise Exception(
            "!Warning: extension of file is not smi or sdf. Assume it is smi. Continue. "
        )
        suppl = Chem.SmilesMolSupplier(inputFile, titleLine=False)

    mol_supply_json = {}
    mol_unique_ID = 0
    for mol in suppl:
        mol_unique_ID += 1
        # print(mol_unique_ID)
        # unique index, mol title, fasta
        # fasta = get_fasta_from_smiles(smi)

        if not mol.HasProp("_Name"):
            mol.SetProp("_Name", "tmpname" + str(mol_unique_ID))

        mol_supply_json[mol_unique_ID] = {
            "mol_name": mol.GetProp("_Name"),
            "mol_obj": mol,
            "fasta": get_fasta_from_mol(mol),
        }

    return mol_supply_json


def read_fasta_file(inputFile):

    filename, ext = os.path.splitext(inputFile)

    # Initialize file reader
    if not ext == ".fasta":
        raise Exception(
            '!Warning: extension of file is not ".fasta". Assuming it is fasta formatted input. Continue. '
        )

    from Bio import SeqIO

    biosuppl = SeqIO.parse(open(inputFile), "fasta")

    mol_supply_json = {}
    mol_unique_ID = 0
    for biofasta in biosuppl:
        mol_unique_ID += 1
        # unique index, mol title, RDkit mol object, mol fasta
        mol_supply_json[mol_unique_ID] = {
            "mol_name": biofasta.id,
            "mol_obj": None,
            "fasta": str(biofasta.seq),
        }

    return mol_supply_json
