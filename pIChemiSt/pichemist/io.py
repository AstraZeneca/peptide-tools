import os
import json

from enum import Enum
from rdkit import Chem
from pichemist.model import InputFormat


class FileExtension(str, Enum):
    SDF = "sdf"
    SMILES = "smi"

    @classmethod
    def repr(cls):
        return str({x: FileExtension[x].value
                    for x in FileExtension.__members__})


class FileFormatError(Exception):
    pass


def generate_input(input_format, input):
    """Produces an input dictionary compatible with the logic."""
    input_dict = dict()
    if input_format == InputFormat.SMILES_FILE:
        input_dict = read_structure_file(input)
    if input_format == InputFormat.SMILES:
        input_dict[1] = {"mol_name": input,
                         "mol_obj": Chem.MolFromSmiles(input),
                         "fasta": None}
    if input_format == InputFormat.JSON:
        input_dict = json.loads(input)
    return input_dict


def read_structure_file(input_filepath):
    """
    Reads a file containing molecule structures.
    It guesses the input type according to the extension
    of the file.

    """
    _, ext = os.path.splitext(input_filepath)
    if not ext:
        raise FileNotFoundError("Something wrong with the file "
                                f"{input_filepath}")

    # Initialize file reader
    if ext[1:] == FileExtension.SDF.value:
        suppl = Chem.SDMolSupplier(input_filepath)
    elif ext[1:] == FileExtension.SMILES.value:
        suppl = Chem.SmilesMolSupplier(input_filepath, titleLine=False)
    else:
        raise FileFormatError("Warning: Only the formats "
                              f"'{FileExtension.repr()} are accepted")

    # Populate input and assign properties
    dict_input = dict()
    uuid = 1
    for mol in suppl:
        if not mol.HasProp("_Name"):
            mol.SetProp('_Name', 'tmp_' + str(uuid))
        dict_input[uuid] = {"mol_name": mol.GetProp('_Name'),
                            "mol_obj": mol,
                            "fasta": None}
        uuid += 1
    return dict_input
