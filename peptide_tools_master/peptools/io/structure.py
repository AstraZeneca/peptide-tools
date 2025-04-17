import os

from peptools.io.input import InputFactory
from peptools.io.input import InputFileExtension
from rdkit import Chem


def _is_input_smi(input_data):
    if Chem.MolFromSmiles(input_data, sanitize=False):
        return True
    return False


def _is_input_sdf(input_data):
    return "$$$$" in input_data


def configure_smi_input(smiles_str, params):
    # Extract name if provided
    smiles_elements = smiles_str.split()
    smiles_str = smiles_elements[0]
    if len(smiles_elements) > 1:
        smiles_str = smiles_elements[0]
        params.mol_name = smiles_elements[1]

    params.input_file_extension = InputFileExtension.SMI
    mol = Chem.MolFromSmiles(smiles_str)
    return {1: InputFactory.new(mol, params.mol_name)}


def read_structure_file(input_filepath):
    RDKIT_NAME_PROP = "_Name"
    filename, ext = os.path.splitext(input_filepath)

    # Initialize file reader
    if ext == InputFileExtension.SDF:
        suppl = Chem.SDMolSupplier(input_filepath)
    elif ext == InputFileExtension.SMI:
        suppl = Chem.SmilesMolSupplier(input_filepath, titleLine=False)

    # Populate input
    mol_supply_json = dict()
    mol_id = 0
    for mol in suppl:
        mol_id += 1
        if not mol.HasProp(RDKIT_NAME_PROP):
            mol.SetProp(RDKIT_NAME_PROP, "tmp_peptools" + str(mol_id))
        mol_supply_json[mol_id] = InputFactory.new(mol, mol.GetProp(RDKIT_NAME_PROP))
    return mol_supply_json
