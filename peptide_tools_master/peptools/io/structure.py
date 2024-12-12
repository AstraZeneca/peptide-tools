from peptools.chem import get_fasta_from_mol
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

    params.calc_extn_coeff = True
    params.calc_pI_fasta = False
    params.calc_pIChemiSt = True
    mol = Chem.MolFromSmiles(smiles_str)
    return {
        1: {
            "mol_name": params.mol_name,
            "mol_obj": mol,
            "fasta": get_fasta_from_mol(mol),
        }
    }
