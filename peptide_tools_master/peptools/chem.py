import os

from rdkit import Chem


class FileExtension:
    SDF = ".sdf"
    SMI = ".smi"


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


def get_fasta_from_mol(mol):
    from smi2scrambledfasta import get_scrambled_fasta_from_mol

    fasta = get_scrambled_fasta_from_mol(mol)
    if len(fasta) == 0:
        raise Exception("ERROR: returned fasta is empry. something is wrong. Exit")
        sys.exit(1)
    return fasta


def get_fasta_from_smiles(smi):
    from smi2scrambledfasta import get_scrambled_fasta_from_smiles

    fasta = get_scrambled_fasta_from_smiles(smi)
    if len(fasta) == 0:
        raise Exception("ERROR: returned fasta is empry. something is wrong. Exit")
        sys.exit(1)
    return fasta
