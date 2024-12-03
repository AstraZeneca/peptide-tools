import os

from peptools.chem import get_fasta_from_mol
from peptools.chem import get_fasta_from_smiles
from peptools.chem import read_structure_file
from rdkit import Chem


def generate_input(input):
    INPUT = input
    mol_name = "none"
    file_extension = None
    out_fext = None
    known_file_types = [".sdf", ".smi", ".smiles", ".fasta"]

    mol_supply_json = {}

    if os.path.exists(INPUT):
        # Assume input is a file.

        lPlot = False

        filename, file_extension = os.path.splitext(INPUT)
        if file_extension not in known_file_types:
            raise Exception(
                "Error! File extention not in supported file types:"
                + str(known_file_types)
            )
            sys.exit(1)
        else:
            # smi = ''
            # mol_name = 'none'
            inputFile = INPUT
            if (
                file_extension == ".smi"
                or file_extension == ".smiles"
                or file_extension == ".csv"
                or file_extension == ".fasta"
            ):
                out_fext = ".csv"
            elif file_extension == ".sdf":
                out_fext = ".sdf"

            outputFile = filename + "_OUTPUT" + out_fext

            if file_extension != ".fasta":
                l_calc_extn_coeff = True
                l_calc_pI_fasta = False
                l_calc_pIChemiSt = True
                mol_supply_json = read_structure_file(inputFile)
            else:
                l_calc_extn_coeff = True
                l_calc_pI_fasta = True
                l_calc_pIChemiSt = False
                mol_supply_json = read_fasta_file(inputFile)

    elif not INPUT.isalpha():
        # print("Input recognised as SMILES")
        # Assume it is smiles, if contains not only letters
        # print("Input is SMILES")
        mol_unique_ID = 1
        smi = INPUT
        # mol_name = 'none'
        fasta = get_fasta_from_smiles(smi)
        inputFile = ""
        outputFile = ""
        l_calc_extn_coeff = True
        l_calc_pI_fasta = False
        l_calc_pIChemiSt = True
        mol = Chem.MolFromSmiles(smi)
        mol_supply_json[mol_unique_ID] = {
            "mol_name": mol_name,
            "mol_obj": mol,
            "fasta": get_fasta_from_mol(mol),
        }

    elif INPUT.isalpha():
        # Assume it is FASTA, if contains only letters
        # print("Input is FASTA")
        mol_unique_ID = 1
        fasta = INPUT

        smi = ""
        # mol_name = 'none'
        inputFile = ""
        outputFile = ""
        l_calc_extn_coeff = True
        l_calc_pI_fasta = True
        l_calc_pIChemiSt = False
        mol_supply_json[mol_unique_ID] = {
            "mol_name": mol_name,
            "mol_obj": None,
            "fasta": fasta,
        }

    elif (
        (INPUT[0:2] == "AZ" and INPUT[2].isdigit())
        or (INPUT[0:2] == "SN" and INPUT[2].isdigit())
        or (INPUT[0:4] == "MEDI" and INPUT[4].isdigit())
    ):
        # A database ID given
        # print("Input is a database ID")
        mol_unique_ID = 1
        smi = get_smiles_from_dbid(INPUT)
        if len(smi) == 0:
            raise Exception(
                "ERROR: could not convert database ID to smiles. Is it corret ID from supported DBs? Exit."
            )
            sys.exit(1)
        mol_name = INPUT
        fasta = get_fasta_from_smiles(smi)
        inputFile = ""
        outputFile = ""
        l_calc_extn_coeff = True
        l_calc_pI_fasta = False
        l_calc_pIChemiSt = True
        mol = Chem.MolFromSmiles(smi)
        mol_supply_json[mol_unique_ID] = {
            "mol_name": mol_name,
            "mol_obj": mol,
            "fasta": fasta,
        }

    else:
        raise Exception(
            "ERROR: input not recongnized: not smiles, not fasta, not a known databaase ID. Must be a bug. Contact developer. Exit."
        )
        sys.exit(1)
    return (
        mol_supply_json,
        inputFile,
        outputFile,
        l_calc_extn_coeff,
        l_calc_pI_fasta,
        l_calc_pIChemiSt,
        file_extension,
        out_fext,
    )


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
