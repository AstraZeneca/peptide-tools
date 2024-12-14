import argparse
import os
import sys

from pichemist.fasta.scrambler import FastaScrambler
from pichemist.molecule import PeptideCutter
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import Recap


class FastaException(Exception):
    pass


def get_scrambled_fasta_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return get_scrambled_fasta_from_mol(mol)


def get_scrambled_fasta_from_mol(mol):
    frags_smi_list = PeptideCutter().break_amide_bonds_and_cap(mol)
    fasta = FastaScrambler().get_scrambled_fasta_from_list(frags_smi_list)
    if not fasta:
        raise FastaException("An empty FASTA was generated.")
    return fasta


def run_main():

    ### Parser
    parser = argparse.ArgumentParser(
        description="Script converts smiles into scrambled fasta. Scrambled mean that the order of aminoacids is not preserved. All non-natural amino acids decoded as X."
    )
    parser.add_argument(
        "-s",
        dest="smiles",
        help="input smiles. if used then single smi is assumed and fasta returned in stdout. filenames are ignored",
        default="",
    )
    parser.add_argument(
        "-i",
        dest="inputFile",
        help="input file with molecule structure. smi or sdf",
        default="",
    )
    parser.add_argument(
        "-o",
        dest="outputFile",
        help="output file with molecule structure. fasta",
        default="",
    )
    args = parser.parse_args()

    # Get options
    if len(args.smiles) != 0:
        # assume smiles input
        suppl = [Chem.MolFromSmiles(args.smiles)]
    elif len(args.inputFile) != 0:
        # Assume filename as input
        inputFile = args.inputFile
        filename, ext = os.path.splitext(inputFile)

        # Initialize file reader
        if ext == ".sdf":
            suppl = Chem.SDMolSupplier(inputFile)
        elif ext == ".smi":
            suppl = Chem.SmilesMolSupplier(inputFile, titleLine=False)
        else:
            raise Exception(
                "!Warning: extension of file is not smi or sdf. Assume it is smi. Continue. "
            )
            suppl = Chem.SmilesMolSupplier(inputFile, titleLine=False)

    else:
        raise Exception("Error: either smiles or input file should be given. Exit. ")
        sys.exit(1)

    ### Run conversion
    fasta_list = []
    for mol in suppl:
        fasta = get_scrambled_fasta_from_mol(mol)
        fasta_list.append("".join(fasta))

    ### Output
    if len(args.smiles) != 0:
        # Write out fasta as string into stdout
        for fasta in fasta_list:
            print("".join(fasta))

    elif len(args.outputFile) != 0:
        # Write into output file
        with open(args.outputFile, "w") as f:
            for i, fasta in enumerate(fasta_list):
                f.write("%s    %s\n" % (fasta, "tmpname" + str(i + 1)))
    else:
        raise Exception(
            "Error: cant decide where ti output. Must be a n error in input. Exit."
        )

    return


if __name__ == "__main__":
    run_main()
