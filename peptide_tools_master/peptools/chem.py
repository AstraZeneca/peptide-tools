import os

from rdkit import Chem


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
