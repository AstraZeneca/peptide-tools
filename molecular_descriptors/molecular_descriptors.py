from Bio.SeqUtils import molecular_weight
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors

from smi2scrambledfasta import get_scrambled_fasta_from_mol


def calc_molecular_descriptors_from_dict(input_dict):
    """
    Calculates the molecular weight and length of
    molecules (or peptides) from the input dictionary.

    Parameters:
    input_dict (dict): A dictionary containing molecular data.
    input_type (str): Defines the input type, either "structure"
        (RDKit Molecule) or "fasta" (Amino Acid Sequence).

    Returns:
    dict_output (dict): A dictionary with molecular weights and
        lengths for each molecule or peptide.
    """
    dict_output = dict()
    for mol_idx, mol_data in input_dict.items():
        fasta = input_dict[mol_idx].get("fasta")
        mol = input_dict[mol_idx].get("mol_obj")
        mol_name = input_dict[mol_idx].get("mol_name", "unnamed molecule")

        # MOL overrides FASTA
        if mol:
            fasta = get_scrambled_fasta_from_mol(mol)

            # Calculate the molecular weight using RDKit
            mol_weight = round(Descriptors.MolWt(mol), 3)

            # Calculate ClogP (logP)
            logp = round(Crippen.MolLogP(mol), 3)

            # Generate SMILES from RDKit mol object
            smiles = Chem.MolToSmiles(mol)

            # Store the result in the output dictionary
            dict_output[mol_idx] = {
                "mol_name": mol_name,
                "molecular_weight": mol_weight,
                "seq_length": len(fasta),  # Length of the molecule fasta
                "logp": logp,
                "fasta": fasta,
                "smiles": smiles,
            }

        elif fasta:
            fasta = input_dict[mol_idx].get("fasta")

            if fasta is None:
                raise ValueError(f"FASTA sequence for {mol_name} is missing.")

            # if nonnatural is given in fasta - skip it to allow further calcualtions.
            if "X" in fasta:
                raise ValueError(
                    f"FASTA sequence for {mol_name} contains X. "
                    "Please skip it or replace by similar amio-acid."
                )
                # fasta = fasta.replace("X","")

            # if nonnatural is given in fasta - skip it to allow further calcualtions.
            if not fasta.isalpha():
                raise ValueError(
                    f"FASTA sequence for {mol_name} contains unknown letters. "
                    "Please replace with known amio-acids."
                )
                # fasta = fasta.replace("X","")

            # Calculate molecular weight for the peptide sequence
            mol_weight = round(molecular_weight(fasta, seq_type="protein"), 3)
            seq_length = len(fasta)  # Length of the peptide sequence

            mol_from_fasta = Chem.MolFromSequence(fasta)
            logp = round(Crippen.MolLogP(mol_from_fasta), 3)

            # Generate SMILES from RDKit mol object
            smiles = Chem.MolToSmiles(mol_from_fasta)

            # Store the result in the output dictionary
            dict_output[mol_idx] = {
                "mol_name": mol_name,
                "molecular_weight": mol_weight,
                "seq_length": seq_length,  # Length of the peptide sequence
                "logp": logp,
                "fasta": fasta,
                "smiles": smiles,
            }

    return dict_output
