from rdkit import Chem
from rdkit.Chem import Descriptors
from Bio.SeqUtils import molecular_weight
# from Bio.SeqUtils import gravy
from rdkit.Chem import Crippen

def calc_molecular_descriptors_from_dict(
    input_dict,
    input_type="structure",
):
    """
    Calculates the molecular weight and length of molecules (or peptides) from the input dictionary.
    
    Parameters:
    input_dict (dict): A dictionary containing molecular data.
    input_type (str): Defines the input type, either "structure" (RDKit Molecule) or "fasta" (Amino Acid Sequence).
    
    Returns:
    dict_output (dict): A dictionary with molecular weights and lengths for each molecule or peptide.
    """
    dict_output = dict()

    if input_type == "structure":
        # If the input is a structure (SMILES or RDKit molecule)
        for mol_idx in input_dict.keys():
            mol_name = input_dict[mol_idx].get("mol_name", "Unnamed Molecule")
            mol_object = input_dict[mol_idx].get("mol_obj")
            fasta = input_dict[mol_idx].get("fasta")

            if mol_object is None:
                raise ValueError(f"Molecule object for {mol_name} is missing.")
            
            # Calculate the molecular weight using RDKit
            mol_weight = Descriptors.MolWt(mol_object)

            # Calculate the GRAVY score
            # gravy_score = gravy(fasta)

            # Calculate ClogP (logP)
            logp = Crippen.MolLogP(mol_object)

            # Generate SMILES from RDKit mol object
            smiles = Chem.MolToSmiles(mol_object)

            # Store the result in the output dictionary
            dict_output[mol_idx] = {
                "mol_name": mol_name,
                "molecular_weight": mol_weight,
                "seq_length": len(fasta),  # Length of the molecule fasta
                # "gravy_score": gravy_score,  # GRAVY score
                "logp":logp,
                "fasta":fasta,
                "smiles":smiles                           
            }
    
    elif input_type == "fasta":
        # If the input is a fasta sequence (protein or peptide sequence)
        for mol_idx in input_dict.keys():
            mol_name = input_dict[mol_idx].get("mol_name", "Unnamed Sequence")
            fasta = input_dict[mol_idx].get("fasta")
            
            if fasta is None:
                raise ValueError(f"FASTA sequence for {mol_name} is missing.")
            
            # if nonnatural is given in fasta - skip it to allow further calcualtions. 
            if "X" in fasta:
                raise ValueError(f"FASTA sequence for {mol_name} contains X. Please skip it or replace by similar amio-acid.")
                #fasta = fasta.replace("X","")

            # if nonnatural is given in fasta - skip it to allow further calcualtions. 
            if not fasta.isalpha():
                raise ValueError(f"FASTA sequence for {mol_name} contains unknown letters. Please replace with known amio-acids.")
                #fasta = fasta.replace("X","")



            # Calculate molecular weight for the peptide sequence
            mol_weight = molecular_weight(fasta, seq_type='protein')
            seq_length = len(fasta)  # Length of the peptide sequence
            
            # Calculate the GRAVY score
            # gravy_score = gravy(fasta)

            mol_from_fasta = Chem.MolFromSequence(fasta)
            logp = Crippen.MolLogP(mol_from_fasta)

            # Generate SMILES from RDKit mol object
            smiles = Chem.MolToSmiles(mol_from_fasta)

            # Store the result in the output dictionary
            dict_output[mol_idx] = {
                "mol_name": mol_name,
                "molecular_weight": mol_weight,
                "seq_length": seq_length,  # Length of the peptide sequence
                # "gravy_score": gravy_score,  # GRAVY score
                "logp":logp,
                "fasta":fasta,
                "smiles":smiles      
            }
    
    else:
        raise ValueError("Unknown input_type. It must be either 'structure' or 'fasta'.")
    
    return dict_output


# def calculate_peptide_molecular_weight(fasta_sequence):
#     """
#     Calculates the molecular weight of a peptide sequence.
    
#     Parameters:
#     fasta_sequence (str): The peptide sequence in FASTA format (e.g., "ACDEFG").
    
#     Returns:
#     mol_weight (float): The molecular weight of the peptide.
#     """
#     # Amino acid molecular weights (average)
#     amino_acid_weights = {
#         'A': 71.08,  # Alanine
#         'C': 103.14, # Cysteine
#         'D': 115.09, # Aspartic acid
#         'E': 129.12, # Glutamic acid
#         'F': 147.18, # Phenylalanine
#         'G': 57.05,  # Glycine
#         'H': 137.14, # Histidine
#         'I': 113.16, # Isoleucine
#         'K': 128.17, # Lysine
#         'L': 113.16, # Leucine
#         'M': 131.19, # Methionine
#         'N': 114.11, # Asparagine
#         'P': 97.12,  # Proline
#         'Q': 128.13, # Glutamine
#         'R': 156.19, # Arginine
#         'S': 87.08,  # Serine
#         'T': 101.11, # Threonine
#         'V': 99.14,  # Valine
#         'W': 186.21, # Tryptophan
#         'Y': 163.18, # Tyrosine
#     }
    
#     mol_weight = 0.0
#     for aa in fasta_sequence:
#         if aa in amino_acid_weights:
#             mol_weight += amino_acid_weights[aa]
#         else:
#             raise ValueError(f"Invalid amino acid: {aa}")
    
#     return mol_weight
