from rdkit import Chem

# Define the liabilities with SMARTS patterns, their corresponding liabilities, pH ranges, and references
liabilities = [
        {
            "smarts": "N[CH1]([CH2][CH0](=[OH0])[OH1,O-])[CH0](=[OH0])[NH0]1[CH2][CH2][CH2][CH1]1[CH0](=[OH0])",
            "fasta_pattern":"DP",
            "liability": "Asp-Pro bond cleavage",
            "pH_stability_range": "6.0–7.5",
            "pH_instability_range": "<5.0 or >8.0",
            "references": "Manning et al. (2010); Cleland et al. (1993)"
        },
        {
            "smarts": "[OH0]=[CH0][CH1](N)([CH2][CH2][SX2][CH3])",
            "fasta_pattern":"M",
            "liability": "Methionine oxidation",
            "pH_stability_range": "6.0–7.5",
            "pH_instability_range": "<5.0 or >8.0",
            "references": "Manning et al. (2010); Cleland et al. (1993); Al Musaimi et al. (2022)"
        },
        {
            "smarts": "N[CH1]([CH2][CH0](=[OH0])[NH2])[CH0](=[OH0])",  
            "fasta_pattern":"N",
            "liability": "Deamidation of Asn",
            "pH_stability_range": "4.5–6.0",
            "pH_instability_range": "<4.0 or >6.0",
            "references": "Wakankar & Borchardt (2006)"
        },

        {
            "smarts": "N[CH1]([CH2][CH0](=[OH0])[NH2])[CH0](=[OH0])[NH1][CH2][CH0](=[OH0])",
            "fasta_pattern":"NG",
            "liability": "Deamidationn of Asn-Gly",
            "pH_stability_range": "4.5–6.5",
            "pH_instability_range": "<4.0 or >7.0",
            "references": "Cleland et al. (1993); Li et al. (1995); Al Musaimi et al. (2022)"
        },

        {
            "smarts": "C[CH2][SX2][SX2][CH2]C",  
            "fasta_pattern":"C",
            "liability": "Disulfide bond scrambling",
            "pH_stability_range": "5.5–8.0",
            "pH_instability_range": "<4.0 or >9.0",
            "references": "Wang & Roberts (2018)"
        },

        {
            "smarts": "[CH3][CH1]([OH1])[CH1]([NH1][CH0](=[OH0])[CH1](N)[CH2][OH1])[CH0](=[OH0])", 
            "fasta_pattern":"ST",
            "liability": "Hydrolysis at Ser-Thr bonds",
            "pH_stability_range": "7.0–8.5",
            "pH_instability_range": "<6.0 or >9.0",
            "references": "Shire et al. (2004)"
        },

        {
            "smarts": "[NH2,N+H3][CH2][CH2][CH2][CH2][CH1](N)[CH0](=[OH0])",
            "liability": "Lys oxidation",
            "fasta_pattern":"K",
            "pH_stability_range": "4.5–6.5",
            "pH_instability_range": "<4.0 or >7.0",
            "references": "Cleland et al. (1993); Li et al. (1995); Al Musaimi et al. (2022)"
        },

        {
            "smarts": "[CH0](=[OH0])[CH1]1[CH2][CH2][CH2]N1",
            "liability": "Hydrolysis near Pro residues",
            "fasta_pattern":"P",
            "pH_stability_range": "6.0–7.0",
            "pH_instability_range": "<5.0 or >8.5",
            "references": "Shire et al. (2004); Nugrahadi et al. (2003)"
        },

        {
            "smarts": "[NH2,N+H3][CH1]([CH2][cH0]1[cH1][nH,n,n+H][cH1][nH,n,n+H]1)[CH0](=[OH0])",
            "liability": "N-terminal His oxidation",
            "fasta_pattern":"H",
            "pH_stability_range": "6.5–8.0",
            "pH_instability_range": "<6.0 or >8.5",
            "references": "Manning et al. (2010); Li et al. (1995); Al Musaimi et al. (2022)"
        }
        
]
    

def calculate_liabilities_from_mol(mol):
    """
    This function calculates the sequence liability based on the SMARTS patterns 
    of specific amino acid degradation pathways. The RDKit molecule is checked 
    for matching SMARTS patterns, and liabilities are flagged accordingly.
    It also counts how many times each liability pattern is matched.

    :param mol: RDKit molecule object
    :return: Dictionary of sequence liabilities with SMARTS patterns, pH ranges, references, and counts
    """

    # Initialize a dictionary to store detected liabilities
    detected_liabilities = {}

    # Check for substructure matches with each SMARTS pattern
    for liability in liabilities:
        smarts_pattern = liability["smarts"]
        smarts_mol = Chem.MolFromSmarts(smarts_pattern)
        
        # Get the list of all matches for the SMARTS pattern
        matches = mol.GetSubstructMatches(smarts_mol)
        match_count = len(matches)  # Count the number of matches
        
        # If there are matches, add the liability to the detected liabilities dictionary
        if match_count > 0:
            detected_liabilities[liability["liability"]] = {
                "liability": liability["liability"],
                "fasta_pattern": liability["fasta_pattern"],
#                "SMARTS_pattern": smarts_pattern,
                "match_count": match_count,
                "pH_stability_range": liability["pH_stability_range"],
                "pH_instability_range": liability["pH_instability_range"],
                "references": liability["references"]
            }

    return detected_liabilities


# # Example usage with RDKit
# from rdkit import Chem

# # Example peptide molecule (you can replace this with any RDKit peptide)
# smiles_string = "CC(=O)NC(=O)C"  # Simple example peptide (Acetyl-Glycine)

# mol = Chem.MolFromSmiles(smiles_string)

# # Calculate liabilities
# liabilities = calculate_liabilities_from_rdkit_mol(mol)

# # Output the identified liabilities
# if liabilities:
#     for liability, details in liabilities.items():
#         print(f"Liability: {liability}")
#         print(f"SMARTS Pattern: {details['SMARTS_pattern']}")
#         print(f"Match Count: {details['match_count']}")
#         print(f"P-H Stability Range: {details['pH_stability_range']}")
#         print(f"P-H Instability Range: {details['pH_instability_range']}")
#         print(f"References: {details['references']}")
# else:
#     print("No liabilities detected in the sequence.")



def calculate_liabilities_from_fasta(fasta_sequence):
    """
    This function calculates the sequence liability based on the presence of specific amino acid sequences
    in the input peptide sequence and includes pH stability and instability ranges.

    :param fasta_sequence: Input peptide sequence as a string (FASTA format)
    :return: Dictionary of sequence liabilities with associated pH ranges
    """

    # Initialize a dictionary to store detected liabilities
    detected_liabilities = {}

    # Convert the sequence to uppercase to handle both L- and D- amino acids
    fasta_sequence_upper = fasta_sequence.upper()

    # Check each liability pattern against the sequence
    for liability in liabilities:
        # case of N-terminal histidine
        if liability["liability"] == "N-terminal His oxidation":
            if liability["fasta_pattern"].upper() == fasta_sequence_upper[0]:
                detected_liabilities[liability["liability"]] = {
                    "liability": liability["liability"],
                    "fasta_pattern": liability["fasta_pattern"],
                    "match_count": 1,
                    "pH_stability_range": liability["pH_stability_range"],
                    "pH_instability_range": liability["pH_instability_range"],
                    "references": liability["references"]
                }
        # general case
        else:
            if liability["fasta_pattern"].upper() in fasta_sequence_upper:
                detected_liabilities[liability["liability"]] = {
                    "liability": liability["liability"],
                    "fasta_pattern": liability["fasta_pattern"],
                    "match_count": fasta_sequence.count(liability["fasta_pattern"]),
                    "pH_stability_range": liability["pH_stability_range"],
                    "pH_instability_range": liability["pH_instability_range"],
                    "references": liability["references"]
                }

    return detected_liabilities


# if __name__ == "__main__"
#     # Example usage
#     fasta_sequence = "MKTPIALDEKDPQCVLNTG"  # Example peptide sequence

#     liabilities = calculate_liabilities_from_fasta(fasta_sequence)

#     # Output the identified liabilities
#     if liabilities:
#         for liability, details in liabilities.items():
#             print(f"Liability: {liability}")
#             print(f"P-H Stability Range: {details['pH_stability_range']}")
#             print(f"P-H Instability Range: {details['pH_instability_range']}")
#             print(f"References: {details['references']}")
#     else:
#         print("No liabilities detected in the sequence.")




def calc_liabilities_from_dict(input_dict, input_type="structure"):
    """
    Calculates liabilities from an input dictionary.
    
    Parameters:
    input_dict (dict): A dictionary containing molecular data.
    input_type (str): Defines the input type, either "structure" (RDKit Molecule) or "fasta" (Amino Acid Sequence).
    
    Returns:
    dict_output (dict): A dictionary with liabilities for each molecule or peptide.
    """
    dict_output = dict()
    for mol_idx, mol_data in input_dict.items():
        fasta = input_dict[mol_idx].get("fasta")
        mol = input_dict[mol_idx].get("mol_obj")

        # MOL overrides FASTA
        if mol:
            mol_name = mol_data.get("mol_name", "Unnamed Molecule")
            mol_object = mol_data.get("mol_obj")

            if mol_object is None:
                raise ValueError(f"Molecule object for {mol_name} is missing.")
            
            # Calculate liabilities based on SMARTS patterns
            liabilities = calculate_liabilities_from_mol(mol_object)
            
            # Store the result in the output dictionary
            dict_output[mol_idx] = {
                "mol_name": mol_name,
                "liabilities": liabilities
            }
    
        elif fasta:
            mol_name = mol_data.get("mol_name", "Unnamed Sequence")
            fasta = mol_data.get("fasta")
            
            if fasta is None:
                raise ValueError(f"FASTA sequence for {mol_name} is missing.")
                       
            # Calculate liabilities based on FASTA sequence
            liabilities = calculate_liabilities_from_fasta(fasta)
            
            # Store the result in the output dictionary
            dict_output[mol_idx] = {
                "mol_name": mol_name,
                "liabilities": liabilities
            }
    return dict_output

