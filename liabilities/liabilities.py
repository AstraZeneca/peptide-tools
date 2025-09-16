import json
import os

from rdkit import Chem

script_dir = os.path.dirname(os.path.realpath(__file__))
liabilities_data = os.path.join(script_dir, "liabilities.json")
with open(liabilities_data) as f:
    liability_defs = json.load(f)


def calculate_liabilities_from_mol(mol):
    """
    This function calculates the sequence liability based on the SMARTS patterns
    of specific amino acid degradation pathways. The RDKit molecule is checked
    for matching SMARTS patterns, and liabilities are flagged accordingly.
    It also counts how many times each liability pattern is matched.

    :param mol: RDKit molecule object
    :return: Dictionary of sequence liabilities with SMARTS patterns,
        pH ranges, references, and counts
    """

    # Initialize a dictionary to store detected liabilities
    detected_liabilities = {}

    # Check for substructure matches with each SMARTS pattern
    for liability in liability_defs:
        smarts_pattern = liability["smarts"]
        smarts_mol = Chem.MolFromSmarts(smarts_pattern)

        # Get the list of all matches for the SMARTS pattern
        matches = mol.GetSubstructMatches(smarts_mol)
        match_count = len(matches)  # Count the number of matches

        # If there are matches, add the liability to the detected liabilities dictionary
        if match_count > 0:
            # correct for the case of symmetric smars patters
            if liability["liability"] == "Oxidation of thioethers":
                match_count = match_count / 2

            detected_liabilities[liability["liability"]] = {
                "liability": liability["liability"],
                "fasta_pattern": liability["fasta_pattern"],
                "match_count": match_count,
                "factors_favoring_stability": liability["factors_favoring_stability"],
                "factors_disfavoring_stability": liability[
                    "factors_disfavoring_stability"
                ],
                "comments": liability["comments"],
                "references": liability["references"],
            }

    return detected_liabilities


def calculate_liabilities_from_fasta(fasta_sequence):
    """
    This function calculates the sequence liability based on the presence
    of specific amino acid sequences in the input peptide sequence
    and includes pH stability and instability ranges.

    :param fasta_sequence: Input peptide sequence as a string (FASTA format)
    :return: Dictionary of sequence liabilities with associated pH ranges
    """

    # Initialize a dictionary to store detected liabilities
    detected_liabilities = {}

    # Convert the sequence to uppercase to handle both L- and D- amino acids
    fasta_sequence_upper = fasta_sequence.upper()

    # Check each liability pattern against the sequence
    for liability in liability_defs:

        # N-terminal special case 1
        if liability["references"] == "Bersin, Patel, and Topp 2021":
            if liability["fasta_pattern"].upper() == fasta_sequence_upper[0]:
                detected_liabilities[liability["liability"]] = {
                    "liability": liability["liability"],
                    "fasta_pattern": liability["fasta_pattern"],
                    "match_count": 1,
                    "factors_favoring_stability": liability[
                        "factors_favoring_stability"
                    ],
                    "factors_disfavoring_stability": liability[
                        "factors_disfavoring_stability"
                    ],
                    "comments": liability["comments"],
                    "references": liability["references"],
                }

        # N-terminal special case 2
        elif liability["references"] in "Capasso and Mazzarella 1999":
            if liability["fasta_pattern"][1:2].upper() == fasta_sequence_upper[1:2]:
                detected_liabilities[liability["liability"]] = {
                    "liability": liability["liability"],
                    "fasta_pattern": liability["fasta_pattern"],
                    "match_count": 1,
                    "factors_favoring_stability": liability[
                        "factors_favoring_stability"
                    ],
                    "factors_disfavoring_stability": liability[
                        "factors_disfavoring_stability"
                    ],
                    "comments": liability["comments"],
                    "references": liability["references"],
                }

        # Skip the 2nd matching of the alternative tautomer
        # liability entry to avoid duplicated entries for Arg
        elif (
            liability["liability"]
            == "Oxidation (carbonylation) of Arg (alternative tautomer)"
        ):
            continue

        # General case
        elif liability["fasta_pattern"].upper() in fasta_sequence_upper:
            detected_liabilities[liability["liability"]] = {
                "liability": liability["liability"],
                "fasta_pattern": liability["fasta_pattern"],
                "match_count": fasta_sequence.count(liability["fasta_pattern"]),
                "factors_favoring_stability": liability["factors_favoring_stability"],
                "factors_disfavoring_stability": liability[
                    "factors_disfavoring_stability"
                ],
                "comments": liability["comments"],
                "references": liability["references"],
            }

    return detected_liabilities


def calc_liabilities_from_dict(input_dict, input_type="structure"):
    """
    Calculates liabilities from an input dictionary.

    Parameters:
    input_dict (dict): A dictionary containing molecular data.
    input_type (str): Defines the input type, either "structure"
        (RDKit Molecule) or "fasta" (Amino Acid Sequence).

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
            matches = calculate_liabilities_from_mol(mol_object)

            # Store the result in the output dictionary
            dict_output[mol_idx] = {"mol_name": mol_name, "liabilities": matches}

        elif fasta:
            mol_name = mol_data.get("mol_name", "Unnamed Sequence")
            fasta = mol_data.get("fasta")

            if fasta is None:
                raise ValueError(f"FASTA sequence for {mol_name} is missing.")

            # Calculate liabilities based on FASTA sequence
            matches = calculate_liabilities_from_fasta(fasta)

            # Store the result in the output dictionary
            dict_output[mol_idx] = {"mol_name": mol_name, "liabilities": matches}
    return dict_output
