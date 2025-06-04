from rdkit import Chem

# Define the liabilities with SMARTS patterns, their corresponding liabilities, pH ranges, and references
# liabilities = [
#     {
#         "smarts": "N[CH1]([CH2][CH0](=[OH0])[OH1,O-])[CH0](=[OH0])[NH0]1[CH2][CH2][CH2][CH1]1[CH0](=[OH0])",
#         "fasta_pattern": "DP",
#         "liability": "Asp-Pro bond cleavage",
#         "pH_stability_range": "6.0–7.5",
#         "pH_instability_range": "<5.0 or >8.0",
#         "references": "Manning et al. (2010); Cleland et al. (1993)",
#     },
#     {
#         "smarts": "[OH0]=[CH0][CH1](N)([CH2][CH2][SX2][CH3])",
#         "fasta_pattern": "M",
#         "liability": "Methionine oxidation",
#         "pH_stability_range": "6.0–7.5",
#         "pH_instability_range": "<5.0 or >8.0",
#         "references": "Manning et al. (2010); Cleland et al. (1993); Al Musaimi et al. (2022)",
#     },
#     {
#         "smarts": "N[CH1]([CH2][CH0](=[OH0])[NH2])[CH0](=[OH0])",
#         "fasta_pattern": "N",
#         "liability": "Deamidation of Asn",
#         "pH_stability_range": "4.5–6.0",
#         "pH_instability_range": "<4.0 or >6.0",
#         "references": "Wakankar & Borchardt (2006)",
#     },
#     {
#         "smarts": "N[CH1]([CH2][CH0](=[OH0])[NH2])[CH0](=[OH0])[NH1][CH2][CH0](=[OH0])",
#         "fasta_pattern": "NG",
#         "liability": "Deamidationn of Asn-Gly",
#         "pH_stability_range": "4.5–6.5",
#         "pH_instability_range": "<4.0 or >7.0",
#         "references": "Cleland et al. (1993); Li et al. (1995); Al Musaimi et al. (2022)",
#     },
#     {
#         "smarts": "C[CH2][SX2][SX2][CH2]C",
#         "fasta_pattern": "C",
#         "liability": "Disulfide bond scrambling",
#         "pH_stability_range": "5.5–8.0",
#         "pH_instability_range": "<4.0 or >9.0",
#         "references": "Wang & Roberts (2018)",
#     },
#     {
#         "smarts": "[CH3][CH1]([OH1])[CH1]([NH1][CH0](=[OH0])[CH1](N)[CH2][OH1])[CH0](=[OH0])",
#         "fasta_pattern": "ST",
#         "liability": "Hydrolysis at Ser-Thr bonds",
#         "pH_stability_range": "7.0–8.5",
#         "pH_instability_range": "<6.0 or >9.0",
#         "references": "Shire et al. (2004)",
#     },
#     {
#         "smarts": "[NH2,N+H3][CH2][CH2][CH2][CH2][CH1](N)[CH0](=[OH0])",
#         "liability": "Lys oxidation",
#         "fasta_pattern": "K",
#         "pH_stability_range": "4.5–6.5",
#         "pH_instability_range": "<4.0 or >7.0",
#         "references": "Cleland et al. (1993); Li et al. (1995); Al Musaimi et al. (2022)",
#     },
#     {
#         "smarts": "[CH0](=[OH0])[CH1]1[CH2][CH2][CH2]N1",
#         "liability": "Hydrolysis near Pro residues",
#         "fasta_pattern": "P",
#         "pH_stability_range": "6.0–7.0",
#         "pH_instability_range": "<5.0 or >8.5",
#         "references": "Shire et al. (2004); Nugrahadi et al. (2003)",
#     },
#     {
#         "smarts": "[NH2,N+H3][CH1]([CH2][cH0]1[cH1][nH,n,n+H][cH1][nH,n,n+H]1)[CH0](=[OH0])",
#         "liability": "N-terminal His oxidation",
#         "fasta_pattern": "H",
#         "pH_stability_range": "6.5–8.0",
#         "pH_instability_range": "<6.0 or >8.5",
#         "references": "Manning et al. (2010); Li et al. (1995); Al Musaimi et al. (2022)",
#     },
# ]



#liabilities = [
#       {
#           "smarts": "[NH2,N+H3][CH1]([CH2][CH0](=[OH0])[OH1,O-])[CH0](=[OH0])[NH0]1[CH2][CH2][CH2][CH1]1[CH0](=[OH0])",
#           "fasta_pattern": "DP (at N-terminus)",
#           "liability": "Hydrolysis of Asp-Pro bond at N-terminus",
#           "ph_stability_range": "4-6",
#           "other_parameters": "Temperature",
#           "references": "Cleland, Powell, and Shire 1993; Manning et al. 2010"
#       },
#       {
#           "smarts": "N[CH1]([CH2][CH0](=[OH0])[NH2])[CH0](=[OH0])",
#           "fasta_pattern": "N",
#           "liability": "Deamidation of Asn",
#           "ph_stability_range": "4-6",
#           "other_parameters": "Temperature",
#           "references": "Wakankar and Borchardt 2006"
#       },
#       {
#           "smarts": "N[CH1]([CH2][CH2][CH0](=[OH0])[NH2])[CH0](=[OH0])",
#           "fasta_pattern": "Q",
#           "liability": "Deamination of Gln",
#           "ph_stability_range": "4-6",
#           "other_parameters": "Temperature",
#           "references": "Ref?"
#       },
#       {
#           "smarts": "N[CH1]([CH2][CH0](=[OH0])[NH2])[CH0](=[OH0])[NH1][CH2][CH0](=[OH0])",
#           "fasta_pattern": "NG",
#           "liability": "Deamidation of Asn in Asn-Gly",
#           "ph_stability_range": "4-5",
#           "other_parameters": "Temperature",
#           "references": "Al Musaimi et al. 2022; Cleland, Powell, and Shire 1993; Li, Schöneich, and Borchardt 1995; Shahrokh et al. 1994"
#       },
#       {
#           "smarts": "C[CH2][SX2][SX2][CH2]C",
#           "fasta_pattern": "C",
#           "liability": "Cys-Cys disulfide bond scrambling & ß-elimination",
#           "ph_stability_range": "6-8",
#           "other_parameters": "Metal ions, oxygen",
#           "references": "Gilbert 1990; Niu et al. 2016; Wang and Roberts 2018"
#       },
#       {
#           "smarts": "[CH3][CH1]([OH1])[CH1]([NH1][CH0](=[OH0])[CH1](N)[CH2][OH1])[CH0](=[OH0])",
#           "fasta_pattern": "ST",
#           "liability": "Hydrolysis at Ser-Thr bonds",
#           "ph_stability_range": "6-7",
#           "other_parameters": "Temperature",
#           "references": "Nick Pace et al. 2004; Shire, Shahrokh, and Liu 2004"
#       },
#       {
#           "smarts": "[NH2,N+H3][CH2][CH2][CH2][CH2][CH1](N)[CH0](=[OH0])",
#           "fasta_pattern": "K",
#           "liability": "Oxidation of Lys",
#           "ph_stability_range": "4-6",
#           "other_parameters": "-",
#           "references": "Al Musaimi et al. 2022; Cleland, Powell, and Shire 1993; Li, Schöneich, and Borchardt 1995"
#       },
#       {
#           "smarts": "[CH0](=[OH0])[CH1]1[CH2][CH2][CH2]N1",
#           "fasta_pattern": "P",
#           "liability": "Hydrolysis near Pro residue",
#           "ph_stability_range": "5-7",
#           "other_parameters": "Temperature",
#           "references": "Nugrahadi et al. 2023; Shire, Shahrokh, and Liu 2004"
#       },
#       {
#           "smarts": "[NH2,N+H3][CH1]([CH2][cH0]1[cH1][nH,n,n+H][cH1][nH,n,n+H]1)[CH0](=[OH0])",
#           "fasta_pattern": "H (at N-terminus)",
#           "liability": "Oxidation of N-terminal His",
#           "ph_stability_range": "5-6",
#           "other_parameters": "Temperature, oxygen and metal ions",
#           "references": "Al Musaimi et al. 2022; Li, Schöneich, and Borchardt 1995; Manning et al. 2010"
#       },
#       {
#           "smarts": "N[CH1]([CH2][cH0]1[cH1][nH1][cH0]2[cH1][cH1][cH1][cH1][cH0]12)[CH0](=[OH0])",
#           "fasta_pattern": "W",
#           "liability": "Oxidation of Trp",
#           "ph_stability_range": "4-6",
#           "other_parameters": "Temperature, light",
#           "references": "Li, Schöneich, and Borchardt 1995"
#       },
#       {
#           "smarts": "N[CH1]([CH2][cH0]1[cH1][cH1][cH1]([OH1,O-H0])[cH1][cH1]1)[CH0](=[OH0])",
#           "fasta_pattern": "Y",
#           "liability": "Oxidation of Tyr",
#           "ph_stability_range": "5-7",
#           "other_parameters": "Metal ions, light, oxygen",
#           "references": "Li, Schöneich, and Borchardt 1995; Wang 1999"
#       },
#       {
#           "smarts": "[C][CH2][SX2,S-X1]",
#           "fasta_pattern": "C",
#           "liability": "Oxidation of Cys & thiols",
#           "ph_stability_range": "6-7",
#           "other_parameters": "Metal ions, light, oxygen",
#           "references": "Stadtman and Levine 2003"
#       },
#       {
#           "smarts": "C[CH2][SX2][CH3]",
#           "fasta_pattern": "M",
#           "liability": "Oxidation of Met",
#           "ph_stability_range": "5-7",
#           "other_parameters": "Temperature, oxygen, metal ions",
#           "references": "Al Musaimi et al. 2022; Cleland, Powell, and Shire 1993; Manning et al. 2010"
#       },
#       {
#           "smarts": "C[CH2][SX2][CH2]C",
#           "fasta_pattern": "-",
#           "liability": "Oxidation of thioethers",
#           "ph_stability_range": "5-7",
#           "other_parameters": "Temperature, oxygen, metal ions",
#           "references": "Al Musaimi et al. 2022; Cleland, Powell, and Shire 1993; Manning et al. 2010"
#       }
#]



liabilities = [
    {
        "liability": "Deamidation of Asn",
        "smarts": "N[CH1]([CH2][CH0](=[OH0])[NH2])C(=[OH0])",
        "fasta_pattern": "N",
        "factors_favoring_stability": "pH 3-6",
        "factors_disfavoring_stability": (
            "Disordered secondary structure & high conformational flexibility, "
            "solvent accessibility, carboxylate buffers: phosphate, bicarbonate, tris"
        ),
        "comments": (
            "Deamidation rate strongly depends on the i+1 residue: very high risk for NG; "
            "high risk for NH, NS, NT, NA; lower risk for bulky and hydrophobic residues "
            "(e.g. NI, NV, NL, NW, NP, etc.)"
        ),
        "references": "N. e. Robinson et al. 2004; N. E. Robinson and Robinson 2004; Wakankar and Borchardt 2006"
    },
    {
        "liability": "Deamidation of Gln",
        "smarts": "N[CH1]([CH2][CH2][CH0](=[OH0])[NH2])C(=[OH0])",
        "fasta_pattern": "Q",
        "factors_favoring_stability": "pH 3-6",
        "factors_disfavoring_stability": (
            "Disordered secondary structure & high conformational flexibility, "
            "solvent accessibility, carboxylate buffers: phosphate, bicarbonate, tris"
        ),
        "comments": (
            "Deamidation rate strongly depends on the i+1 residue: high risk for QG; "
            "medium risk for QC; low risk for bulky and acidic residues "
            "(e.g. QW, QF, QD, QE, QP, etc.)"
        ),
        "references": "N. e. Robinson et al. 2004; N. E. Robinson and Robinson 2004; Wakankar and Borchardt 2006"
    },
    {
        "liability": "β-isomerization of Asp",
        "smarts": "N[CH1]([CH2][CH0](=[OH0])[OH,O-])C(=[OH0])",
        "fasta_pattern": "D",
        "factors_favoring_stability": "pH 6-8, addition of electrolyte",
        "factors_disfavoring_stability": (
            "Disordered secondary structure & high conformational flexibility, solvent accessibility, "
            "carboxylate buffers: phosphate, bicarbonate, addition of glycerol, sucrose"
        ),
        "comments": "High risk",
        "references": "Geiger and Clarke 1987; Wakankar and Borchardt 2006"
    },
    {
        "liability": "Racemization of Ser",
        "smarts": "N[CH1]([CH2][OH1])C(=[OH0])",
        "fasta_pattern": "S",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": (
            "Disordered secondary structure & high conformational flexibility, solvent accessibility"
        ),
        "comments": "-",
        "references": "Cloos and Christgau 2002"
    },
    {
        "liability": "Racemization of Ala",
        "smarts": "N[CH1]([CH3])C(=[OH0])",
        "fasta_pattern": "A",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": (
            "Disordered secondary structure & high conformational flexibility, solvent accessibility"
        ),
        "comments": "-",
        "references": "Cloos and Christgau 2002"
    },
    {
        "liability": "Backbone hydrolysis near Asp",
        "smarts": "N[CH1]([CH2][CH0](=[OH0])[OH,O-])C(=[OH0])",
        "fasta_pattern": "D",
        "factors_favoring_stability": "pH >5",
        "factors_disfavoring_stability": (
            "Disordered secondary structure & high conformational flexibility, high solvent accessibility"
        ),
        "comments": "Hydrolysis rate of DP motif is 8-10 times higher than the rest.",
        "references": "Lauer et al. 2016; Vlasak and Ionescu 2011"
    },
    {
        "liability": "Backbone hydrolysis near Ser & racemization of Ser",
        "smarts": "N[CH1]([CH2][OH1])C(=[OH0])",
        "fasta_pattern": "S",
        "factors_favoring_stability": "Transition metals",
        "factors_disfavoring_stability": "-",
        "comments": "-",
        "references": "Lyons, Jamie, and Truscott 2011, 2014; Wezynfeld, Frączyk, and Bal 2016"
    },
    {
        "liability": "Backbone hydrolysis near Thr",
        "smarts": "N[CH1]([CH]([CH3])[OH])C(=[OH0])",
        "fasta_pattern": "T",
        "factors_favoring_stability": "Transition metals",
        "factors_disfavoring_stability": "-",
        "comments": "-",
        "references": "Lyons, Jamie, and Truscott 2011, 2014; Wezynfeld, Frączyk, and Bal 2016"
    },
    {
        "liability": "Backbone hydrolysis between sequence positions 2 and 3 via diketopiperazine formation",
        "smarts": "[NH2,NH3+]C[CH0](=[OH0])[NH0]1[CH2][CH2][CH2][CH1]1[CH0](=[OH0])[NH1][CH2]C(=[OH0])",
        "fasta_pattern": "XPG",
        "factors_favoring_stability": "pH < ~5",
        "factors_disfavoring_stability": "Polar or charged amino acids in the 1st position, alkylammonium carboxylates",
        "comments": "Low risk",
        "references": "Capasso and Mazzarella 1999; Goolcharran and Borchardt 1998; Zhang, Hales, and Clemmer 2022"
    },
    {
        "liability": "Cyclization of N-terminal Gln (pGlu formation)",
        "smarts": "[NH2,NH3+][CH1]([CH2][CH2][CH0](=[OH0])[NH2])C(=[OH0])",
        "fasta_pattern": "Q",
        "factors_favoring_stability": "pH 5.5 – 6.5",
        "factors_disfavoring_stability": "-",
        "comments": "-",
        "references": "Bersin, Patel, and Topp 2021; Chelius et al. 2006"
    },
    {
        "liability": "Oxidation (carbonylation) of Lys",
        "smarts": "N[CH1]([CH2][CH2][CH2][CH2][NH2,NH3+])C(=[OH0])",
        "fasta_pattern": "K",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": (
            "Transition metal ions, light, ascorbic acid, proximity to RKPT-enriched regions, solvent exposure"
        ),
        "comments": "-",
        "references": "Kryndushkin and Rao 2016; Maisonneuve et al. 2009; Yang et al. 2014"
    },
    {
        "liability": "Oxidation (carbonylation) of Arg",
        "smarts": "N[CH1]([CH2][CH2][CH2][NH][CH0](=[NH1,NH2+])[NH2])C(=[OH0])",
        "fasta_pattern": "R",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": (
            "Transition metal ions, light, ascorbic acid, proximity to RKPT-enriched regions, solvent exposure"
        ),
        "comments": "-",
        "references": "Kryndushkin and Rao 2016; Maisonneuve et al. 2009; Yang et al. 2014"
    },
    {
        "liability": "Oxidation (carbonylation) of Arg (alternative tautomer)",
        "smarts": "N[CH1]([CH2][CH2][CH2][NH0,NH1+]=[CH0]([NH2])[NH2])C(=[OH0])",
        "fasta_pattern": "R",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": (
            "Transition metal ions, light, ascorbic acid, proximity to RKPT-enriched regions, solvent exposure"
        ),
        "comments": "-",
        "references": "Kryndushkin and Rao 2016; Maisonneuve et al. 2009; Yang et al. 2014"
    },
    {
        "liability": "Oxidation (carbonylation) of Pro",
        "smarts": "N1[CH2][CH2][CH2][CH1]1C(=[OH0])",
        "fasta_pattern": "P",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": (
            "Transition metal ions, light, ascorbic acid, proximity to RKPT-enriched regions, solvent exposure"
        ),
        "comments": "-",
        "references": "Kryndushkin and Rao 2016; Maisonneuve et al. 2009; Yang et al. 2014"
    },
    {
        "liability": "Oxidation (carbonylation) of Thr",
        "smarts": "N[CH1]([CH]([CH3])[OH])C(=[OH0])",
        "fasta_pattern": "T",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": (
            "Transition metal ions, light, ascorbic acid, proximity to RKPT-enriched regions, solvent exposure"
        ),
        "comments": "-",
        "references": "Kryndushkin and Rao 2016; Maisonneuve et al. 2009; Yang et al. 2014"
    },
    {
        "liability": "Oxidation of Phe",
        "smarts": "N[CH1]([CH2][cH0]1[cH1][cH1][cH1][cH1][cH1]1)C(=[OH0])",
        "fasta_pattern": "F",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": "Light",
        "comments": "-",
        "references": "Kerwin and Remmele 2007"
    },
    {
        "liability": "Oxidation of His",
        "smarts": "N[CH1]([CH2][cH0]1[nH0,nH1,nH1+][cH1][nH0,nH1,nH1+][cH1]1)C(=[OH0])",
        "fasta_pattern": "H",
        "factors_favoring_stability": "pH significantly lower than the apparent pKa of His",
        "factors_disfavoring_stability": "Transition metal ions, light, oxygen",
        "comments": "Low risk",
        "references": "Li, Schöneich, and Borchardt 1995; Sluyterman 1962"
    },
    {
        "liability": "Oxidation of Trp",
        "smarts": "N[CH1]([CH2][cH0]1[cH1][nH1][cH0]2[cH1][cH1][cH1][cH1][cH0]12)C(=[OH0])",
        "fasta_pattern": "W",
        "factors_favoring_stability": "Weak pH dependence",
        "factors_disfavoring_stability": "light, oxygen, transition metal ions",
        "comments": "Low risk",
        "references": "Li, Schöneich, and Borchardt 1995; Sluyterman 1962"
    },
    {
        "liability": "Oxidation of Tyr",
        "smarts": "N[CH1]([CH2][cH0]1[cH1][cH1][cH1]([OH1,OH0-])[cH1][cH1]1)C(=[OH0])",
        "fasta_pattern": "Y",
        "factors_favoring_stability": "pH significantly lower than apparent pKa of Tyr",
        "factors_disfavoring_stability": "Transition metal ions, light, oxygen",
        "comments": "High risk",
        "references": "Li, Schöneich, and Borchardt 1995"
    },
    {
        "liability": "Oxidation of Cys & thiols",
        "smarts": "C[CH2][SX2H1,S-X1H0]",
        "fasta_pattern": "C",
        "factors_favoring_stability": "pH significantly lower than apparent pKa of Cys",
        "factors_disfavoring_stability": (
            "Proximity to other Cys, lower thiol’s pKa, transition metal ions, light, oxygen"
        ),
        "comments": "High risk",
        "references": "Li, Schöneich, and Borchardt 1995"
    },
    {
        "liability": "Oxidation of Cys-Cys disulfide bonds",
        "smarts": "C[CH2][SX2][SX2][CH2]C",
        "fasta_pattern": "C",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": "Oxygen",
        "comments": "High risk",
        "references": "Karimi et al. 2016"
    },
    {
        "liability": "Oxidation of Met",
        "smarts": "C[CH2][SX2][CH3]",
        "fasta_pattern": "M",
        "factors_favoring_stability": "Weak/no pH dependence",
        "factors_disfavoring_stability": "Oxygen, solvent exposure & degree of hydration, metal ions",
        "comments": "High risk",
        "references": "Chu et al. 2004; Sluyterman 1962; Thakkar et al. 2023"
    },
    {
        "liability": "Cys-Cys disulfide bond scrambling & ß-elimination",
        "smarts": "C[CH2][SX2][SX2][CH2]C",
        "fasta_pattern": "C",
        "factors_favoring_stability": "-",
        "factors_disfavoring_stability": "Proximity to and accessibility by thiols (e.g. Cys), alkaline pH",
        "comments": "High risk",
        "references": "Kolšek, Aponte-Santamaría, and Gräter 2017; Sung et al. 2016"
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
            # correct for the case of symmetric smars patters
            if liability["liability"] == "Oxidation of thioethers": match_count = match_count / 2

            detected_liabilities[liability["liability"]] = {
                "liability": liability["liability"],
                "fasta_pattern": liability["fasta_pattern"],
                "match_count": match_count,
                "factors_favoring_stability": liability["factors_favoring_stability"],
                "factors_disfavoring_stability": liability["factors_disfavoring_stability"],
                "comments": liability["comments"],
                "references": liability["references"],
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
        
   #{
   #    "liability": "Cyclization of N-terminal Gln (pGlu formation)",
   #    "smarts": "[NH2,NH3+][CH1]([CH2][CH2][CH0](=[OH0])[NH2])C(=[OH0])",
   #    "fasta_pattern": "Q (at N-terminus)",
   #    "factors_favoring_stability": "pH 5.5 – 6.5",
   #    "factors_disfavoring_stability": "-",
   #    "comments": "-",
   #    "references": "Bersin, Patel, and Topp 2021; Chelius et al. 2006"
   #},
        # N-terminal case 1
        if liability["liability"] == "Cyclization of N-terminal Gln (pGlu formation)":
            if liability["fasta_pattern"].upper() == fasta_sequence_upper[0]:
                detected_liabilities[liability["liability"]] = {
                    "liability": liability["liability"],
                    "fasta_pattern": liability["fasta_pattern"],
                    "match_count": 1,
                    "factors_favoring_stability": liability["factors_favoring_stability"],
                    "factors_disfavoring_stability": liability["factors_disfavoring_stability"],
                    "comments": liability["comments"],
                    "references": liability["references"],
                }
   #{
   #    "liability": "Backbone hydrolysis between sequence positions 2 and 3 via diketopiperazine formation",
   #    "smarts": "[NH2,NH3+]C[CH0](=[OH0])[NH0]1[CH2][CH2][CH2][CH1]1[CH0](=[OH0])[NH1][CH2]C(=[OH0])",
   #    "fasta_pattern": "XPG (at N-terminus)",
   #    "factors_favoring_stability": "pH < ~5",
   #    "factors_disfavoring_stability": "Polar or charged amino acids in the 1st position, alkylammonium carboxylates",
   #    "comments": "Low risk",
   #    "references": "Capasso and Mazzarella 1999; Goolcharran and Borchardt 1998; Zhang, Hales, and Clemmer 2022"
   #},
        # N-terminal case 2
        elif liability["liability"] == "Backbone hydrolysis between sequence positions 2 and 3 via diketopiperazine formation":
            if liability["fasta_pattern"][1:2].upper() == fasta_sequence_upper[1:2]:
                detected_liabilities[liability["liability"]] = {
                    "liability": liability["liability"],
                    "fasta_pattern": liability["fasta_pattern"],
                    "match_count": 1,
                    "factors_favoring_stability": liability["factors_favoring_stability"],
                    "factors_disfavoring_stability": liability["factors_disfavoring_stability"],
                    "comments": liability["comments"],
                    "references": liability["references"],
                }           
        # skip the 2nd matching of the alternative tautomer liability entry to avoid duplicated entries for Arg
        elif liability["liability"] == "Oxidation (carbonylation) of Arg (alternative tautomer)":
                continue
        # general case
        elif liability["fasta_pattern"].upper() in fasta_sequence_upper:
                detected_liabilities[liability["liability"]] = {
                    "liability": liability["liability"],
                    "fasta_pattern": liability["fasta_pattern"],
                    "match_count": fasta_sequence.count(liability["fasta_pattern"]),
                    "factors_favoring_stability": liability["factors_favoring_stability"],
                    "factors_disfavoring_stability": liability["factors_disfavoring_stability"],
                    "comments": liability["comments"],
                    "references": liability["references"],
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
            dict_output[mol_idx] = {"mol_name": mol_name, "liabilities": liabilities}

        elif fasta:
            mol_name = mol_data.get("mol_name", "Unnamed Sequence")
            fasta = mol_data.get("fasta")

            if fasta is None:
                raise ValueError(f"FASTA sequence for {mol_name} is missing.")

            # Calculate liabilities based on FASTA sequence
            liabilities = calculate_liabilities_from_fasta(fasta)

            # Store the result in the output dictionary
            dict_output[mol_idx] = {"mol_name": mol_name, "liabilities": liabilities}
    return dict_output
