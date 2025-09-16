# Chemical Liability Detection

This repository provides a Python script for detecting chemical liabilities in molecules.  
The detection is implemented using **RDKit**.

---

## Requirements
- Python 3.8+
- [RDKit](https://www.rdkit.org/)

## Usage
```python
from rdkit import Chem
from liabilities import calculate_liabilities_from_mol

# Example peptide molecule (replace with your own RDKit molecule)
smiles_string = "CC(C)C[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@H](C)N)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](C)C(=O)NCC(=O)NCC(=O)O"
mol = Chem.MolFromSmiles(smiles_string)

# Calculate liabilities
matches = calculate_liabilities_from_mol(mol)

# Output the identified liabilities
if matches:
    for liability, details in matches.items():
        for k, v in details.items():
            print(f"{k}: {v}")
        print("***")
else:
    print("No liabilities detected in the sequence.")
>
liability: Racemization of Ala
fasta_pattern: A
match_count: 2
factors_favoring_stability: -
factors_disfavoring_stability: Disordered secondary structure & high conformational flexibility, solvent accessibility
comments: -
references: Cloos and Christgau 2002
***
liability: Oxidation of Trp
fasta_pattern: W
match_count: 1
factors_favoring_stability: Weak pH dependence
factors_disfavoring_stability: light, oxygen, transition metal ions
comments: Low risk
references: Li, Schöneich, and Borchardt 1995; Sluyterman 1962
***
```
