# pIChemiSt

## Description
The program calculates the isoelectric point of proteins or peptides based on their 2D structure. The input structure is cut into monomers by targeting its amide bonds, and then each monomer's pKa value is determined using different methods: Natural amino acids pKa values are matched against a dictionary; non-natural amino acid values are calculated using either pKaMatcher (built-in tool based on SMARTS patterns) or ACD perceptabat GALAS algorithm (a commercial tool that requires licence). For natural amino acids the following sets of amino-acid pKa values are implemented: 'IPC2_peptide', 'IPC_peptide', 'ProMoST', 'Gauci', 'Rodwell', 'Grimsley', 'Thurlkill', 'Solomon', 'Lehninger', 'EMBOSS' as described in http://isoelectric.org. The mean value and variation between different sets are also calculated as well as the total charge at pH 7.4. The program can also plot the corresponding pH/Q curves for each input structure.

## How to install the software
- Clone the repository
- Ensure that you have Python version >=3.8
- Enter the package folder `cd peptide-tools/pIChemiSt`
- Run `pip install .` to install the Python library
- Run `sh setup.cli` to configure the CLI (effects take change only when a new terminal is started)
- (optional) - To use ACD for the prediction of non-natural amino acid pKa, make sure that the command `perceptabat` points to its binary

## Examples of usage (CLI)
```bash
# Run the predictor against a SMILES file using pKaMatcher and output results to console
pichemist -i test/examples/payload_1.smi --method pkamatcher

>
======================================================================================================================================================
pI
---------------------------------
     pI mean  9.02
         err  0.61
         std  1.72
IPC2_peptide  8.05
 IPC_peptide  9.81
     ProMoST  8.38
       Gauci  8.69
    Grimsley  8.94
   Thurlkill  9.06
   Lehninger  9.86
    Toseland  9.41


======================================================================================================================================================
Q at pH7.4
---------------------------------
Q at pH7.4 mean  0.73
         err  0.24
         std  0.67
IPC2_peptide  0.63
 IPC_peptide  0.99
     ProMoST  0.26
       Gauci  0.55
    Grimsley  0.66
   Thurlkill  0.8
   Lehninger  0.99
    Toseland  0.95


pH interval with charge between -0.2 and  0.2 and prediction tool: pkamatcher
pI interval:  8.6 -  9.4
```

Other flavours of flags and arguments can be configured to feed different inputs or produce different outputs:
```bash
# Use SMILES string as input
pichemist -i "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=O)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=O)N)C(=O)O" -if smiles_stdin

# Use SMILES string as input and output JSON to console
pichemist -i "C([C@@H](C(=O)O)N)SSC[C@@H](C(=O)O)N" -if smiles_stdin -of json

# Use FASTA as input
# Note that FASTA cannot be used to indicate if the N- and C- termini are capped or not.
# Most of natural peptides have an acid on the C-terminus, however, synthetic peptides
# may have an amide (not ionizable) at the C-terminus.
pichemist -i "MNSERSDVTLY" -if fasta_stdin

# Use SDF as input
pichemist -i test/examples/payload_4.sdf -if sdf

# Output SDF
pichemist -i test/examples/payload_1.smi -o results.sdf -of sdf

# Output JSON to console
pichemist -i test/examples/payload_1.smi -of json

# Plot pH/Q curve
pichemist -i test/examples/payload_1.smi --plot_ph_q_curve

# Plot pH/Q curve (with custom prefix 'plot')
pichemist -i test/examples/payload_2.smi --plot_ph_q_curve -pp "plot"

# Print the pKa values of fragments
pichemist -i test/examples/payload_1.smi --print_fragment_pkas

# Use ACD instead of pKaMatcher
pichemist -i test/examples/payload_1.smi --method acd

# Use ACD with SMILES string
pichemist -i "NCCC(=O)N[C@@H](Cc1c[nH]cn1)C(=O)O" --method acd -if smiles_stdin
```

## Examples of usage (Python API)
```python
from pichemist.io import generate_input
from pichemist.api import pichemist_from_list

smiles = "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O"

args = {
        "input_data": smiles,
        "input_format": "smiles_stdin",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "pkamatcher",
    }

input_dict = generate_input(args["input_format"], args["input_data"])
output = pichemist_from_list(
    input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
)

print(output)
>
{1: {'mol_name': 'C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O', 'pI': {'IPC2_peptide': 8.046875, 'IPC_peptide': 9.8125, 'ProMoST': 8.375, 'Gauci': 8.6875, 'Grimsley': 8.9375, 'Thurlkill': 9.0625, 'Lehninger': 9.859375, 'Toseland': 9.40625, 'pI mean': 9.0234375, 'std': 1.721588565104915, 'err': 0.6086734743994516}, 'QpH7': {'IPC2_peptide': 0.6314906212267486, 'IPC_peptide': 0.9915539516610472, 'ProMoST': 0.26174063515548607, 'Gauci': 0.5540630760817584, 'Grimsley': 0.6645409545014482, 'Thurlkill': 0.797542965316429, 'Lehninger': 0.9932283675959863, 'Toseland': 0.9515959465104951, 'Q at pH7.4 mean': 0.7307195647561748, 'std': 0.6749606913955383, 'err': 0.23863464096007284}, 'pI_interval': (8.624999999999998, 9.362499999999997), 'pI_interval_threshold': 0.2, 'pKa_set': 'IPC2_peptide'}}
```

## Contributions
- Andrey I. Frolov (https://orcid.org/0000-0001-9801-3253)
- Gian Marco Ghiandoni (https://orcid.org/0000-0002-2592-2939)
- Jonas Boström - contributed the function of splitting molecules by peptide bonds (https://orcid.org/0000-0002-9719-9137)
- Johan Ulander - contributed the idea to do SMARTS pattern matching Smiles into Amino Acid single letter code (https://orcid.org/0009-0004-7655-2212)

## For developers
- To create a new build, the package version first needs to be configured inside `setup.py` then the command `python setup.py sdist` will build the distribution. The command `bdist_wheel` should not be used since this mode in `setup.py` skips including the required static files in the wheel distribution
- The code can be automatically tested using `python setup.py test` which requires `pytest` to be installed
- Tests can also be run using the `Makefile` in the root of the repository. The file allows granular testing as follows:
  - `make test_core` runs only the `core` tests including pKaMatcher and plots
  - `make test_acd` only runs `acd` tests (which require an ACD license)
  - `make test` runs both `core` and `acd` tests
- Before committing new code, please always check that the formatting is consistent using `flake8`
