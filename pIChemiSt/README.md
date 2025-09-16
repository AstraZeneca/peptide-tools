# pIChemiSt

## Description
The program calculates the isoelectric point of proteins or peptides based on their 2D molecular structure. The input structure is cut into monomers by targeting its amide bonds, and then each monomer's pKa value is determined using different methods: natural amino acids pKa values are matched against a dictionary; non-natural amino acid values are calculated using either pKaMatcher (built-in tool based on SMARTS patterns) or ACD perceptabat GALAS algorithm (a commercial tool that requires licence). For natural amino acids the following predefined sets of amino-acid pKa values are implemented: 'IPC2_peptide', 'IPC_peptide', 'ProMoST', 'Gauci', 'Rodwell', 'Grimsley', 'Thurlkill', 'Solomon', 'Lehninger', 'EMBOSS'. The mean value and variation between different sets are also calculated as well as the total charge at pH 7.4. The program can also plot the corresponding pH/Q curves for each input structure. Please refer to pIChemiSt publication for more details: https://pubs.acs.org/doi/10.1021/acs.jcim.2c01261

## How to install the software via pypi
- Ensure that you have Python version >=3.9
- Run `pip install pichemist`
- (optional) - To use ACD for the prediction of non-natural amino acid pKa, make sure that the command `perceptabat` points to its binary

### How to install the software via Github
- Clone the repository
- Ensure that you have Python version >=3.9
- Enter the package folder `cd peptide-tools/pIChemiSt`
- Run `pip install .` to install the Python library and the CLI command
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
# Note that FASTA assumes that C- and N- termini are ionisable by default
pichemist -i "MNSERSDVTLY" -if fasta_stdin

# Use FASTA with capped termini, i.e., non-ionisable, for both C- and N-.
# These can be configured as preferred by removing the corresponding flags.
# This configuration can also be used as a 'trick' for feeding cyclic peptides
# as linear FASTA sequences as their termini will be assumed to be non ionisable.
pichemist -i "MNSERSDVTLY" -if fasta_stdin --ionizable_nterm false --ionizable_cterm false

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
import pprint
from pichemist.io import generate_input
from pichemist.api import pichemist_from_dict

smiles = "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O"

args = {
        "input_data": smiles,
        "input_format": "smiles_stdin",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "pkamatcher",
    }

input_dict = generate_input(args["input_format"], args["input_data"])
output = pichemist_from_dict(
    input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
)

pp = pprint.PrettyPrinter(depth=4)
pp.pprint(output)
>
{1: {'QpH7': {'Gauci': 0.5541,
              'Grimsley': 0.6645,
              'IPC2_peptide': 0.6315,
              'IPC_peptide': 0.9916,
              'Lehninger': 0.9932,
              'ProMoST': 0.2617,
              'Q at pH7.4 mean': 0.7307,
              'Thurlkill': 0.7975,
              'Toseland': 0.9516,
              'err': 0.2386,
              'std': 0.675},
     'mol_name': 'C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O',
     'pI': {'Gauci': 8.6875,
            'Grimsley': 8.9375,
            'IPC2_peptide': 8.046875,
            'IPC_peptide': 9.8125,
            'Lehninger': 9.859375,
            'ProMoST': 8.375,
            'Thurlkill': 9.0625,
            'Toseland': 9.40625,
            'err': 0.6087,
            'pI mean': 9.0234,
            'std': 1.7216},
     'pI_interval': (8.625, 9.3625),
     'pI_interval_threshold': 0.2,
     'pKa_set': 'IPC2_peptide'}}
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
- We strongly recommend using `pre-commit` when contributing to this repo. The root folder of peptide-tools contains a `.pre-commit-config.yaml` which can be used to set a `pre-commit` hook and automatically run a series of validations. 
