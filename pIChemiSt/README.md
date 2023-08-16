# pIChemiSt

## Description
The program calculates the isoelectric point of proteins or peptides based on their 2D structure. The input structure is cut into monomers by targeting its amide bonds, and then each monomer's pKa value is determined using different methods: Natural amino acids pKa values are matched against a dictionary; non-natural amino acid values are calculated using either pKaMatcher (built-in tool based on SMARTS patterns) or ACD perceptabat GALAS algorithm (a commercial tool that requires licence). For natural amino acids the following sets of amino-acid pKa values are implemented: 'IPC2_peptide', 'IPC_peptide', 'ProMoST', 'Gauci', 'Rodwell', 'Grimsley', 'Thurlkill', 'Solomon', 'Lehninger', 'EMBOSS' as described in http://isoelectric.org. The mean value and variation between different sets are also calculated as well as the total charge at pH 7.4. The program can also plot the corresponding pH/Q curves for each input structure.

## How to install the software
- Clone the repository
- Ensure that you have Python version >=3.8
- Enter the package folder `cd peptide-tools/pIChemiSt`
- Run `pip install dist/pichemist-*` to install the Python library
- Run `sh setup.cli` to configure the CLI (effects take change only when the a new terminal is started)
- (optional) - To use ACD for the prediction of non-natural amino acid pKa, make sure that the command `perceptabat` points to its binary

## Examples of usage
```bash
# Run the predictor against a SMILES file using pKaMatcher and output results to console
pichemist -i test/examples/payload_1.smi --method pkamatcher

>>>
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
pichemist -i "C([C@@H](C(=O)O)N)SSC[C@@H](C(=O)O)N" -if smiles_stdin

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
```

## Contributions
- Andrey I. Frolov (https://orcid.org/0000-0001-9801-3253)
- Gian Marco Ghiandoni (https://orcid.org/0000-0002-2592-2939)
- Jonas Boström - contributed function of splitting molecule by peptide bonds (https://orcid.org/0000-0002-9719-9137)
- Johan Ulander - contributed idea to do SMARTS pattern matching Smiles into Amino Acid single letter code (https://orcid.org/0009-0004-7655-2212)

## For developers
- The package can be installed from the wheel in the `dist/` folder. When a new version needs to be released, a new build must be created. That can be done by changing the version of the package inside `setup.py` then calling `python setup.py sdist` which will build the release. TODO: Do not use wheels 
- The code can be automatically tested using `python setup.py test` which requires `pytest` to be installed
- Tests can also be run using the `Makefile` in the root of the repository. The file allows granular testing as follows:
  - `make test_core` runs only the `core` tests including pKaMatcher and plots
  - `make test_acd` only runs `acd` tests (which require an ACD license)
  - `make test` runs both `core` and `acd` tests
- Before committing new code, please always check that the formatting is consistent using `flake8`