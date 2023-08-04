# pIChemiSt

## Description
The program calculates the isoelectric point of proteins or peptides based on their 2D structure. The input structure is cut into monomers by targeting its amide bonds, and then each monomer's pKa value is determined using different methods: Natural amino acids pKa values are matched against a dictionary; non-natural amino acid values are calculated using either pKaMatcher (built-in tool based on SMARTS patterns) or ACD perceptabatch (a commercial tool that requires licence). For natural amino acids the following sets of amino-acid pKa values are implemented: 'IPC2_peptide', 'IPC_peptide', 'ProMoST', 'Gauci', 'Rodwell', 'Grimsley', 'Thurlkill', 'Solomon', 'Lehninger', 'EMBOSS' as described in http://isoelectric.org. The mean value and variation between different sets are also calculated as well as the total charge at pH 7.4. The program can also plot the corresponding titration curves for each input structure.

## How to install the software
- Clone the repository
- Ensure that you have Python version >=3.8
- Enter the package folder `cd peptide-tools/pIChemiSt`
- Run `pip install dist/pichemist-*.whl`
- (optional) - To use ACD for the prediction of non-natural amino acid pKa, make sure that the command `perceptabat` points to its binary

## Examples of usage
### TODO: Examples divided per sections
```bash
python3 pichemist/cli.py -i "C([C@@H](C(=O)O)N)SSC[C@@H](C(=O)O)N" -if smiles_stdin --plot_titration_curve --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/examples/payload_1.smi --plot_titration_curve --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/examples/payload_1.smi -if smiles_file --plot_titration_curve --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/examples/payload_2.smi --plot_titration_curve -tfp plot --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/validation_set_modified_peptides.smi -if smiles_file --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/validation_set_modified_peptides.sdf -if sdf --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/examples/payload_1.smi -of json --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/examples/payload_1.smi --print_fragment_pkas --method acd
python3 pichemist/cli.py -i test/examples/payload_1.smi -o results.sdf -of sdf --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/examples/payload_2.smi -of json --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i test/examples/payload_1.smi --print_fragment_pkas --method pkamatcher
python3 pichemist/cli.py -i "CC(=O)NCC(=O)NC" -if smiles --print_fragment_pkas --method pkamatcher
```

## Contributions
Andrey I. Frolov (https://orcid.org/0000-0001-9801-3253)
Gian Marco Ghiandoni (https://orcid.org/0000-0002-2592-2939)
Jonas Boström - contributed function of splitting molecule by peptide bonds (https://orcid.org/0000-0002-9719-9137)
Johan Ulander - contributed idea to do SMARTS pattern matching Smiles into Amino Acid single letter code (https://orcid.org/0009-0004-7655-2212)
