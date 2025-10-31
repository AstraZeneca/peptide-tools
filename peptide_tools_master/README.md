# Peptide Tools Master

## Description
CLI layer for accessing multiple tools with a single command. The tools include pIChemiSt, the prediction of the extinction coefficient, the identification of chemical liabilities, and the calculation of molecular descriptors. The interface accepts different types of inputs including SMILES, SMILES files, SDF, FASTA, FASTA sequences (including D amino acids), and FASTA files. Some logic for the recognition of the input type is implemented. This was done to facilitate the integrate the CLI with a front end, removing the need for users to specify their input types.

## How to set up a temporary environment to run the tool
- Ensure that you have Python version >=3.8
- Clone the repository and access the folder `peptide_tools_master`
```bash
pip install pichemist
> Successfully installed pichemist-0.2.0

pwd
> /../peptide-tools/peptide_tools_master
PEPTIDE_TOOLS_PATH=`cd .. && echo $PWD`
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/smi2scrambledfasta
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/extn_coeff_fasta
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/pIChemiSt
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/liabilities
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/molecular_descriptors

python peptide_tools_master.py --input FPYVAE
> {
  "output_descriptors": {
    "1": {
      "mol_name": "none",
      "molecular_weight": 724.8,
      "seq_length": 6,
      "logp": 0.06,
      ...
}
```

Other examples of commands as provided as follows.

```bash
python peptide_tools_master.py --input 'C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N'

python peptide_tools_master.py --input 'C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N' --print_fragment_pkas

python peptide_tools_master.py --input '>sp|P43220|GLP1R_HUMAN Glucagon-like peptide 1 receptor OS=Homo sapiens OX=9606 GN=GLP1R PE=1 SV=2
MAGAPGPLRLALLLLGMVGRAGPRPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDL
FCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPW
RDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIH
LNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAAN
YYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRN
SNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLL
GTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERW
RLEHLHIQRDSSMKPLKCPTSSLSSGATAGSSMYTATCQASCS'

python peptide_tools_master.py --input FPYKPAE --print_fragment_pkas --ionizable_nterm false --ionizable_cterm false
```

## Contributions
- Andrey I. Frolov (https://orcid.org/0000-0001-9801-3253)
- Gian Marco Ghiandoni (https://orcid.org/0000-0002-2592-2939)

## For developers
- The code can be automatically tested using `python setup.py test` which requires `pytest` to be installed
- Tests can also be run using the `Makefile` in the root of the repository
- We strongly recommend using `pre-commit` when contributing to this repo. The root folder of peptide-tools contains a `.pre-commit-config.yaml` which can be used to set a `pre-commit` hook and automatically run a series of validations. 
