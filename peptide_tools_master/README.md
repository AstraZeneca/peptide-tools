# peptide_tools_master

## Description
Wrapper CLI of all peptide tools. TODO... The interface accepts different types of inputs including SMILES, SMILES files, SDF, FASTA, FASTA sequences, and FASTA files. Some logic for the recognition of the input type is implemented. This was done to facilitate the integrate the CLI with a front end, removing the need for users to specify their input types.

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
> 
```
