
### Set path to dependencies
#---------------------------

ml rdkit
ml Biopython

PEPTIDE_TOOLS_PATH=`echo $PWD`

export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/smi2scrambledfasta
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/extn_coeff_fasta
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/pI_fasta
export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/rdkit_pI

### Example peptide structure
#----------------------------
INPUT_SMILES='C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N'    # FPYVAE - peptide

### Example usage
#----------------
#python peptide_tools_master.py --input "${INPUT_SMILES}"


