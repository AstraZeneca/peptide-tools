![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

peptide_tools_mater.py

Wrapper script to run all peptide tools programs in one go. 


HOW TO RUN


    ### Set path to dependencies
    #---------------------------

    ml rdkit
    ml Biopython

    PEPTIDE_TOOLS_PATH=`echo $PWD`

    export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/smi2scrambledfasta
    export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/extn_coeff_fasta
    export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/pI_fasta
    export PYTHONPATH=${PYTHONPATH}:${PEPTIDE_TOOLS_PATH}/pIChemiSt

    ### Example peptide structure
    #----------------------------
    INPUT_SMILES='C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N'    # FPYVAE - peptide



    ### Example usage
    #----------------
    cd TEST
    python peptide_tools_master.py --input "${INPUT_SMILES}" # SMILES input
    python peptide_tools_master.py --input  FPYVAE           # Fasta sequence input
    python peptide_tools_master.py --input  myfile.smi       # Smiles file input 
    python peptide_tools_master.py --input  myfile.fasta     # FASTA file input 
    python peptide_tools_master.py --input  myfile.sdf       # SDF file input 


DEPENDENCIES 

    python3 or later 
    rdkit/2021.03.1 (tested with)
    matplotlib/3.0.3 (tested with) 
    Biopython/1.73 (tested with)
    
    if --use_acdlabs enabled:
        acdperceptabatch     # "perceptabat" executable should be callable

    if --use_pkamatcher enabled:
        uses internal pKaMatcher to assign missing pKas


PLATFORM

    Tested on linux CentOS

