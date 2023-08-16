![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

smi2scrambledfasta.py

Helper script to generate a scrambled fasta from smiles. Useful for calculating phys-chem properties of peptides which are not very sensitive to the amino acid order in teh sequence. E.g. sequence-based calculation o fthe extinction coefficients and sequence-based calculation of the isoelectric point.


HOW TO RUN

    module load rdkit
    export PYTHONPATH=<your_path>/peptide_tools/pIChemiSt_v3.1:${PYTHONPATH}
    cd TEST
    python ../smi2scrambledfasta.py -i FPYVAE.smi -o FPYVAE.fasta

DEPENDENCIES 

    python2 or later 
    rdkit/2020.03.1 (tested with)
    pIChemiSt v3.1 or later

PLATFORM

    Tested on linux CentOS

