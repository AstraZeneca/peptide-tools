![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

peptide-tools

OVERVIEW
    Set of programs to calcaulte phys-chem properties of synthetic peptides and proteins: isoelectic point and extinction coefficients.
    In total there are 5 scripts/programs there, each can serve as stand alone:

    1) rdkit_pI_v3.2  (depends on 5)
    program to calcualte isoelectric point of natural and modified peptides. 

    2) smi2scrambledfasta_v1.0  (depends on 1)
    program to generate the scrambled FASTA sequence from Smiles.

    3) pI_fasta_v1.4  (depends on 2)
    program to calcualte isoelecric point using FASTA sequence as an input.

    4) extn_coeff_fasta_v2.1  (depends on 2)
    program to calculate extinction coefficients.

    5) dimorphite_dl-1.2.4_for_rdkit_pI  
    slightly modified version of the Dimorphite_DL-1.2.4 code, avaialble on GitHub, to enable pKa predition of unknown nonatural amino-acids side-chains.  


HOW TO RUN & DEPENDENCIES 

    Check the documentaion for each script. 

PLATFORM

    Tested on linux CentOS

