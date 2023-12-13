![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

peptide-tools-refactor

OVERVIEW
    Set of programs to calcaulte phys-chem properties of synthetic peptides and proteins: isoelectic point and extinction coefficients.
    In total there are 5 scripts/programs there, each can serve as stand alone:

    1) pIChemiSt  
    program to calcualte isoelectric point of natural and modified peptides. 

    2) smi2scrambledfasta  (depends on 1)
    program to generate the scrambled FASTA sequence from Smiles.

    3) pI_fasta  (depends on 2)
    program to calcualte isoelecric point using FASTA sequence as an input.

    4) extn_coeff_fasta  (depends on 2)
    program to calculate extinction coefficients.

    5) peptide_tools_master  (depends on 1,2,3,4)
    wrapper script for all the programs above. Can run for one compound and for batch submissions.


HOW TO RUN & DEPENDENCIES 

    Check the documentaion for each script. 

PLATFORM

    Tested on linux CentOS

