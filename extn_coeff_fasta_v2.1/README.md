![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

extn_coeff_fasta.py

Script calculates extinction (molar absorption) coefficients for proteins/peptides based on the FASTA sequence at three different wavelength: 205, 214, 280 nm. The methods are described in the following publications:

    Pace, Vajdos, Fee, Grimsley "How to measure and predict the molar absorption coefficient of a protein", Protein Science 1995, 4, 2411-2423
    Kuipers, Gruppen, "Prediction of molar extinction coefficients of proteins and peptides ...", J. Agric. Food Chem. 2007, 55, 5445
    Anthis, Clore, "Sequence-specific determination of protein and peptide concentrations by absorbance at 2015 nm", Protein Science 2013, 22, 851


HOW TO RUN

    python3 extn_coeff_fasta.py -i TEST/P43220.fasta 

DEPENDENCIES 

    python3 or later 
    Biopython/1.73 (tested with)

PLATFORM

    Tested on linux CentOS

