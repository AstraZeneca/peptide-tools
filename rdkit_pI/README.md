![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

rdkit_pI.py

Program calculates isoelectic point of protein/peptide based on 2D structure (SMILES string). For nonatural amino acids pKa values are predicted on the fly either using pKaMatcher (built-in tool based on SMARTS patters) orACDlabs perceptabatch pKa GALAS method (commercial tool that requires licence). For natural amino acids the following sets of amino-acid pKa values are implemented:'IPC_peptide','ProMoST','Gauci_calib','Bjellqvist','Rodwell','Grimsley','Thurlkill','Solomon','Lehninger','EMBOSS' as described in http://isoelectric.org. The mean value and variation between different sets are also calculated. The program can plot the corresponding titration curves. Also the total charge at pH 7.4 is reported. 


HOW TO RUN

    module load matplotlib
    module load rdkit

    # Use internal pKaMatcher for pKa calculation of nonatural amino acids
    python3 ../rdkit_pI.py -i Phe-Ornithine-aMeAsp-Lys-dAla.smi --plot_titration_curve --print_fragment_pkas --use_pkamatcher

    or

    # Use ACDlabs for pKa calculation of nonatural amino acids
    module load acdperceptabatch     # "perceptabat" executable should be callable
    python3 ../rdkit_pI.py -i Phe-Ornithine-aMeAsp-Lys-dAla.smi --plot_titration_curve --print_fragment_pkas --use_acdlabs
    
    

DEPENDENCIES 

    python3 or later 
    rdkit/2021.03.1 (tested with)
    matplotlib/3.0.3 (tested with) 
    Biopython/1.73 (tested with)
    
    if --use_acdlabs enabled:
        acdperceptabatch     # "perceptabat" executable should be callable


PLATFORM

    Tested on linux CentOS

