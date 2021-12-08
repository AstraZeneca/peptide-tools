![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

pI_fasta.py

Program calculates isoelectic point of protein/peptide based on the FASTA sequence. The following sets of amino-acid pKa values are implemented:'IPC_peptide','ProMoST','Gauci_calib','Bjellqvist','Rodwell','Grimsley','Thurlkill','Solomon','Lehninger','EMBOSS' as described in http://isoelectric.org 
The mean value and variation between different sets are also calculated. The program can plot the corresponding titration curves. Also the total charge at pH 7.4 is reported.


HOW TO RUN

    module load matplotlib
    python3 ../pI_fasta.py -i P43220.fasta -x

DEPENDENCIES 

    python3 or later 
    matplotlib/3.0.3 (tested with) 
    Biopython/1.73 (tested with)

PLATFORM

    Tested on linux CentOS

