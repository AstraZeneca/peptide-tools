###
### Last  update: Andrey Frolov, AstraZeneca, Molndal. 08/01/2020
### First verion: Andrey Frolov, AstraZeneca, Molndal. 11/02/2016
###

import sys
import optparse

#import numpy as np
#import pylab

#import json
#from json import encoder
#encoder.FLOAT_REPR = lambda o: format(o, '.2f')


######
###### Some global definitions
######

all_known_pKa_sets=['ProMoST',
'IPC_peptide',
'IPC2_peptide',
'Gauci',
'Bjellqvist',
'Rodwell',
'Grimsley',
'Thurlkill',
'EMBOSS',
'DTASelect',
'Solomon',
'Sillero',
'Lehninger',
'Toseland',
'Nozaki',
'Dawson']


def list_to_comma_seprated_string(l):
	s=""
	for v in l: s+=str(v)+","
	return s[:-1]

### Preselected set of pKa to display
#pKa_sets_to_use=['IPC_peptide','ProMoST','Gauci_calib','Bjellqvist','Rodwell','Grimsley','Thurlkill','Solomon','Lehninger','EMBOSS']
#pKa_sets_to_use=['IPC_peptide','ProMoST','Gauci','Bjellqvist','Grimsley','Thurlkill','Lehninger','Toseland']
pKa_sets_to_use=['IPC2_peptide','IPC_peptide','ProMoST','Gauci','Grimsley','Thurlkill','Lehninger','Toseland']

known_basic_res=['K','R','H']
known_acidic_res=['D','E','C','Y','U']
known_res=['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'X', 'Z', 'B', 'U']


def FillMissingAAtopKa_TerminusIonizableGroup(pKa_TerminusIonizableGroup):
    
    # Calc average
    sumNterm=0
    sumCterm=0
    #for k,v in pKa_TerminusIonizableGroup.iteritems():
    for k in pKa_TerminusIonizableGroup.keys():
        v = pKa_TerminusIonizableGroup[k]
        sumNterm += v[0]
        sumCterm += v[1]
    avNterm = sumNterm / len(pKa_TerminusIonizableGroup.keys())
    avCterm = sumCterm / len(pKa_TerminusIonizableGroup.keys())

    for R in known_res:
        if R not in pKa_TerminusIonizableGroup.keys():
            if   R == 'X': pKa_TerminusIonizableGroup[R] = [ avNterm, avCterm ]
            elif R == 'Z': pKa_TerminusIonizableGroup[R] = [ (pKa_TerminusIonizableGroup['E'][0]+pKa_TerminusIonizableGroup['Q'][0])/2,  (pKa_TerminusIonizableGroup['E'][1]+pKa_TerminusIonizableGroup['Q'][1])/2  ]
            elif R == 'B': pKa_TerminusIonizableGroup[R] = [ (pKa_TerminusIonizableGroup['N'][0]+pKa_TerminusIonizableGroup['D'][0])/2,  (pKa_TerminusIonizableGroup['N'][1]+pKa_TerminusIonizableGroup['D'][1])/2  ]
            elif R == 'U': 
                # copy of X
                pKa_TerminusIonizableGroup[R] = [ avNterm, avCterm ]
            else:
                print("---!Error: data for specific -NH2 and -COOH termini pKa values for residue "+R+" is not given in the "+SetName+" pKa set. Set this residue identical to X (average of all available). Check set. Exit.")
                sys.exit(1)
    
    return pKa_TerminusIonizableGroup
    





###------------------------------------------------------------------------------------------
### Sets of pKa 

pKa_sets={}
pKa_sets_short={}

pKa_sets_short['EMBOSS']={
 'K':       10.8 ,
 'R':       12.5 ,
 'H':       6.5  ,
 'D':       3.9  ,
 'E':       4.1  ,
 'C':       8.5  ,
 'Y':       10.1 ,
 'Nterm':   8.6  ,
 'Cterm':   3.6  
}


pKa_sets_short['IPC2_peptide']={
 'K':       8.165 ,
 'R':       11.493 ,
 'H':       6.439 ,
 'D':       3.969 ,
 'E':       4.507 ,
 'C':       9.439 ,
 'Y':       9.153 ,
 'Nterm':   7.947 ,
 'Cterm':   2.977  
}


pKa_sets_short['IPC_peptide']={
 'K':       10.517 ,
 'R':       12.503 ,
 'H':       6.018  ,
 'D':       3.887  ,
 'E':       4.317  ,
 'C':       8.297  ,
 'Y':       10.071 ,
 'Nterm':   9.564  ,
 'Cterm':   2.383  
}

#Amino acid	NH2	COOH	C	D	E	H	K	R	Y
pKa_sets_short['DTASelect']={
'Nterm':	8.0  ,
'Cterm':	3.1  ,
'C'    :	8.5  ,
'D'    :	4.4  ,
'E'    : 	4.4  ,
'H'    : 	6.5  ,
'K'    : 	10.0 ,
'R'    : 	12.0 ,
'Y'    : 	10.0
}

pKa_sets_short['Bjellqvist']={	
 'Nterm':       7.5   ,
 'Cterm':       3.55  ,
 'C'    :       9.0   ,
 'D'    :       4.05  ,
 'E'    :       4.45  ,
 'H'    :       5.98  ,
 'K'    :       10.0  ,
 'R'    :       12.0  ,
 'Y'    :       10.0  
}

pKa_sets_short['Solomon']={
 'Nterm':      9.6   ,
 'Cterm':      2.4   ,
 'C'    :      8.3   ,
 'D'    :      3.9   ,
 'E'    :      4.3   ,
 'H'    :      6.0   ,
 'K'    :      10.5  ,
 'R'    :     12.5   ,
 'Y'    :     10.1   
}

pKa_sets_short['Sillero']={
 'Nterm':      8.2  ,
 'Cterm':      3.2  ,
 'C'    :      9.0  ,
 'D'    :      4.0  ,
 'E'    :      4.5  ,
 'H'    :      6.4  ,
 'K'    :      10.4 ,
 'R'    :     12.0  ,
 'Y'    :     10.0  
}


pKa_sets_short['Rodwell']={
 'Nterm':      8.0   , 
 'Cterm':      3.1   ,
 'C'    :      8.33  ,
 'D'    :     3.68   ,
 'E'    :     4.25   ,
 'H'    :     6.0    ,
 'K'    :      11.5  ,
 'R'    :     11.5   ,
 'Y'    :     10.07  
}

pKa_sets_short['Lehninger']={
  'Nterm':     9.69  ,
  'Cterm':     2.34  ,
  'C'    :     8.33  ,
  'D'    :     3.86  ,
  'E'    :     4.25  ,
  'H'    :     6.0   ,
  'K'    :      10.5 ,
  'R'    :     12.4  ,
  'Y'    :     10.0  
}

pKa_sets_short['Grimsley']={
  'Nterm':     7.7   ,
  'Cterm':     3.3   ,
  'C'    :     6.8   ,
  'D'    :     3.5   ,
  'E'    :     4.2   ,
  'H'    :     6.6   ,
  'K'    :     10.5  ,
  'R'    :    12.04  ,
  'Y'    :    10.3   
}

pKa_sets_short['Toseland']={
   'Nterm':   8.71   ,
   'Cterm':   3.19   ,
   'C'    :   6.87   ,
   'D'    :   3.6    ,
   'E'    :   4.29   ,
   'H'    :   6.33   ,
   'K'    :   10.45  ,
   'R'    :   12.0   ,
   'Y'    :    9.61  
}


pKa_sets_short['Thurlkill']={
   'Nterm':   8.0   ,
   'Cterm':   3.67  ,
   'C'    :   8.55  ,
   'D'    :   3.67  ,
   'E'    :   4.25  ,
   'H'    :   6.54  ,
   'K'    :   10.4  ,
   'R'    :   12.0  ,
   'Y'    :   9.84  
}


pKa_sets_short['Nozaki']={
   'Nterm':     7.5   ,
   'Cterm':     3.8   ,
   'C'    :     9.5   ,
   'D'    :     4.0   ,
   'E'    :     4.4   ,
   'H'    :     6.3   ,
   'K'    :     10.4  ,
   'R'    :     12.0  ,
   'Y'    :     9.6   
}


pKa_sets_short['Dawson']={
    'Nterm':     8.2 ,
    'Cterm':     3.2 ,
    'C'    :     8.3 ,
    'D'    :     3.9 ,
    'E'    :     4.3 ,
    'H'    :     6   ,
    'K'    :    10.5 ,
    'R'    :    12.0 ,
    'Y'    :    10.0 
}


def ConvertpKaSetIntoProMoSTformat(pKaset):
    pKa_basic1={}
    pKa_acidic1={}
    pKa_TerminusIonizableGroup1={}

    for R in known_basic_res:
        if R in pKa_sets_short[pKaset].keys(): 
            pKa=pKa_sets_short[pKaset][R]
            pKa_basic1[R]=[pKa,pKa,pKa]
    
    for R in known_acidic_res:
        if R in pKa_sets_short[pKaset].keys(): 
            pKa=pKa_sets_short[pKaset][R]
            pKa_acidic1[R]=[pKa,pKa,pKa]
    
    for R in known_res:
            pKa_Cterm=pKa_sets_short[pKaset]['Cterm']
            pKa_Nterm=pKa_sets_short[pKaset]['Nterm']
            pKa_TerminusIonizableGroup1[R]=[pKa_Nterm,pKa_Cterm]

    pKa_sets[pKaset]={
     'pKa_acidic': pKa_acidic1,
     'pKa_basic': pKa_basic1,
     'pKa_TerminusIonizableGroup': pKa_TerminusIonizableGroup1
    }
    return



for pKaset in pKa_sets_short.keys():
#for pKaset in ['EMBOSS']:
    ConvertpKaSetIntoProMoSTformat(pKaset)


### Calibrated ExPASY - from Gauci et al. Proteomics 2008, 8, 4898 as implemented in pIR
SetName='Gauci'

# Acidic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
pKa_acidic1 = { 
 "D": [ 4.05, 4.05, 4.05 ],
 "C": [ 9.0 , 9.0 , 9.0  ],
 "E": [ 4.45, 4.45, 4.45 ],
 "Y": [ 10.0, 10.0, 10.0 ], 
 'U': [ 5.43, 5.20, 5.60 ]  # pK for U was taken from Byun et al. Biopolymers 2011, 95, 345
}


# Basic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
pKa_basic1 = { 
 "R": [ 12.0, 12.0, 12.0 ],
 "H": [ 5.98, 5.98, 5.98 ],
 "K": [ 10.0, 10.0, 10.0 ]
}


# Terminal_Amino_Acids
# AA N-term  C-Term
pKa_TerminusIonizableGroup1 = { 
 "A": [  7.59,   3.55 ], 
 "R": [  7.5,    3.55 ], 
 "N": [  6.7,    3.55 ], 
 "D": [  7.5,    4.55 ], 
 "C": [  6.5,    3.55 ], 
 "E": [  7.7,    4.75 ], 
 "Q": [  7.5,    3.55 ], 
 "G": [  7.5,    3.55 ], 
 "H": [  7.5,    3.55 ], 
 "I": [  7.5,    3.55 ], 
 "L": [  7.5,    3.55 ], 
 "K": [  7.5,    3.55 ], 
 "M": [  7.0,    3.55 ], 
 "F": [  7.5,    3.55 ], 
 "P": [  8.3599, 3.55 ],
 "S": [  6.93,   3.55 ], 
 "T": [  6.82,   3.55 ], 
 "W": [  7.5,    3.55 ], 
 "Y": [  7.5,    3.55 ], 
 "V": [  7.44,   3.55 ]   
}


pKa_sets[SetName]={
 'pKa_acidic': pKa_acidic1,
 'pKa_basic': pKa_basic1,
 'pKa_TerminusIonizableGroup': FillMissingAAtopKa_TerminusIonizableGroup(pKa_TerminusIonizableGroup1)
}



###
### Set ProMoST From http://proteomics.mcw.edu/promost_adv.html
###

# Acidic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
pKa_acidic1 = { 'D': [ 4.07,  3.57,  4.57 ],
                'E': [ 4.45,  4.15,  4.75 ],
                'C': [ 8.28,  8.00,  9.00 ],
                'Y': [ 9.84,  9.34, 10.34 ],
                'U': [ 5.43,  5.20,  5.60 ] } # pK for U was taken from Byun et al. Biopolymers 2011, 95, 345

# Basic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
pKa_basic1 = {'K':  [  9.8,  10.00,  10.30 ],
              'R':  [ 12.5,  11.50,  11.50 ],
              'H':  [ 6.08,   4.89,   6.89 ] }

# Terminal_Amino_Acids
# AA N-term  C-Term
pKa_TerminusIonizableGroup1 = { 
 'G': [ 7.50,  3.70 ],
 'A': [ 7.58,  3.75 ],
 'S': [ 6.86,  3.61 ],
 'P': [ 8.36,  3.40 ],
 'V': [ 7.44,  3.69 ],
 'T': [ 7.02,  3.57 ],
 'C': [ 8.12,  3.10 ],
 'I': [ 7.48,  3.72 ],
 'L': [ 7.46,  3.73 ],
 'N': [ 7.22,  3.64 ],
 'D': [ 7.70,  3.50 ],
 'Q': [ 6.73,  3.57 ],
 'K': [ 6.67,  3.40 ],
 'E': [ 7.19,  3.50 ],
 'M': [ 6.98,  3.68 ],
 'H': [ 7.18,  3.17 ],
 'F': [ 6.96,  3.98 ],
 'R': [ 6.76,  3.41 ],
 'Y': [ 6.83,  3.60 ],
 'W': [ 7.11,  3.78 ],
 'X': [ 7.26,  3.57 ],
 'U': [ 7.26,  3.57 ], ### copy of X
 'Z': [ 6.96,  3.54 ],
 'B': [ 7.46,  3.57 ]  }

pKa_sets['ProMoST']={
 'pKa_acidic': pKa_acidic1,
 'pKa_basic': pKa_basic1,
 'pKa_TerminusIonizableGroup': pKa_TerminusIonizableGroup1
}


### Noncanonical AAs. These values used for all sets of pKa available for standard AAs.  
# PTM. Not complete... to exdend upon request
pKa_noncanonical = {'pKa1_phosphate':1.2,
                    'pKa2_phosphate':6.9,
                    'dpKa_alkylLys': 0.15, # data from ACD lab: pKa of amine: 10.69. The delta for methylated amine compared to amine.  ### Zhang, Vogel, J. Bio. Chem. 1993, 268, 30, 22420 (Table III, Lys75) pKas of methylated 10.87, dimethylated 10.12, 
                    'dpKa_dialkylLys': 0.15 - 0.75
                   } # data from Zhang, Vogel et al.  (ACD lab: pKa of dimethylamine: 9.83 +- 0.28 - error too high. The delta for methylated amine compared to amine.  ### Zhang, Vogel, J. Bio. Chem. 1993, 268, 30, 22420 (Table III, Lys75) pKas of methylated 10.87, dimethylated 10.12, 





#Ala	A	Alanine
#Arg	R	Arginine
#Asn	N	Asparagine
#Asp	D	Aspartic acid
#Cys	C	Cysteine
#Gln	Q	Glutamine
#Glu	E	Glutamic acid
#Gly	G	Glycine
#His	H	Histidine
#Ile	I	Isoleucine
#Leu	L	Leucine
#Lys	K	Lysine
#Met	M	Methionine
#Phe	F	Phenylalanine
#Pro	P	Proline
#Pyl	O	Pyrrolysine
#Ser	S	Serine
#Sec	U	Selenocysteine
#Thr	T	Threonine
#Trp	W	Tryptophan
#Tyr	Y	Tyrosine
#Val	V	Valine
#Asx	B	Aspartic acid or Asparagine
#Glx	Z	Glutamic acid or Glutamine
#Xaa	X	Any amino acid
#Xle	J	Leucine or Isoleucine
#TERM		termination codon




