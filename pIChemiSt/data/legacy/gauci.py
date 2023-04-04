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