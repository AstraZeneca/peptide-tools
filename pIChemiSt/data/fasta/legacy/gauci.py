### Calibrated ExPASY - from Gauci et al. Proteomics 2008, 8, 4898 as implemented in pIR
SetName='Gauci'

# Acidic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
acidic1 = { 
 "D": [ 4.05, 4.05, 4.05 ],
 "C": [ 9.0 , 9.0 , 9.0  ],
 "E": [ 4.45, 4.45, 4.45 ],
 "Y": [ 10.0, 10.0, 10.0 ], 
 'U': [ 5.43, 5.20, 5.60 ]  # pK for U was taken from Byun et al. Biopolymers 2011, 95, 345
}


# Basic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
basic1 = { 
 "R": [ 12.0, 12.0, 12.0 ],
 "H": [ 5.98, 5.98, 5.98 ],
 "K": [ 10.0, 10.0, 10.0 ]
}


# Terminal_Amino_Acids
# AA N-term  C-Term
terminus_ionizable1 = { 
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


known_res=['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'X', 'Z', 'B', 'U']


def FillMissingAAtoterminus_ionizable(terminus_ionizable):
    
    # Calc average
    sumNterm=0
    sumCterm=0
    #for k,v in terminus_ionizable.iteritems():
    for k in terminus_ionizable.keys():
        v = terminus_ionizable[k]
        sumNterm += v[0]
        sumCterm += v[1]
    avNterm = sumNterm / len(terminus_ionizable.keys())
    avCterm = sumCterm / len(terminus_ionizable.keys())

    for R in known_res:
        if R not in terminus_ionizable.keys():
            if   R == 'X': terminus_ionizable[R] = [ avNterm, avCterm ]
            elif R == 'Z': terminus_ionizable[R] = [ (terminus_ionizable['E'][0]+terminus_ionizable['Q'][0])/2,  (terminus_ionizable['E'][1]+terminus_ionizable['Q'][1])/2  ]
            elif R == 'B': terminus_ionizable[R] = [ (terminus_ionizable['N'][0]+terminus_ionizable['D'][0])/2,  (terminus_ionizable['N'][1]+terminus_ionizable['D'][1])/2  ]
            elif R == 'U': 
                # copy of X
                terminus_ionizable[R] = [ avNterm, avCterm ]
            else:
                print("---!Error: data for specific -NH2 and -COOH termini pKa values for residue "+R+" is not given in the "+SetName+" pKa set. Set this residue identical to X (average of all available). Check set. Exit.")
                # sys.exit(1)
    
    return terminus_ionizable