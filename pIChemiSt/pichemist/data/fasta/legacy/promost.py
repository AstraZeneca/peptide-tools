
###
### Set ProMoST From http://proteomics.mcw.edu/promost_adv.html
###
SetName='ProMoST'

# Acidic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
acidic = { 'D': [ 4.07,  3.57,  4.57 ],
                'E': [ 4.45,  4.15,  4.75 ],
                'C': [ 8.28,  8.00,  9.00 ],
                'Y': [ 9.84,  9.34, 10.34 ],
                'U': [ 5.43,  5.20,  5.60 ] } # pK for U was taken from Byun et al. Biopolymers 2011, 95, 345

# Basic_Amino_Acids
#             AA    Primary  N-Terminal  C-Terminal
basic = {'K':  [  9.8,  10.00,  10.30 ],
              'R':  [ 12.5,  11.50,  11.50 ],
              'H':  [ 6.08,   4.89,   6.89 ] }

# Terminal_Amino_Acids
# AA N-term  C-Term
terminus_ionizable = { 
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

# PKA_SETS[SetName]={
#  'acidic': acidic,
#  'basic': basic,
#  'terminus_ionizable': terminus_ionizable
# }
