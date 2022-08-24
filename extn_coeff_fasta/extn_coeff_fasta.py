#!/usr/bin/env python
###
### Last  update: Andrey Frolov, AstraZeneca, Molndal. 29/01/2021
### First verion: Andrey Frolov, AstraZeneca, Molndal. 11/02/2016
###
import sys, os
import argparse
import json
from copy import copy
#from json import encoder
#encoder.FLOAT_REPR = lambda o: format(o, '.2f')


######
###### Some global definitions
######

known_basic_res=['K','R','H']
known_acidic_res=['D','E','C','Y','U']
known_res=['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'X', 'Z', 'B', 'U']

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

#####======================================================================================================================================================
##### START
##### Calc Molar Absorption coefficient

# Turns a dictionary into a class 
class Dict2Class(object): 
    def __init__(self, my_dict): 
        for key in my_dict: 
            setattr(self, key, my_dict[key]) 








def calc_extn_coeff(options={}):

    args = Dict2Class(options)


    # Get options
    if len(args.seq)!=0:
        # assume single fasta input
        mol_unique_ind = 1
        mol_name='unknown'
        fasta = args.seq
        mol_supply_json={}
        mol_supply_json[mol_unique_ind] = {'mol_name': mol_name, 'mol_obj':None, 'fasta':fasta}

    elif len(args.inputFile)!=0:
        # Assume filename as input
        inputFile = args.inputFile
        mol_supply_json = read_fasta_file(inputFile)

    elif len(args.inputJSON)!=0:
        # Assume molecule JSON supply as input
        mol_supply_json = json.loads(args.inputJSON)
    elif args.inputDict: # if not an empty dictionary
        # Assume Dict molecule supply as input
        mol_supply_json = args.inputDict
    else:
        raise Exception('Error: either fasta, input file *.fasta or JSON should be given. Exit. ')
        sys.exit(1)



    dict_out = {}
    #molid_list = []
    #molid_ind_list = []
    #molid_ind = 0
    #for molid,molfasta in suppl:
    for mol_unique_ind in mol_supply_json.keys():
        #molid_ind += 1
        options_single = copy(options)
        #options_single["seq"] = molfasta
        #dict_out_single = calc_extn_coeff_single_sequence(options_single)
        dict_out_single = calc_extn_coeff_single_sequence(mol_supply_json[mol_unique_ind]['fasta'])
        dict_out_single['mol_name'] = mol_supply_json[mol_unique_ind]['mol_name']
        dict_out[mol_unique_ind] = dict_out_single
        #molid_list.append(molid)
        #molid_ind_list += [molid_ind]

    #dict_out['molid_ind_list'] = molid_ind_list

    return dict_out







def calc_extn_coeff_single_sequence(sequence):

    for R in sequence:
        if R not in known_res: 
            raise Exception("---!Error: residue "+R+" is not known. Please use X if this is a noncaninical residue. Exiting.")
            sys.exit(1)
 

    ### Pace, Vajdos, Fee, Grimsley "How to measure and predict the molar absorption coefficient of a protein", Protein Sciecnce 1995, 4, 2411-2423

    nW = sequence.count('W')
    nF = sequence.count('F')
    nY = sequence.count('Y')
    nH = sequence.count('H')
    nM = sequence.count('M')
    nR = sequence.count('R')
    nC = sequence.count('C')
    nN = sequence.count('N')
    nQ = sequence.count('Q')
    nA = sequence.count('A')
    nD = sequence.count('D')
    nE = sequence.count('E')
    nG = sequence.count('G')
    nI = sequence.count('I')
    nL = sequence.count('L')
    nK = sequence.count('K')
    nP = sequence.count('P')
    nS = sequence.count('S')
    nT = sequence.count('T')
    nV = sequence.count('V')
    nPepBond = len(sequence) -1


    e205 =( 20400*nW + 
    8600*nF +
    6080*nY +
    5200*nH +
    1830*nM +
    1350*nR +
    690*nC +
    400*nN +
    400*nQ +
    2780*nPepBond )


    e214 =( 29050*nW +
    5200*nF +
    5375*nY +
    5125*nH +
    980*nM +
    102*nR +
    225*nC +
    136*nN +
    142*nQ +
    32*nA +
    58*nD +
    78*nE +
    21*nG +
    45*nI +
    45*nL +
    41*nK +
    2675*sequence[1:].count('P') + 30*sequence[0].count('P') +
    34*nS +
    41*nT +
    43*nV +
    923*nPepBond )

    e280 = 5500*nW  +  1490*nY  + 0.5*125*nC 

    return {'fasta':sequence,'e205':e205,'e214':e214,'e280':e280}



def print_stdout(dict_in):
    #for molid_ind in dict_in['molid_ind_list']:
    for molid_ind in dict_in.keys():
        dict_extn_coeff = dict_in[molid_ind]
        print(dict_extn_coeff.keys())
        molid = dict_extn_coeff['mol_name']
        print("======================================================================================================================================================")
        print("--- Molar absorption coefficient at different wavelength")
        print("mol ID: "+molid)
        print("Sequence: "+dict_extn_coeff['fasta'])
        print( "e(205 nm) = "+str(dict_extn_coeff['e205'])+" (M*cm)^-1 ")
        print( "e(214 nm) = "+str(dict_extn_coeff['e214'])+" (M*cm)^-1 ")
        print( "e(280 nm) = "+str(dict_extn_coeff['e280'])+" (M*cm)^-1 ")
        print("")
        print( "Pace, Vajdos, Fee, Grimsley \"How to measure and predict the molar absorption coefficient of a protein\", Protein Science 1995, 4, 2411-2423")
        print( "Kuipers, Gruppen, \"Prediction of molar extinction coefficients of proteins and peptides ...\", J. Agric. Food Chem. 2007, 55, 5445")
        print( "Anthis, Clore, \"Sequence-specific determination of protein and peptide concentrations by absorbance at 215 nm\", Protein Science 2013, 22, 851")
        return




if __name__ == '__main__':

    # Parse options
    usage = "extn_coeff_fasta.py is the program for calculation of peptide extinction coefficient  (molar absorption coefficient) based on FASTA sequernce.\n"+\
    "\n"+\
    "Usage:          python extn_coeff.py -s GGKGD\n"+\
    ""
    parser = argparse.ArgumentParser(description="extn_coeff_fasta.py is the program for calculation of peptide extinction coefficient  (molar absorption coefficient) based on FASTA sequernce.\n Usage:          python extn_coeff.py -s GGKGD")
    parser.add_argument("-s", action="store", dest="seq", help="peptide sequence", default='')
    parser.add_argument("-i", dest="inputFile", help="input file with molecule structure. smi or sdf",default='')
    parser.add_argument("-o", dest="outputFile", help="output file with molecule structure. fasta",default='')
    parser.add_argument("--json",default=False, action='store_true',dest="l_json", help="Print output as JSON")
    args = parser.parse_args()

    #dict_extn_coeff = calc_extn_coeff(args.sequence)
    dict_extn_coeff = calc_extn_coeff(args.__dict__)

    if not args.l_json:
        print_stdout(dict_extn_coeff)
    else:
        print(json.dumps(dict_extn_coeff))


    






