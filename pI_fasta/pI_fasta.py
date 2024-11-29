#!/usr/bin/env python
###
### Last  update: Andrey Frolov, AstraZeneca, Molndal. 29/01/2021
### First verion: Andrey Frolov, AstraZeneca, Molndal. 11/02/2016
###

import sys, os
import optparse
import math
from copy import copy

import numpy as np 
from matplotlib.pyplot import *
from itertools import cycle

import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')

from pka_sets_fasta import *

import csv

# Turns a dictionary into a class 
class Dict2Class(object): 
    def __init__(self, my_dict): 
        for key in my_dict: 
            setattr(self, key, my_dict[key]) 


def list_to_comma_seprated_string(l):
	s=""
	for v in l: s+=str(v)+","
	return s[:-1]


# http://www.petercollingridge.co.uk/sites/files/peter/predictPI.txt
#def calculateAminoAcidCharge(amino_acid, pH, pKa):
#    q = charges[amino_acid] 
#    if q>0:
#        return q / (1 + 10**(pH - pKa[amino_acid]))
#    else:
#        return q / (1 + 10**(pKa[amino_acid] - pH))


def calculateBasicAminoAcidCharge(pH, pKa):
        return 1 / (1 + 10**(pH - pKa))

def calculateAcidicAminoAcidCharge(pH, pKa):
        return -1 / (1 + 10**(pKa - pH))

def calculatePhosphateCharge(pH, pKa1, pKa2):
        Ka1=10**(-pKa1)
        Ka2=10**(-pKa2)
        H=10**(-pH)

        f1 = (H*Ka1)/(H**2+H*Ka1+Ka1*Ka2)  # fraction of [AH-]
        f2 = f1 * Ka2 / H                  # fraction of [A2-]
        
        return -2*f2 + (-1)*f1     # average charge of phosphate group




#def calculateProteinCharge(IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes, pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups):
def calculateProteinCharge(pH):

    global seq,IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups, na, nb, lCyclic, lPrintpKaSets, lIgnoreC, tolerance

    protein_charge = 0

    # Middle sequence
    for AA in pKa_basic.keys():
        pKa=pKa_basic[AA][0]
        protein_charge += MiddleSeq.count(AA) * calculateBasicAminoAcidCharge(pH, pKa)
    for AA in pKa_acidic.keys():
        pKa=pKa_acidic[AA][0]
        protein_charge += MiddleSeq.count(AA) * calculateAcidicAminoAcidCharge(pH, pKa)

    # Terminus residues
    for AA in NTermRes + CTermRes:
        if AA in pKa_basic.keys():
            pKa=pKa_basic[AA][1]
            protein_charge += calculateBasicAminoAcidCharge(pH, pKa)
        if AA in pKa_acidic.keys():
            pKa=pKa_acidic[AA][1]
            protein_charge += calculateAcidicAminoAcidCharge(pH, pKa)

			
			
    # Ionizable terminus groups
    for AA in IonizableTerminiOfNTermRes:
        pKa=pKa_TerminusIonizableGroup[AA][0]
        protein_charge += calculateBasicAminoAcidCharge(pH, pKa)

    for AA in IonizableTerminiOfCTermRes:
        pKa=pKa_TerminusIonizableGroup[AA][1]
        protein_charge += calculateAcidicAminoAcidCharge(pH, pKa)
		
		
    # PTMs
    ### Now all phosphorilated AAs have the same pKa s for phosphate group 
    if NPhosphateGroups != 0: protein_charge += NPhosphateGroups * calculatePhosphateCharge(pH, pKa_noncanonical['pKa1_phosphate'], pKa_noncanonical['pKa2_phosphate'])
    if NAlkylLysGroups != 0: protein_charge += NAlkylLysGroups * calculateBasicAminoAcidCharge(pH, pKa_basic['K'][0] + pKa_noncanonical['dpKa_alkylLys'])
    if NDiAlkylLysGroups != 0: protein_charge += NDiAlkylLysGroups * calculateBasicAminoAcidCharge(pH, pKa_basic['K'][0] + pKa_noncanonical['dpKa_dialkylLys'])

    return protein_charge


# Define pH span tocalcualte itration curve and where to search for pI.
def define_pH_span():
    pH_llim=-1
    pH_hlim=15
    return [pH_llim,pH_hlim]


#def calculateIsoelectricPoint(IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups):   
def calculateIsoelectricPoint():   

    global seq,IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups, na, nb, lCyclic, lPrintpKaSets, lIgnoreC, tolerance

    #tolerance=0.01
    #charge_tol=0.05
    #min_pH, max_pH = 0 , 14 
    pH_lim = define_pH_span()
    min_pH0 = pH_lim[0]
    max_pH0 = pH_lim[1]

    min_pH = min_pH0
    max_pH = max_pH0

    while True:
        mid_pH = 0.5 * (max_pH + min_pH)
        #protein_charge = calculateProteinCharge(IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups)
        protein_charge = calculateProteinCharge(mid_pH)
        
        if na == 0 and nb != 0:
            refcharge = charge_tol * nb

        elif nb == 0 and na != 0:
            refcharge = -charge_tol * na

        else:
            refcharge = 0.0

        if protein_charge > refcharge + tolerance:
            min_pH = mid_pH
        elif protein_charge < refcharge - tolerance:
            max_pH = mid_pH
        else:
            return mid_pH
            
        if mid_pH <= min_pH0:
            return min_pH0
        elif mid_pH >= max_pH0:
            return max_pH0
            

#def CalcChargepHCurve(IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes):
#def CalcChargepHCurve():
#   from numpy import arange
#   global seq,IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups, na, nb, lCyclic, lPrintpKaSets, lIgnoreC, tolerance
#   
#   pH_a=arange(0,14,0.1)
#   Q_a=pH_a*0.0    

#   for i in range(len(pH_a)):
#       #Q = calculateProteinCharge(IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,pH_a[i])
#       Q = calculateProteinCharge(pH_a[i])
#       Q_a[i]=Q
#       
#   return pH_a, Q_a


def CalcChargepHCurve():
    global seq,IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups, na, nb, lCyclic, lPrintpKaSets, lIgnoreC, tolerance
    pH_lim = define_pH_span()
    pH_a=np.arange(pH_lim[0],pH_lim[1],0.1)
    Q_a=pH_a*0.0    
    for i in range(len(pH_a)):
#        Q = calculateMolCharge(base_pkas, acid_pkas, diacid_pkas, pH_a[i],constant_q=constant_q)
        Q = calculateProteinCharge(pH_a[i])
        Q_a[i]=Q
    pH_Q = np.vstack((pH_a,Q_a))
    return pH_Q


def separateTerminalRes(sequence):
    NTermRes = sequence[0]
    CTermRes = sequence[-1]
    MiddleSeq = sequence[1:-1]
    return NTermRes,MiddleSeq,CTermRes
    

def split_sequence(sequence):
    
        #global IonizableTerminiOfNTermRes, NTermRes, CTermRes, IonizableTerminiOfCTermRes
        global seq,IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups, na, nb, lCyclic, lPrintpKaSets, lIgnoreC, tolerance
        global titAdd
        titAdd=""

        if not lCyclic:
	### Assume linear, not branched peptide sequence

                ### If custom 
#                if NTermRes == '_' and CTermRes == '_' and (IonizableTerminiOfNTermRes == '_' or IonizableTerminiOfNTermRes == '') and (IonizableTerminiOfCTermRes == '_' or IonizableTerminiOfCTermRes == ''):
#                	MiddleSeq = sequence[1:-1]
#                else:
#                        print("---NOTE! custom termini residues and/or ionizable termini given. No termini residues would be identified from the given sequence. Custom given residues will be used instead.")
#                        if NTermRes == '_': 
#                            print('---Error! custom termini specified => NTermRes thus must be specified explicitly ')
#                            sys.exit(1)
#                       if CTermRes == '_': 
#                            print('---Error! custom termini specified => CTermRes thus must be specified explicitly ')
#                            sys.exit(1)
#                        if IonizableTerminiOfNTermRes == '_':
#                            print('---Error! custom termini specified => IonizableTerminiOfNTermRes thus must be specified explicitly ')
#                            sys.exit(1)
#                        if IonizableTerminiOfCTermRes == '_': 
#                            print('---Error! custom termini specified => IonizableTerminiOfCTermRes thus must be specified explicitly ')
#                            sys.exit(1)
#                        MiddleSeq = sequence

                MiddleSeq = sequence[1:-1]
                
                if IonizableTerminiOfNTermRes == '_': IonizableTerminiOfNTermRes = sequence[0]
                else: 
                    if IonizableTerminiOfNTermRes == '': titAdd+=", capped N-terminus"	
                    else: titAdd+=", custom ionizable termini -NH2: "+IonizableTerminiOfNTermRes	

                if NTermRes == '_': NTermRes = sequence[0]
                else: titAdd+=", custom N termini residues: "+NTermRes

                if CTermRes == '_': CTermRes = sequence[-1]
                else: titAdd+=", custom C termini residues: "+CTermRes

                if IonizableTerminiOfCTermRes == '_': IonizableTerminiOfCTermRes = sequence[-1]
                else: 
                    if IonizableTerminiOfCTermRes == '': titAdd+=", capped C-terminus"	
                    else: titAdd+=", custom ionizable termini -COOH: "+IonizableTerminiOfCTermRes	

        else:
	### Assume cyclic peptide: no ionizable termini, no terminal residues 
            MiddleSeq = sequence
            titAdd+=", cyclic"
            NTermRes=''
            CTermRes=''
            IonizableTerminiOfNTermRes=''
            IonizableTerminiOfCTermRes=''	

        if NPhosphateGroups > 0: titAdd += ', '+str(NPhosphateGroups)+' phosphorilated res'
        if NAlkylLysGroups > 0: titAdd += ', '+str(NAlkylLysGroups)+' monoalkyl Lys'
        if NDiAlkylLysGroups > 0: titAdd += ', '+str(NDiAlkylLysGroups)+' dialkyl Lys'

        return IonizableTerminiOfNTermRes,NTermRes,MiddleSeq,CTermRes,IonizableTerminiOfCTermRes
    



def mean(lst):
    """calculates mean"""
    return sum(lst) / len(lst)

def stddev(lst):
    """returns the standard deviation of lst"""
    mn = mean(lst)
    variance = sum([(e-mn)**2 for e in lst])
    return math.sqrt(variance)

def stddev(lst):
    """returns the standard deviation of lst"""
    mn = mean(lst)
    variance = sum([(e-mn)**2 for e in lst])
    return math.sqrt(variance)

def stderr(lst):
    """returns the standard error of the mean of lst"""
    mn = mean(lst)
    variance = sum([(e-mn)**2 for e in lst])
    return math.sqrt(variance) / math.sqrt(len(lst))

def print_pka_set():
        print()
        print()
        print('----------------------------------------------------------------------------')
        print('--- Used pKa values for each set: (http://isoelectric.ovh.org/theory.html) ')
        print(json.dumps(pKa_sets_short, indent=2))
        print() 
        print('----------------------------------------------------------------------------')
        print('For the following sets the pKa depends on the Residue poistion in the sequence, also the termini group pKa change among different terminal Residue.')
        print('See references details. ')
#        print 'Therefore 3 pKa values are given for each residue: middle, N-terminal, C-terminal positions in the sequence'
#        print 'Also, the pKas of the termini -NH2 and -COOH depend on the residue. Thus for each residue 2 values are given: pKa of N-term and pKa of C-term. '
#        print 
        print('ProMoST set: see http://isoelectric.ovh.org/theory.html and http://proteomics.mcw.edu/promost_adv.html  for details')
#        print 
#        print json.dumps(pKa_sets['ProMoST'], indent=2)
    
#        print
        print('Gauci_calibrated set: see Gauci et al. Proteomics 2008, 8, 4898  and  https://github.com/ypriverol/pIR   for details')
#        print json.dumps(pKa_sets['Gauci_calib'], indent=2)
        
    
        print('--------------------------------------------------------------' )
        print('Supported nonatural aminoacids (teh same for all sets of pKa):')
        print()
        print('   Phosphate pKa1  '+str(pKa_noncanonical['pKa1_phosphate']) )
        print ('   Phosphate pKa2  '+str(pKa_noncanonical['pKa2_phosphate']) )
        print() 
        print('   Alkylated Lys: addition to Lys pKa  '+str(pKa_noncanonical['dpKa_alkylLys'])+ "      data from ACD lab: pKa of amine: 10.69. The delta for methylated amine compared to amine.  ### Zhang, Vogel, J. Bio. Chem. 1993, 268, 30, 22420 (Table III, Lys75) pKas of methylated 10.87, dimethylated 10.12" )
        print() 
        print('   Dialkylated Lys: addition to Lys pKa  '+str(pKa_noncanonical['dpKa_dialkylLys'])+"       data from Zhang, Vogel et al.  (ACD lab: pKa of dimethylamine: 9.83 +- 0.28 - error too high. The delta for methylated amine compared to amine.  ### Zhang, Vogel, J. Bio. Chem. 1993, 268, 30, 22420 (Table III, Lys75) pKas of methylated 10.87, dimethylated 10.12" )
        print()

        return




def print_output(dict_pI_fasta,args):
    for molid_ind in dict_pI_fasta.keys():
        dict_single = dict_pI_fasta[molid_ind]
        print_output_dict(dict_single['pI'],'pI',dict_single['plot_title_info']) 
        print_output_dict(dict_single['QpH7'],'Q at pH7.4',dict_single['plot_title_info']) 
        if args.lPrintpKa: print_pka_set()
    return
    

def print_output_dict(out_dict,prop,title_info):
    #global tit
    lj=12
    keys=list(out_dict.keys())
    keys.remove('std'); keys.insert(0, 'std')
    keys.remove('err'); keys.insert(0, 'err')
    keys.remove(prop + ' mean'); keys.insert(0, prop+' mean')

    #tit="sequence: "+seq+titAdd
    p=out_dict
    print(" ")
    print("======================================================================================================================================================")
    print(prop+" for " + title_info  )
    print( "---------------------------------")
    for k in keys:
        print(k.rjust(lj)  + "  " +  str(round(p[k],2)).ljust(lj) )
    print(" ")

    return








def options_parser():
    # Parse options
    usage = """pI_fasta.py is the program for calculation of peptide isoelectric points using Henderson-Hasselbalch equations. Various sets of pKa data are supported (see below).

Sequence input:         python pI_fasta.py -s GGKGD
FASTA file input:       python pI_fasta.py -i example.fasta
With plot:              python pI_fasta.py -s GGKGD -x 
Cyclic peptide:         python pI_fasta.py -s GGKGD -r -x
Capped N terminus:      python pI_fasta.py -s GGKGD -b \"\" -x
Capped C terminus:      python pI_fasta.py -s GGKGD -a \"\" -x
Phosphorylated residue: python pI_fasta.py -s GXD -p 1 -x
Monoalkylated Lys:      python pI_fasta.py -s GXD -l 1 -x
Dialkylated Lys:        python pI_fasta.py -s GXD -d 1 -x
Branched peptide (custom terminal residues):
                        python pI_fasta.py -s GGKGD   -c AX -n E -a AX -b E -x
Use custom pKa set:     python pI_fasta.py -s GGKGD -m IPC2_peptide,ProMoST,Gauci -x
Help:                   python pI_fasta.py -h
Most extended. All defaults listed. Mind: \"_\" indicates that the residue is automatically deduced from the given sequence:
                        python pI_fasta.py  -s GGKGD  -t 0.001  -c _  -n _  -a _  -b _  -p 0 -l 0 -d 0 -m IPC2_peptide,IPC_peptide,ProMoST,Gauci,Grimsley,Thurlkill,Lehninger,Toseland

--- JSON fomatting
JSON formated input:    python pI_fasta.py -g \'{\"1\":{\"mol_name\":\"name1\",\"fasta\":\"GGKGD\"}}\'
JSON formated output:   python pI_fasta.py -s GGKGD -j

Nested structure of JSON, top level - unique IDs of the molecules, bottom level data fields for each molecule "mol_name" and "fasta"

See help for all the options.

--- For theory and pKa sets see: http://isoelectric.ovh.org/theory.html
For ProMoST and Gauchi sets pKa depends on whether the residue sits at the termini or in the middle of the sequence
Also, the pKa of ionizable termini depends on the type of termini residue.
The rest of the sets position invariant.
for ProMoST set: see http://proteomics.mcw.edu/promost_adv.html
for Gauci set: see Gauci et al. Proteomics 2008, 8, 4898  and  https://github.com/ypriverol/pIR
Supported nonatural aminoacids (the same for all sets of pKa): Monoalkylated Lys, DiAlkylated Lys, phosphorilations
derived from literature, ACDlab predictions, see -z for more details.

Andrey Frolov, AstraZeneca, Molndal. 04/03/2016
Last updated 25/10/2022

"""

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-s", action="store", dest="seq", help="peptide sequence", default="")
    parser.add_option("-i", action="store", dest="inputFile", help="input file name in fasta format", default="")
    parser.add_option("-g", action="store", dest="inputJSON", help="input file name in JSON format", default="")
    parser.add_option("-o", action="store", dest="outputFile", help="output file name in csv format", default="")
    parser.add_option("-t", action="store", type="float", dest="tol", help="tolerance on total protein charge. default = 0.001", default=0.001)

    parser.add_option("-c", action="store", type='string',dest="CTermRes", help="Custom list of C terminus residues. By default it is set to the last residues of the given sequence. This option is useful if you have a branched peptide with several terminus residues", default='_')
    parser.add_option("-n", action="store", type='string',dest="NTermRes", help="Custom list of N terminus residues. By default it is set to the first residues of the given sequence. This option is useful if you have a branched peptide with several terminus residues", default='_')

    parser.add_option("-a", action="store", type='string',dest="IonizableTerminiOfCTermRes", help="Custom list of residues with ionizable C terminus (-COOH). By default it is set to the last residue of the given sequence. This option is useful if you want to cap the termini (exclude them from calculations) or you have a branched peptide with several terminus residues", default='_')
    parser.add_option("-b", action="store", type='string',dest="IonizableTerminiOfNTermRes", help="Custom list of residues with ionizable N terminus (-NH2). By default it is set to the first residue of the given sequence. This option is useful if you want to cap the termini (exclude them from calculations) or you have a branched peptide with several terminus residues", default='_')


    parser.add_option("-p", action="store", type='int',dest="NPhosphateGroups", help="Number of phosphorilated residues. Phosphorilated residues must be denoted as X in the sequence. default = 0", default=0)
    parser.add_option("-l", action="store", type='int',dest="NAlkylLysGroups", help="Number of monoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0", default=0)
    parser.add_option("-d", action="store", type='int',dest="NDiAlkylLysGroups", help="Number of dialkylated Lys residues. These residues should be denoted as X in the sequence. default = 0", default=0)

    parser.add_option("-m", action="store", type='string',dest="pka_set_list", help="List of pKa sets to use in calculation (comma separated). default = "+list_to_comma_seprated_string(pKa_sets_to_use), default='')
    parser.add_option("-r", action="store_true", dest="lCyclic", help="Is it cyclic? No termini residues are derived from sequence and also no ionizable termini (-NH2 and -COOH) of the terminal residues are derived from sequence (However, custom values are not overwritten). default = False", default=False)
    parser.add_option("-x", action="store_true", dest="lPlot", help="plot charge/pH curves. Requires NumPy and Matplotlib. default = False", default=False)
    parser.add_option("-q", action="store_true", dest="lIgnoreC", help="ignore cysteins, assume they are protected. default = False", default=False)
    parser.add_option("-z", action="store_true", dest="lPrintpKa", help="print used pKa values in the output. default = False", default=False)
    parser.add_option("-j", action="store_true", dest="l_json", help="use JSON as an output format. default = False", default=False)

    (options, args) = parser.parse_args()

    return options.__dict__
 




### PLOT titration curve
def plot_titration_curve(pH_Q_dict,figFileName):
    matplotlib.rcParams.update({'font.size': 16})
    lines = ["-","--","-.",":"]
    w1=4.0 ; w2=3.0 ; w3=2.0 ; w4=1.0
    linew = [w1,w1, w2,w2, w3,3, w4,w4]
    linecycler = cycle(lines)
    linewcycler = cycle(linew)

    figure(figsize=(8,6))
    i=0
    for pKaset in pKa_sets_to_use:
        i+=1
        pH_Q = pH_Q_dict[pKaset] 
        l=plot(pH_Q[:,0],pH_Q[:,1],next(linecycler),label=pKaset,linewidth=next(linewcycler)) 
        if pKaset == 'IPC2_peptide': 
            setp(l,linewidth=8,linestyle='-',color='k')

        # Store data for output
        if i==1: 
            pH = pH_Q[:,0]
            Q_M = pH_Q[:,1]
        else:
            Q_M = np.column_stack([Q_M,pH_Q[:,1]])

    plot(pH,pH*0,'k-')
    plot([7,7],[np.min(Q_M),np.max(Q_M)],'k-')
    #xlim([np.min(pH),np.max(pH)])
    xlim([2,12])
    ylim([np.min(Q_M),np.max(Q_M)])

    legend(loc="center right", bbox_to_anchor=[1.1, 0.5],ncol=1, shadow=True, fontsize=10).get_frame().set_alpha(1)
    ylabel('peptide charge')
    xlabel('pH')
	
    title('Titration curve')    
	
    minorticks_on()
    grid(True)

    #show()
    savefig(figFileName)
    return




#def plot_titration_curve(fig_file_name='OUT_titration_curve.png'):
#       from numpy import arange, column_stack
#       from matplotlib.pyplot import plot, figure,setp,legend,ylabel,xlabel,title,minorticks_on,grid,savefig
#       #matplotlib.rcParams.update({'font.size': 16})
#       from itertools import cycle
#       global seq,IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups, na, nb, lCyclic, lPrintpKaSets, lIgnoreC, tolerance, tit
#       lines = ["-","--","-.",":"]
#       w1=4.0
#       w2=3.0
#       w3=2.0
#       w4=1.0
#       linew = [w1,w1, w2,w2, w3,3, w4,w4]
#       linecycler = cycle(lines)
#       linewcycler = cycle(linew)

#       figure()

#       i=0
#       for pKaset in pKa_sets_to_use:
#           i+=1
#           pKa_basic=pKa_sets[pKaset]['pKa_basic']
#           pKa_acidic=pKa_sets[pKaset]['pKa_acidic']
#           pKa_TerminusIonizableGroup=pKa_sets[pKaset]['pKa_TerminusIonizableGroup']

#           #pH,Q = CalcChargepHCurve(IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes) 
#           pH,Q = CalcChargepHCurve() 
#  
#           l=plot(pH,Q,next(linecycler),label=pKaset,linewidth=next(linewcycler)) 
#           if pKaset == 'ProMoST': 
#               setp(l,linewidth=8,linestyle='-',color='k')

#           # Store data for output
#           if i==1: Q_M = pH
#           Q_M=column_stack([Q_M,Q])


#       plot(pH,pH*0,'k-')
#       plot([7,7],[min(Q),max(Q)],'k-')

#       legend(loc="center right", bbox_to_anchor=[1.1, 0.5],ncol=1, shadow=True).get_frame().set_alpha(1)
#       ylabel('peptide charge')
#       xlabel('pH')
#       
#       title(tit)    
#       
#       minorticks_on()
#       grid(True)

#       #show()
#       #savefig("OUT_titration_curve.png")
#       savefig(fig_file_name)
#       #pltsave("OUT_titration_curve", ext="png", close=True, verbose=True)

#       return



def read_fasta_file(inputFile):

    filename, ext = os.path.splitext(inputFile)

    # Initialize file reader
    if not ext == '.fasta': raise Exception('!Warning: extension of file is not ".fasta". Assuming it is fasta formatted input. Continue. ')

    from Bio import SeqIO
    biosuppl = SeqIO.parse(open(inputFile),'fasta')

    mol_supply_json={}
    mol_unique_ID = 0
    for biofasta in biosuppl:
        mol_unique_ID += 1
        # unique index, mol title, RDkit mol object, mol fasta
        mol_supply_json[mol_unique_ID] = {'mol_name':biofasta.id, 'mol_obj':None, 'fasta':str(biofasta.seq)}
        
    return mol_supply_json



def calc_pI_fasta(options={"inputFile":"","inputDict":{},"inputJSON":"","outputFile":"","sequenceList":"","seq":"", "tol": 0.001, "CTermRes": "_", "NTermRes": "_", "IonizableTerminiOfCTermRes": "_", "IonizableTerminiOfNTermRes": "_", "lCyclic": False, "NPhosphateGroups": 0, "NAlkylLysGroups": 0, "NDiAlkylLysGroups": 0, "lPrintpKa": False, "lPlot": False, "lIgnoreC": False,"plot_filename":"OUT_titration_curve.png","l_json":False}):

    args = Dict2Class(options)

    # Get options
    if len(args.seq)!=0:
        # assume single fasta input
        mol_unique_ind=1
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
        raise Exception('Error: either fasta, input file *.fasta or JSON should be given for pI_fasta.py. Exit. ')
        sys.exit(1)



    dict_out_pI_fasta = {}
    for mol_unique_ind in mol_supply_json.keys():
        options_single = copy(options)
        options_single["seq"] = mol_supply_json[mol_unique_ind]['fasta']
        dict_pI_fasta_single = calc_pI_fasta_single_sequence(options_single)
        dict_pI_fasta_single["mol_name"] = mol_supply_json[mol_unique_ind]['mol_name']
        dict_out_pI_fasta[mol_unique_ind] = dict_pI_fasta_single

    return dict_out_pI_fasta



def calc_pI_fasta_single_sequence(options={"seq":"", "tol": 0.001, "CTermRes": "_", "NTermRes": "_", "IonizableTerminiOfCTermRes": "_", "IonizableTerminiOfNTermRes": "_", "lCyclic": False, "NPhosphateGroups": 0, "NAlkylLysGroups": 0, "NDiAlkylLysGroups": 0, "lPrintpKa": False, "lPlot": False, "lIgnoreC": False,"plot_filename":"OUT_titration_curve.png","l_json":False,"pka_set_list":''}):

    global seq,IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes,  mid_pH,pKa_basic,pKa_acidic,pKa_TerminusIonizableGroup, NPhosphateGroups, NAlkylLysGroups, NDiAlkylLysGroups, na, nb, lCyclic, lPrintpKaSets, lIgnoreC, tolerance, tit, pKa_sets_to_use

    if len(options['seq']) == 0: raise Exception("no sequence given. See --help for more information")

    if options['pka_set_list'] != '': pKa_sets_to_use = options['pka_set_list'].split(',')

    for set in pKa_sets_to_use:
            if set not in all_known_pKa_sets: raise Exception("---Error! pKa set "+set+" is not known. Check your -m option. Exit.")

    # Get options
    orig_seq = options['seq']
    seq=orig_seq.upper()
    sequence=seq
    IonizableTerminiOfNTermRes=options['IonizableTerminiOfNTermRes']
    IonizableTerminiOfCTermRes=options['IonizableTerminiOfCTermRes']
    NTermRes=options['NTermRes']
    CTermRes=options['CTermRes']
    tolerance=options['tol']
    NPhosphateGroups=options['NPhosphateGroups']
    NAlkylLysGroups=options['NAlkylLysGroups']
    NDiAlkylLysGroups=options['NDiAlkylLysGroups']
    lIgnoreC=options['lIgnoreC']
    lPrintpKaSets=options['lPrintpKa']
    lPlot=options['lPlot']
    lCyclic=options['lCyclic']
    
    lsaveplotdata=False
    lCalc=True
    

    if len(seq) < 2:
        raise Exception("---!Error: number of residues in sequence is less than 2. Not yet supported. Exiting.")
        sys.exit(1)

    for R in seq:
        if R not in known_res: 
            raise Exception("---!Error: residue "+R+" is not known. Please use X if this is a noncaninical residue. Exiting.")
            sys.exit(1)
 

    if lIgnoreC:
        sequence = sequence.replace('C','X')



    ##### Calc pI

    IonizableTerminiOfNTermRes, NTermRes, MiddleSeq, CTermRes, IonizableTerminiOfCTermRes = split_sequence(sequence)

    nb=0
    na=0
    for R in NTermRes + MiddleSeq + CTermRes:
        if R in known_basic_res: nb+=1
        if R in known_acidic_res: na+=1
    for R in IonizableTerminiOfNTermRes: 
            nb+=1
    for R in IonizableTerminiOfCTermRes: 
            na+=1
	
    if na == 0 and nb == 0:
        print("---!Warning: no ionizable groups in the sequence. pI is not defined or = any pH. Exiting.")
        sys.exit(0)

    charge_tol=0.05
    if na == 0 and nb != 0:
        print("---!Warning: no acidic ionizable groups, only basic groups present in the sequence. pI is not defined and thus won't be calculated. However, will estimate pH when peptide has charge less than "+str(charge_tol*100)+"% of the peptide maximum possible charge (by absolute value). Continue.")
        #lPlot = True
        #lCalc = True

    if nb == 0 and na != 0:
        print("---!Warning: no basic ionizable groups, only basic groups present in the sequence. pI is not defined and thus won't be calculated. However, will estimate pH when peptide has charge less than "+str(charge_tol*100)+"% of the peptide maximum possible charge (by absolute value). Continue.")
        #lPlot = True
        #lCalc = True




    tit="sequence: "+orig_seq+titAdd
    ### Calculate pI
    if lCalc:

        seq_dict={}
        pI_dict={}
        Q_dict={}
        pH_Q_dict={}

        #for pKaset in pKa_sets.keys():
        for pKaset in pKa_sets_to_use:
    
            pKa_basic=pKa_sets[pKaset]['pKa_basic']
            pKa_acidic=pKa_sets[pKaset]['pKa_acidic']
            pKa_TerminusIonizableGroup=pKa_sets[pKaset]['pKa_TerminusIonizableGroup']

            pI = calculateIsoelectricPoint()
            pI_dict[pKaset] = pI   

            Q_pH74 = calculateProteinCharge(7.4)
            Q_dict[pKaset] = Q_pH74

            #pH_Q = CalcChargepHCurve(all_base_pkas, all_acid_pkas, all_diacid_pkas, constant_q = molecule_constant_charge)
            pH_Q = CalcChargepHCurve()
            pH_Q = pH_Q.T

            pH_Q_dict[pKaset] = pH_Q


        # pI 
        pIl=[]
        for k in pI_dict.keys(): pIl += [pI_dict[k]]
        pI_dict['pI mean']=mean(pIl)
        pI_dict['std']=stddev(pIl)
        pI_dict['err']=stderr(pIl)

        # charge at pH=7.4
        Ql=[]
        for k in Q_dict.keys(): Ql += [Q_dict[k]]
        Q_dict['Q at pH7.4 mean']=mean(Ql)
        Q_dict['std']=stddev(Ql)
        Q_dict['err']=stderr(Ql)

        # print isoelectric interval
        pKaset='IPC2_peptide'
        int_tr = 0.2    # TODO define it elsewhere 
        pH_Q = pH_Q_dict[pKaset]
        Q=pH_Q[:,1]
        pH=pH_Q[:,0]
        pH_int = ( pH[(Q>-int_tr) & (Q<int_tr)] )
        pH_Q = pH_Q_dict[pKaset]

        # isoelectric interval - pH range where the charge is within the given threshold. If molecule permanently has a charge the interval is not defined and NaN are provided. 
        if len(pH_int) > 1:
            interval = (pH_int[0], pH_int[-1])
        else:
            interval = (float('NaN'), float('NaN'))

    ### Plot titration curve
    plot_filename = ''
    if lPlot:
        if "plot_filename" in options.keys():
            plot_filename = options["plot_filename"]
        else:
            plot_filename = "OUT_titration_curve_pI_fasta.png"
            
        #plot_titration_curve(pH_Q_dict,fig_file_name=plot_filename)
        plot_titration_curve(pH_Q_dict,plot_filename)
            

    dict_pI_fasta = {'sequence':orig_seq,'pI':pI_dict,'QpH7':Q_dict,'pI_interval':interval,'plot_filename':plot_filename,'plot_title_info':tit}

    return dict_pI_fasta




#####======================================================================================================================================================

if __name__ == "__main__":

    options = options_parser()
    args = Dict2Class(options)
    
    dict_out_pI_fasta = calc_pI_fasta(options)
    
    ### ----------------------------------------------------------------------
    # Output 
    if args.outputFile == '': # output plain text
        if args.l_json:
            print(json.dumps(dict_out_pI_fasta, indent=2))
        else:
            print_output(dict_out_pI_fasta,args)    

    else: # output file

        known_out_file_types = ['.csv']
        filename, out_fext = os.path.splitext(args.outputFile)
        if out_fext not in known_out_file_types:
            raise Exception('Error! Output file extention not in supported file types:'+str(known_file_types))
            sys.exit(1)

        elif out_fext == '.csv':
            with open(args.outputFile,'w') as csv_f:
                csv_w = csv.writer(csv_f)
                count=0
                for mi in dict_out_pI_fasta.keys():
                    count+=1
                    if count == 1:
                        header = ['mol_name','fasta','pI mean','pI std','pI interval']
                        csv_w.writerow(header)

                    row = []
                    row += [dict_out_pI_fasta[mi]['mol_name']] 
                    row += [dict_out_pI_fasta[mi]['sequence']] 
                    row += [dict_out_pI_fasta[mi]['pI']['pI mean']] 
                    row += [dict_out_pI_fasta[mi]['pI']['std']]         
                    row += [ ' - '.join([ "%.2f" % x for x in dict_out_pI_fasta[mi]['pI_interval'] ]) ] 
                    csv_w.writerow(row)
                        
        # print info 
        dict_file = {'outputFile':args.outputFile,'outputInfo':'Number of molecules processed:'+str(len(dict_out_pI_fasta.keys()))}
        print(json.dumps(dict_file))
   


