#!/usr/bin/python
# -*- coding: utf-8 -*

import sys, os, re, string, random, time, urllib
from time import strftime, gmtime
from operator import itemgetter
import argparse
import json
#import commands


from extn_coeff_fasta import calc_extn_coeff
from pI_fasta import calc_pI_fasta
from rdkit_pI import calc_rdkit_pI


currentdir     = os.getcwd()

contact = "If you think this might be a bug. Please contact <a href=\"mailto:andrey.frolov@astrazeneca.com\">Andrey Frolov</a>"



def get_fasta_from_smiles(smi):
    from smi2scrambledfasta import get_scrambledfasta_from_smiles
    fasta = get_scrambledfasta_from_smiles(smi)
 
    if len(fasta)==0: 
        raise Exception('ERROR: returned fasta is empry. something is wrong. Exit')
        sys.exit(1)
    
    return fasta



def arg_parser():

    parser = argparse.ArgumentParser(description="")

    ### common input
    parser.add_argument("--input", dest="input", help="input of a molecule structure: Smiles, fasta, database ID, filename (smi, sdf, fasta)", default='',required=True)

    ### rdkit_pI.py keys
    parser.add_argument("--print_fragment_pkas", dest="l_print_fragment_pkas", help="Print out fragments with corresponding pKas used in pI calcution", default='')
    parser.add_argument("--print_pka_set", dest="l_print_pka_set", help="Print out stored pka sets explicitly.", default='')

    ### pI_fasta.py keys
    parser.add_argument("--ionized_Cterm", dest="ionized_Cterm", help="is C-terminus ionized [COO-]?", default=True)
    parser.add_argument("--ionized_Nterm", dest="ionized_Nterm", help="is N-terminus ionized [N+]?", default=True)
    parser.add_argument("-p", action="store", dest="NPhosphateGroups", help="Number of phosphorilated residues. Phosphorilated residues must be denoted as X in the sequence. default = 0", default=0,type=int)
    parser.add_argument("-l", action="store", dest="NAlkylLysGroups", help="Number of monoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0", default=0,type=int)
    parser.add_argument("-d", action="store", dest="NDiAlkylLysGroups", help="Number of dinoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0", default=0,type=int)

    args = parser.parse_args()

    return args





if __name__ == "__main__":

    args = arg_parser()
    INPUT=args.input

    if not INPUT.isalpha():
        # Assume it is smiles, if contains not only letters
        #print("Input is SMILES")
        smi = INPUT
        molid = 'none'
        fasta = get_fasta_from_smiles(smi)
        print(smi,fasta)
        #dict_out_extn_coeff = calc_extn_coeff(fasta)
        l_calc_extn_coeff=True
        l_calc_pI_fasta=False
        l_calc_rdkit_pI=True
        #l_calc_peptide_cutter_fasta=True

    elif INPUT.isalpha():
        # Assume it is FASTA, if contains only letters
        #print("Input is FASTA")
        fasta = INPUT
        smi = 'none'
        molid = 'none'
        l_calc_extn_coeff=True
        l_calc_pI_fasta=True
        l_calc_rdkit_pI=False
        #l_calc_peptide_cutter_fasta=True

    elif ( (INPUT[0:2] == 'AZ'  and INPUT[2].isdigit()) or (INPUT[0:2] == 'SN'  and INPUT[2].isdigit())  or (INPUT[0:4] == 'MEDI'  and INPUT[4].isdigit())  ):
        # A database ID given
        #print("Input is a database ID")
        smi = get_smiles_from_dbid(INPUT)
        if len(smi) == 0:
            raise Exception("ERROR: could not convert database ID to smiles. Is it corret ID from supported DBs? Exit.")
            sys.exit(1)
        molid = INPUT
        fasta = get_fasta_from_smiles(smi)
        l_calc_extn_coeff=True
        l_calc_pI_fasta=False
        l_calc_rdkit_pI=True
        #l_calc_peptide_cutter_fasta=True


    else:
        raise Exception("ERROR: input not recongnized: not smiles, not fasta, not a known databaase ID. Must be a bug. Contact developer. Exit.")
        sys.exit(1)


    #print("molID: "+molid)
    #print("FASTA: "+fasta)
    #print("SMILES: "+smi)

    dict_out_extn_coeff_fasta = {}
    if l_calc_extn_coeff:
        extn_coeff_options={"seq":fasta,
                         "inputFile":"",
                         "outputFile":"",
                         "l_json":True
                         }
        dict_out_extn_coeff = calc_extn_coeff(extn_coeff_options)

    dict_out_pI_fasta = {}
    if l_calc_pI_fasta:
        # prepare pI_fasta predictor 
    
        if not args.ionized_Cterm: IonizableTerminiOfCTermRes = "\'\'"
        else: IonizableTerminiOfCTermRes = "_"

        if not args.ionized_Nterm: IonizableTerminiOfNTermRes = "\'\'"
        else: IonizableTerminiOfNTermRes = "_"

        pI_fasta_options={"seq":fasta,
                         "inputFile":"",
                         "outputFile":"",
                         "tol": 0.001,
                         "CTermRes": "_",
                         "NTermRes": "_",
                         "IonizableTerminiOfCTermRes": IonizableTerminiOfCTermRes,
                         "IonizableTerminiOfNTermRes": IonizableTerminiOfNTermRes,
                         "lCyclic": False,
                         "NPhosphateGroups": args.NPhosphateGroups,
                         "NAlkylLysGroups": args.NAlkylLysGroups,
                         "NDiAlkylLysGroups": args.NDiAlkylLysGroups,
                         "lPrintpKa": False,
                         "lPlot": True,
                         "lIgnoreC": False,
                         "plot_filename":"OUT_titration_curve.png",
                         "l_json":True
                         }

        dict_out_pI_fasta = calc_pI_fasta(pI_fasta_options)

    dict_out_rdkit_pI = {}
    if l_calc_rdkit_pI:

        rdkit_pI_options={'smiles':smi,
                        'inputFile':'',
                        'outputFile':'',
                        'use_acdlabs':False,
                        'use_dimorphite':True,
                        'l_print_fragments':args.l_print_fragment_pkas,
                        'l_plot_titration_curve':True,
                        'l_print_pka_set':args.l_print_pka_set,
                        'l_json':True
                         }


        dict_out_rdkit_pI = calc_rdkit_pI(rdkit_pI_options)

    dict_out_peptide_tools_master = {'output_extn_coeff':dict_out_extn_coeff,
                                    'output_pI_fasta':dict_out_pI_fasta,
                                    'output_rdkit_pI':dict_out_rdkit_pI
                                    }

    print(json.dumps(dict_out_peptide_tools_master, indent=2))

