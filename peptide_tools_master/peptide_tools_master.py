#!/usr/bin/python
# -*- coding: utf-8 -*

import sys, os, re, string, random, time, urllib
from time import strftime, gmtime
from operator import itemgetter
import argparse
import json
#import commands
import csv

from rdkit import Chem

from extn_coeff_fasta import calc_extn_coeff
from pI_fasta import calc_pI_fasta
from rdkit_pI import calc_rdkit_pI

import tempfile

currentdir     = os.getcwd()

contact = "If you think this might be a bug. Please contact <a href=\"mailto:andrey.frolov@astrazeneca.com\">Andrey Frolov</a>"



def get_fasta_from_smiles(smi):
    from smi2scrambledfasta import get_scrambledfasta_from_smiles
    fasta = get_scrambledfasta_from_smiles(smi)
    if len(fasta)==0: 
        raise Exception('ERROR: returned fasta is empry. something is wrong. Exit')
        sys.exit(1)
    return fasta


def get_fasta_from_mol(mol):
    from smi2scrambledfasta import get_scrambledfasta_from_mol
    fasta = get_scrambledfasta_from_mol(mol)
    if len(fasta)==0: 
        raise Exception('ERROR: returned fasta is empry. something is wrong. Exit')
        sys.exit(1)
    return fasta



def arg_parser():

    parser = argparse.ArgumentParser(description="")

    ### common input
    parser.add_argument("--input", dest="input", help="input of a molecule structure: Smiles, fasta, database ID, filename (smi, sdf, fasta)", default='',required=True)

    ### rdkit_pI.py keys
    #parser.add_argument("--print_fragment_pkas", dest="l_print_fragment_pkas", help="Print out fragments with corresponding pKas used in pI calcution", default='')
    #parser.add_argument("--print_pka_set", dest="l_print_pka_set", help="Print out stored pka sets explicitly.", default='')
    parser.add_argument("--print_fragment_pkas",default="no", action='store',dest="l_print_fragment_pkas", help="Print out fragments with corresponding pKas used in pI calcution")
    parser.add_argument("--print_pka_set",default="no", action='store',dest="l_print_pka_set", help="Print out stored pka sets explicitly.")

    ### pI_fasta.py keys
    parser.add_argument("--ionized_Cterm", dest="ionized_Cterm", help="is C-terminus ionized [COO-]?", default=True)
    parser.add_argument("--ionized_Nterm", dest="ionized_Nterm", help="is N-terminus ionized [N+]?", default=True)
    parser.add_argument("-p", action="store", dest="NPhosphateGroups", help="Number of phosphorilated residues. Phosphorilated residues must be denoted as X in the sequence. default = 0", default=0,type=int)
    parser.add_argument("-l", action="store", dest="NAlkylLysGroups", help="Number of monoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0", default=0,type=int)
    parser.add_argument("-d", action="store", dest="NDiAlkylLysGroups", help="Number of dinoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0", default=0,type=int)

    args = parser.parse_args()

    if args.l_print_pka_set == "yes": args.l_print_pka_set = True
    else: args.l_print_pka_set = False
    
    if args.l_print_fragment_pkas == "yes": args.l_print_fragment_pkas = True
    else: args.l_print_fragment_pkas = False


    return args




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



def read_structure_file(inputFile):

    filename, ext = os.path.splitext(inputFile)

    # Initialize file reader
    if   ext == '.sdf':
        suppl = Chem.SDMolSupplier(inputFile) 
    elif ext == '.smi':
        suppl = Chem.SmilesMolSupplier(inputFile,titleLine=False) 
    else:
        raise Exception('!Warning: extension of file is not smi or sdf. Assume it is smi. Continue. ')
        suppl = Chem.SmilesMolSupplier(inputFile,titleLine=False) 

    mol_supply_json={}
    mol_unique_ID = 0
    for mol in suppl:
        mol_unique_ID += 1
        # unique index, mol title, fasta
        #fasta = get_fasta_from_smiles(smi)

        if not mol.HasProp('_Name'): mol.SetProp('_Name','tmpname'+str(mol_unique_ID))

        mol_supply_json[mol_unique_ID] = {'mol_name':mol.GetProp('_Name'), 'mol_obj':mol, 'fasta':get_fasta_from_mol(mol)}

    return mol_supply_json



if __name__ == "__main__":

    args = arg_parser()
    INPUT=args.input
    INPUT=INPUT.strip()

    mol_name = 'none'
    IN_lines = INPUT.split('\n')

    if len(IN_lines) == 1:
        IN_vals = IN_lines[0].split()
        if len(IN_vals) == 1: 
            INPUT = IN_vals[0]
        else:
            INPUT = IN_vals[0]
            mol_name = IN_vals[1]

    elif len(IN_lines) == 0 :
        # Houston, we have a problem
        # No INPUT, nothing to do
        raise Exception('ERROR: no input in '+sys.argv[0])
        sys.exit(1)

    else:
        # Multiple string input
        if IN_lines[0][0] == '>' or IN_lines[0][0] == ';': # check https://en.wikipedia.org/wiki/FASTA_format
        # Assuming the multiple rows are the FASTA file input 
            tf = tempfile.NamedTemporaryFile(prefix='tmp_peptide_tools_master',suffix='.fasta',delete=True)
            INPUT = tf.name
            with open(INPUT,'w') as f:
                for line in IN_lines:
                    f.write(line+'\n')

        elif '$$$$' in INPUT: 
        # Assuming the multiple rows are SDF file
            tf = tempfile.NamedTemporaryFile(prefix='tmp_peptide_tools_master',suffix='.sdf',delete=True)
            INPUT = tf.name
            with open(INPUT,'w') as f:
                for line in IN_lines:
                    f.write(line+'\n')

#       elif not IN_lines[0].split()[0].isalpha(): 
#       # Assuming the multiple rows are the smiles 
#           tf = tempfile.NamedTemporaryFile(prefix='tmp_peptide_tools_master',suffix='.smi',delete=True)
#           INPUT = tf.name
#           with open(INPUT,'w') as f:
#               for line in IN_lines:
#                   f.write(line+'\n')

        else:
        # Assuming the multiple rows are the smiles 
            tf = tempfile.NamedTemporaryFile(prefix='tmp_peptide_tools_master',suffix='.smi',delete=True)
            INPUT = tf.name
            with open(INPUT,'w') as f:
                for line in IN_lines:
                    f.write(line+'\n')


    known_file_types = ['.sdf','.smi','.smiles','.fasta']

    
    lPlot=True
    mol_supply_json = {}

    if os.path.exists(INPUT):
        # Assume input is a file.

        lPlot=False

        filename, file_extension = os.path.splitext(INPUT)
        if file_extension not in known_file_types:
            raise Exception('Error! File extention not in supported file types:'+str(known_file_types))
            sys.exit(1)
        else:
            #smi = ''
            #mol_name = 'none'
            inputFile = INPUT
            if file_extension == '.smi' or file_extension == '.smiles' or file_extension == '.csv' or file_extension == '.fasta':
                out_fext = '.csv'
            elif file_extension == '.sdf':
                out_fext = '.sdf'
                
            outputFile = filename + '_OUTPUT' + out_fext
            if file_extension != '.fasta':
                l_calc_extn_coeff=True
                l_calc_pI_fasta=False
                l_calc_rdkit_pI=True
                mol_supply_json = read_structure_file(inputFile)
            else:
                l_calc_extn_coeff=True
                l_calc_pI_fasta=True
                l_calc_rdkit_pI=False
                mol_supply_json = read_fasta_file(inputFile)

    elif not INPUT.isalpha():
        # Assume it is smiles, if contains not only letters
        #print("Input is SMILES")
        mol_unique_ID = 1
        smi = INPUT
        #mol_name = 'none'
        fasta = get_fasta_from_smiles(smi)
        inputFile = ''
        outputFile = ''
        l_calc_extn_coeff=True
        l_calc_pI_fasta=False
        l_calc_rdkit_pI=True
        mol = Chem.MolFromSmiles(smi)
        mol_supply_json[mol_unique_ID] = {'mol_name': mol_name, 'mol_obj':mol, 'fasta':get_fasta_from_mol(mol)}

    elif INPUT.isalpha():
        # Assume it is FASTA, if contains only letters
        #print("Input is FASTA")
        mol_unique_ID = 1
        fasta = INPUT
        smi = ''
        #mol_name = 'none'
        inputFile = ''
        outputFile = ''
        l_calc_extn_coeff=True
        l_calc_pI_fasta=True
        l_calc_rdkit_pI=False
        mol_supply_json[mol_unique_ID] = {'mol_name': mol_name, 'mol_obj':None, 'fasta':fasta}

    elif ( (INPUT[0:2] == 'AZ'  and INPUT[2].isdigit()) or (INPUT[0:2] == 'SN'  and INPUT[2].isdigit())  or (INPUT[0:4] == 'MEDI'  and INPUT[4].isdigit())  ):
        # A database ID given
        #print("Input is a database ID")
        mol_unique_ID = 1
        smi = get_smiles_from_dbid(INPUT)
        if len(smi) == 0:
            raise Exception("ERROR: could not convert database ID to smiles. Is it corret ID from supported DBs? Exit.")
            sys.exit(1)
        mol_name = INPUT
        fasta = get_fasta_from_smiles(smi)
        inputFile = ''
        outputFile = ''
        l_calc_extn_coeff=True
        l_calc_pI_fasta=False
        l_calc_rdkit_pI=True
        mol = Chem.MolFromSmiles(smi)
        mol_supply_json[mol_unique_ID] = {'mol_name': mol_name, 'mol_obj':mol, 'fasta':fasta}

    else:
        raise Exception("ERROR: input not recongnized: not smiles, not fasta, not a known databaase ID. Must be a bug. Contact developer. Exit.")
        sys.exit(1)


    #print("mol_name: "+mol_name)
    #print("FASTA: "+fasta)
    #print("SMILES: "+smi)

    dict_out_extn_coeff_fasta = {}
    if l_calc_extn_coeff:
        extn_coeff_options={"seq":'',
                         "inputDict":mol_supply_json,
                         "inputJSON":'',
                         "inputFile":'',
                         "outputFile":'',
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

        pI_fasta_options={"seq":'',
                         "inputDict":mol_supply_json,
                         "inputJSON":'',
                         "inputFile":'',
                         "outputFile":'',
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
                         "lPlot": lPlot,
                         "lIgnoreC": False,
                         "plot_filename":"OUT_titration_curve.png",
                         "l_json":True
                         }

        dict_out_pI_fasta = calc_pI_fasta(pI_fasta_options)

    dict_out_rdkit_pI = {}
    if l_calc_rdkit_pI:

        rdkit_pI_options={'smiles':'',
                        'inputDict':mol_supply_json,
                        'inputJSON':'',
                        'inputFile':inputFile,
                        'outputFile':'',
                        'use_acdlabs':False,
                        'use_pkamatcher':True,
                        'l_print_fragments':args.l_print_fragment_pkas,
                        'l_plot_titration_curve':lPlot,
                        'l_print_pka_set':args.l_print_pka_set,
                        'l_json':True
                         }
                        #'use_dimorphite':True,


        dict_out_rdkit_pI = calc_rdkit_pI(rdkit_pI_options)

    dict_out_peptide_tools_master = {'output_extn_coeff':dict_out_extn_coeff,
                                    'output_pI_fasta':dict_out_pI_fasta,
                                    'output_rdkit_pI':dict_out_rdkit_pI
                                    }

    ### ----------------------------------------------------------------------
    # Output 
    if outputFile == '': # output JSON
        print(json.dumps(dict_out_peptide_tools_master, indent=2))

    else: # output file
        if file_extension != '.fasta':

            #for mi in mol_supply_json.keys():
            mol_list=[]
            for mi in mol_supply_json.keys():
                    mol = mol_supply_json[mi]['mol_obj']
                    
                    if l_calc_rdkit_pI:
                        mol.SetProp('pI mean',"%.2f" % dict_out_rdkit_pI[mi]['pI']['pI mean'])
                        mol.SetProp('pI std',"%.2f" % dict_out_rdkit_pI[mi]['pI']['std'])
                        mol.SetProp('pI interval',' - '.join([ "%.2f" % x for x in dict_out_rdkit_pI[mi]['pI_interval'] ] ))
                        mol.SetProp('pI interval threshold',"%.2f" % dict_out_rdkit_pI[mi]['pI_interval_threshold'])

                    if l_calc_extn_coeff:
                        mol.SetProp('mol_name',dict_out_extn_coeff[mi]['mol_name'])
                        mol.SetProp('Sequence(FASTA)',dict_out_extn_coeff[mi]['fasta'])
                        mol.SetProp('e205(nm)',"%i" % dict_out_extn_coeff[mi]['e205'])
                        mol.SetProp('e214(nm)',"%i" % dict_out_extn_coeff[mi]['e214'])
                        mol.SetProp('e280(nm)',"%i" % dict_out_extn_coeff[mi]['e280'])

                    mol_list.append(mol)

            if out_fext == '.sdf':
                with Chem.SDWriter(outputFile) as sdf_w:
                    for mol in mol_list:
                        sdf_w.write(mol)

            elif out_fext == '.csv':
                with open(outputFile,'w') as csv_f:
                    csv_w = csv.writer(csv_f)
                    count=0
                    for mol in mol_list:
                        props = mol.GetPropsAsDict()

                        count+=1
                        if count == 1:
                            header = ['SMILES'] + list(props.keys())
                            csv_w.writerow(header)

                        row=[Chem.MolToSmiles(mol)]
                        for p in header[1:]:
                            row += [props[p]] 
                        csv_w.writerow(row)
                        
                        

        else:

            dict_list=[]
            for mi in mol_supply_json.keys():
                    fasta = mol_supply_json[mi]['fasta']
                    
                    D = {}

                    if l_calc_pI_fasta:
                        D['pI mean'] = "%.2f" % dict_out_pI_fasta[mi]['pI']['pI mean']
                        D['pI std'] = "%.2f" % dict_out_pI_fasta[mi]['pI']['std']

                    if l_calc_extn_coeff:
                        D['mol_name'] = dict_out_extn_coeff[mi]['mol_name']
                        D['Sequence(FASTA)'] = dict_out_extn_coeff[mi]['fasta']
                        D['e205(nm)'] = "%i" % dict_out_extn_coeff[mi]['e205']
                        D['e214(nm)'] = "%i" % dict_out_extn_coeff[mi]['e214']
                        D['e280(nm)'] = "%i" % dict_out_extn_coeff[mi]['e280']

                    dict_list.append(D)

            if out_fext == '.csv':
                with open(outputFile,'w') as csv_f:
                    csv_w = csv.writer(csv_f)
                    count=0
                    for props in dict_list:

                        count+=1
                        if count == 1:
                            header = list(props.keys())
                            csv_w.writerow(header)

                        row=[]
                        for p in header:
                            row += [props[p]] 
                        csv_w.writerow(row)
                        


        dict_file = {'outputFile':outputFile,'outputInfo':'Number of molecules processed:'+str(len(mol_supply_json.keys()))}
        print(json.dumps(dict_file))
       

 
        
        

