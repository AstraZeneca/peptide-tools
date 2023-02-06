import sys, os, re, string, csv
import argparse
import json
import os.path
import math
from rdkit import Chem
from rdkit.Chem import AllChem, Recap, Descriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize

from numpy import *
from matplotlib.pyplot import *

import subprocess

from pka_sets_fasta import *
from smarts_matcher_aminoacids import *

from itertools import cycle
import bisect 

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.info")



# Turns a dictionary into a class 
class Dict2Class(object): 
    def __init__(self, my_dict): 
        for key in my_dict: 
            setattr(self, key, my_dict[key]) 

def clean_up_output(text):
    txt=text.split('\n')
    text_new=''
    for line in txt: 
        if "mkdir: cannot create directory '/.local': Permission denied" in line: continue
        text_new += line+'\n'
    return text_new

def get_status_output(*args, **kwargs):
    #from subprocess import Popen
    p = subprocess.Popen(*args, **kwargs)
    stdout, stderr = p.communicate()
    return p.returncode, stdout, stderr

def run_exe(exe):
    #if sys.version_info[0] < 3:
    #    from commands import getstatusoutput
    #    status, output = getstatusoutput(exe)
    #else:
    status, stdout, stderr = get_status_output(exe.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output =  stdout.decode("utf-8")
    if status != 0: raise Exception("ERROR: happend while executing "+exe+" : " + output)
    output = clean_up_output(output)
    return output

#def calc_pkas_dimorphite_dl(unknown_fragments):
#   from rdkit import Chem
#   from dimorphite_dl import Protonate
#   from smarts_matcher_nonnaturals_dimorphite import D_dimorphite_dl_type_pka

#   skip_site_names = ['TATA','*Amide']

#   base_pkas=[]
#   acid_pkas=[]
#   diacid_pkas=[]

#   for smiles in unknown_fragments:
#       protonation_sites = Protonate({'smiles':smiles}).get_protonation_sites()

#       for sites in protonation_sites:
#           #(3, 'BOTH', '*Amide') 
#           site_name = sites[2]

#           if site_name in skip_site_names: continue

#           if site_name in D_dimorphite_dl_type_pka.keys():
#               site_data = D_dimorphite_dl_type_pka[site_name]
#               if site_data['type']=='base':
#                       base_pkas.append((site_data['pka'],smiles))
#               if site_data['type']=='acid':
#                       acid_pkas.append((site_data['pka'],smiles))
#               if site_data['type']=='diacid':
#                       diacid_pkas.append(((site_data['pka1'],site_data['pka2']),smiles))
#               
#           else:
#               print('Error:  not known dimorphite site type ' + site_type)
#               sys.exit(1)
#   return (base_pkas,acid_pkas,diacid_pkas)


def calc_pkas_pkamatcher(unknown_fragments):
    from rdkit import Chem
    from smarts_pKaMatcher import list_smarts_pka

    pka_lim_base_1=2
    pka_lim_base_2=15
    pka_lim_acid_1=-5
    pka_lim_acid_2=12

    #skip_smarts_names = ['TATA','*Amide']
    skip_smarts_names = []

    base_pkas=[]
    acid_pkas=[]
    diacid_pkas=[]

    for smiles in unknown_fragments:
        mol = Chem.MolFromSmiles(smiles)
        matched_indices = set()
        for pka_dict_list in list_smarts_pka:
          for smarts_pka in pka_dict_list:
            pat = Chem.MolFromSmarts(smarts_pka['smarts'])
 
            if smarts_pka['name'] in skip_smarts_names: continue

            #print(mol.GetSubstructMatches(pat))
            used_ind_l = []
            used_ind_local = [] 
            for match in mol.GetSubstructMatches(pat):
                used_ind_l += match
                mat_ind = set(match)
                
                available_ind = mat_ind.difference(matched_indices)
               
                if smarts_pka['type']=='base':
                    if match[smarts_pka['ind']-1] in available_ind and match[smarts_pka['ind']-1] not in used_ind_local:
                        pka = smarts_pka['pka']
                        if pka > pka_lim_base_1 and pka < pka_lim_base_2: 
                            base_pkas.append((pka,smiles))

                elif smarts_pka['type']=='acid':
                    if match[smarts_pka['ind']-1] in available_ind and match[smarts_pka['ind']-1] not in used_ind_local:
                        pka = smarts_pka['pka']
                        if pka > pka_lim_acid_1 and pka < pka_lim_acid_2: 
                            acid_pkas.append((pka,smiles))
                else:
                    raise Exception('Error:  not known site type ' + smarts_pka['type'])
                    sys.exit(1)

                # indices of the ionizable centres in the match. used to exclude same ionizable centre be couned multiple times if the same SMARTS matches multiple times. 
                used_ind_local.append(match[smarts_pka['ind']-1])

          matched_indices = matched_indices.union(set(used_ind_l))

    return (base_pkas,acid_pkas,diacid_pkas)
    

def calc_pkas_acdlabs(smi_list):
    # use ACDlabs.
    # perceptabat executable should be in the path.

    pka_lim_base_1=2
    pka_lim_base_2=15
    pka_lim_acid_1=-5
    pka_lim_acid_2=12

    tmpfile='TMP_SMI_FOR_PKA.smi'
    #tmpfile4='TMP_SMI_FOR_PKA.json'
    if os.path.isfile(tmpfile): os.remove(tmpfile)
    with open(tmpfile,'w') as smifile:
        i=0
        for smi in smi_list:
            i+=1
            smifile.write(smi+" tmpname"+str(i)+'\n')

    tmpfile2='TMP_CLAB_OUTPUT.txt'
    if os.path.isfile(tmpfile2): os.remove(tmpfile2)
    #exe = 'perceptabat -TFNAME'+tmpfile2+' -MPKAAPP -TPKA ' + tmpfile
    exe = 'perceptabat -TFNAME'+tmpfile2+' -MPKAAPPGALAS -TPKA ' + tmpfile
    tmpoutput = run_exe(exe)
    if not os.path.isfile(tmpfile2): 
        print("ERROR: no tmpfile generated by acdlabs "+tmpfile2)
        sys.exit(1)


    with open(tmpfile2,'r') as f:
        base_pkas=[]
        acid_pkas=[]
        diacid_pkas=[]
        D={}
        f.readline() # skip first line
        for line in f.readlines():
            ln=line.split()
            mol_idx = int(ln[0])
            if mol_idx not in D.keys(): D[mol_idx]={}
            #if 'ACD_pKa_Apparent' in ln[1]: D[mol_idx]['pka']=float(ln[2])
            if 'ACD_pKa_Apparent' in ln[1]: pka = float(ln[2])
            if 'ACD_pKa_DissType_Apparent' in ln[1]: 
                if ln[2] in ['MB','B']:
                    if pka > pka_lim_base_1 and pka < pka_lim_base_2: 
                        base_pkas.append((pka,smi_list[mol_idx-1]))
                if ln[2] in ['MA','A']:
                    if pka > pka_lim_acid_1 and pka < pka_lim_acid_2: 
                        acid_pkas.append((pka,smi_list[mol_idx-1]))

### Examle of perceptabatch output 
#PKA: Calculate apparent values using classic algorithm
#1	 ID: 1
#1	 ACD_pKa_IonicForm_Apparent: L
#1	 ACD_pKa_Apparent_1: 10.689
#1	 ACD_pKa_Error_Apparent_1: 0.10
#1	 ACD_pKa_DissAtom_Apparent_1: 5
#1	 ACD_pKa_DissType_Apparent_1: MB
#1	 ACD_pKa_Equation_Apparent_1: HL/H+L
#1	 ACD_pKa_All_Apparent_1: pKa(HL/H+L; 5) = 10.69+/-0.10
#2	 ID: 2
#2	 ACD_pKa_IonicForm_Apparent: L
#...
#4	 ID: 4
#4	 ACD_pKa_Caution_Apparent: The structure does not contain ionization centers calculated by current version of ACD/pKa
#5	 ID: 5
#5	 ACD_pKa_IonicForm_Apparent: L
#5	 ACD_pKa_Apparent_1: 10.679
#5	 ACD_pKa_Error_Apparent_1: 0.10
#5	 ACD_pKa_DissAtom_Apparent_1: 1
#5	 ACD_pKa_DissType_Apparent_1: MB
#5	 ACD_pKa_Equation_Apparent_1: HL/H+L
#5	 ACD_pKa_All_Apparent_1: pKa(HL/H+L; 1) = 10.68+/-0.10
#5	 ACD_pKa_Apparent_2: 9.334
#5	 ACD_pKa_Error_Apparent_2: 0.10
#5	 ACD_pKa_DissAtom_Apparent_2: 6
#5	 ACD_pKa_DissType_Apparent_2: B 
#5	 ACD_pKa_Equation_Apparent_2: H2L/H+HL
#5	 ACD_pKa_All_Apparent_2: pKa(H2L/H+HL; 6) = 9.33+/-0.10

    lClean=False
    if lClean:     
            if os.path.isfile(tmpfile): os.remove(tmpfile)
            if os.path.isfile(tmpfile2): os.remove(tmpfile2)

    return (base_pkas,acid_pkas,diacid_pkas)
       


# https://www.rdkit.org/docs/Cookbook.html
# Replace all ionized centers by the corresponsing neutral form. Used to track constanty ionized fragments, like quaternary amines
def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def calc_molecule_constant_charge(net_Qs):
    q = [ v[0] for v in net_Qs]
    if len(q)>0:
        constant_q = float(sum(q))/float(len(q))
    else:
        constant_q = 0.0

    return constant_q
    

def standardize_molecule(mol):
    clean_mol = rdMolStandardize.Cleanup(mol)
    m = neutralize_atoms(clean_mol)
    return m
    

def calc_net_Qs(smi_list):
    net_Qs = []
    for smi in smi_list:
        mol = Chem.MolFromSmiles(smi)
        m = standardize_molecule(mol)

        # positively charged N, but not connected to any negatively charged atom. To avoid azido and nitro groups being counted.
        pattern = Chem.MolFromSmarts("[#7+;!H1;!H2;!H3;!H4!$([#7+]~[*-])]")
        at_matches = m.GetSubstructMatches(pattern)

        at_matches_list = [y[0] for y in at_matches]

        for v in at_matches_list:
            net_Qs.append( (1,smi) )
    
    return net_Qs
        
        

# calcualtes pKa for the list of smiles. 
def calc_pkas(smi_list, use_acdlabs=False, use_pkamatcher=True):

    if use_acdlabs and use_pkamatcher: raise Exception('Error: requested to use both ACDlabs and pKaMatcher for pka calculation. Should be only one. ')
    if not use_acdlabs and not use_pkamatcher: raise Exception('Error: requested to use none of ACDlabs or pKaMatcher for pka calculation. Should be at least one. ')

    if use_acdlabs:
        base_pkas,acid_pkas,diacid_pkas = calc_pkas_acdlabs(smi_list)

    if use_pkamatcher:
        base_pkas,acid_pkas,diacid_pkas = calc_pkas_pkamatcher(smi_list)

    #if use_opera:
    #    base_pkas,acid_pkas,diacid_pkas = calc_pkas_opera(smi_list)

    #if use_dimorphite:
    #    base_pkas,acid_pkas,diacid_pkas = calc_pkas_dimorphite_dl(smi_list)

    net_Qs = calc_net_Qs(smi_list) 
    
    return (base_pkas,acid_pkas,diacid_pkas,net_Qs)


        
def calculateBasicCharge(pH, pKa):
        return 1 / (1 + 10**(pH - pKa))

def calculateAcidicCharge(pH, pKa):
        return -1 / (1 + 10**(pKa - pH))

def calculateDiacidCharge(pH, pKa1, pKa2):
        Ka1=10**(-pKa1)
        Ka2=10**(-pKa2)
        H=10**(-pH)
        f1 = (H*Ka1)/(H**2+H*Ka1+Ka1*Ka2)  # fraction of [AH-]
        f2 = f1 * Ka2 / H                  # fraction of [A2-]
        return -2*f2 + (-1)*f1     # average charge of phosphate group


def calculateMolCharge(base_pkas, acid_pkas, diacid_pkas, pH, constant_q=0):
    charge = constant_q
    for pka in base_pkas:
        charge += calculateBasicCharge(pH, pka)

    for pka in acid_pkas:
        charge += calculateAcidicCharge(pH, pka)

    for pkas in diacid_pkas:
        charge += calculateDiacidCharge(pH, pkas)

    #print(pH,charge)
    return charge


# Define pH span tocalcualte itration curve and where to search for pI.
def define_pH_span():
    pH_llim=-1
    pH_hlim=15
    return [pH_llim,pH_hlim]


def calculateIsoelectricPoint(base_pkas, acid_pkas, diacid_pkas, constant_q=0):   
    tolerance=0.01
    charge_tol=0.05
    na=len(acid_pkas)+len(diacid_pkas)
    nb=len(base_pkas)
    
    pH_lim = define_pH_span()
    lower_pH = pH_lim[0] 
    higher_pH = pH_lim[1] 

    while True:
        mid_pH = 0.5 * (higher_pH + lower_pH)
        charge = calculateMolCharge(base_pkas, acid_pkas, diacid_pkas, mid_pH, constant_q=constant_q)
        
        if na == 0 and nb != 0:
            #print "---!Warning: no acidic ionizable groups, only basic groups present in the sequence. pI is not defined and thus won't be calculated. However, you can still plot the titration curve. Continue."
            refcharge = charge_tol * nb

        elif nb == 0 and na != 0:
            #print "---!Warning: no basic ionizable groups, only acidic groups present in the sequence. pI is not defined and thus won't be calculated. However, you can still plot the titration curve. Continue."
            refcharge = -charge_tol * na

        else:
            refcharge = 0.0


        if charge > refcharge + tolerance:
            lower_pH = mid_pH
        elif charge < refcharge - tolerance:
            higher_pH = mid_pH
        else:
            return mid_pH
            
        if mid_pH <= pH_lim[0]:
            return pH_lim[0]
        elif mid_pH >= pH_lim[1]:
            return pH_lim[1]
            

def CalcChargepHCurve(base_pkas, acid_pkas, diacid_pkas, constant_q=0):
    pH_lim = define_pH_span()
    dpH=0.1
    pH_a=arange(pH_lim[0],pH_lim[1]+dpH,dpH)
    Q_a=pH_a*0.0    
    for i in range(len(pH_a)):
        Q = calculateMolCharge(base_pkas, acid_pkas, diacid_pkas, pH_a[i],constant_q=constant_q)
        Q_a[i]=Q
    pH_Q = np.vstack((pH_a,Q_a))
    return pH_Q



def fragmentMol(mol, breakableBonds, atomNum):
    # Function from breaking and capping
    readwritemol = Chem.RWMol(mol)
    for bond in breakableBonds: 
        atm1 = bond[0]
        atm2 = bond[1]
        readwritemol.RemoveBond(atm1, atm2)
        amineC = readwritemol.AddAtom(Chem.Atom(6))
        amineCO = readwritemol.AddAtom(Chem.Atom(8))
        amineCOC = readwritemol.AddAtom(Chem.Atom(6))
        readwritemol.AddBond(atm1, amineC, Chem.BondType.SINGLE) 
        readwritemol.AddBond(amineC, amineCO, Chem.BondType.DOUBLE) 
        readwritemol.AddBond(amineC, amineCOC, Chem.BondType.SINGLE) 
        newatom2 = readwritemol.AddAtom(Chem.Atom(atomNum))
        readwritemol.AddBond(atm2, newatom2, Chem.BondType.SINGLE) 
        fragmentedSmiles = Chem.MolToSmiles(readwritemol)
    return fragmentedSmiles

def break_amide_bonds_and_cap(mol):
    ### SMARTS pattern
    amideSMARTS = '[NX3,NX4;H0,H1][CX3](=[OX1])' # secondary and tertiary amide bonds
    amidePattern = Chem.MolFromSmarts(amideSMARTS)
    atomNum = 6

    breakableBonds = []
    #molname=mol.GetProp("_Name")
    smiles = Chem.MolToSmiles(mol)
    if mol.HasSubstructMatch(amidePattern):
        #smiles = Chem.MolToSmiles(mol)
        #results.write("%s," % molname)
        atomIDs = mol.GetSubstructMatches(amidePattern)
        for bond in atomIDs:
            breakableBonds.append((bond[0],bond[1]))
        fragmentedSmiles = fragmentMol(mol, breakableBonds, atomNum)
    else:
        fragmentedSmiles = smiles

    fragmentedSmilesList = fragmentedSmiles.split(".")
    numFrags = len(fragmentedSmilesList)
    #print("%s fragmented into %s pieces" % (molname,numFrags))
    return fragmentedSmilesList


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

def print_output_prop_dict(prop_dict,prop,l_print_pka_set=False):
    global tit
    lj=12
    #keys = prop_dict.keys()
    keys = list(prop_dict.keys())
    keys.remove('std'); keys.insert(0, 'std')
    keys.remove('err'); keys.insert(0, 'err')
    keys.remove(prop + ' mean'); keys.insert(0, prop+' mean')
    tit="sequence"
    print(" ")
    print("======================================================================================================================================================")
    print(prop)
    print( "---------------------------------")
    for k in keys:
        p = prop_dict[k]
        print(k.rjust(lj)  + "  " +  str(round(p,2)).ljust(lj) )
    print(" ")

    if l_print_pka_set: print_pka_set()

    return


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
            Q_M = column_stack([Q_M,pH_Q[:,1]])

    plot(pH,pH*0,'k-')
    plot([7,7],[np.min(Q_M),np.max(Q_M)],'k-')
    #xlim([np.min(pH),np.max(pH)])
    xlim([2,12])
    ylim([np.min(Q_M),np.max(Q_M)])

    legend(loc="center right", bbox_to_anchor=[1.1, 0.5],ncol=1, shadow=True, fontsize=10).get_frame().set_alpha(1)
    ylabel('peptide charge')
    xlabel('pH')
	
    minorticks_on()
    grid(True)

    #show()
    savefig(figFileName)
    return


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
            if not mol.HasProp('_Name'): mol.SetProp('_Name','tmpname'+str(mol_unique_ID))
            mol_supply_json[mol_unique_ID] = {'mol_name':mol.GetProp('_Name'), 'mol_obj':mol, 'fasta':''}

        return mol_supply_json



def calc_pIChemiSt(options={'smiles':'','inputDict':{},'inputJSON':'','inputFile':'','outputFile':'','use_acdlabs':False,'use_pkamatcher':True,'l_print_fragments':False,'l_plot_titration_curve':False,'l_print_pka_set':False,'l_json':False}):

    global mol_supply_json

    args = Dict2Class(options)

    # Get options
    if len(args.smiles)!=0:
        # assume smiles input
        mol_unique_ID = 1
        mol = Chem.MolFromSmiles(args.smiles)
        mol_supply_json={}
        mol_supply_json[mol_unique_ID] = {'mol_name': 'none', 'mol_obj':mol, 'fasta':'none'}
    elif len(args.inputFile)!=0:
        # Assume filename as input
        inputFile = args.inputFile
        mol_supply_json = read_structure_file(inputFile)
    elif len(args.inputJSON)!=0:
        # Assume molecule JSON supply as input
        mol_supply_json = json.loads(args.inputJSON)
    elif args.inputDict: # if not an empty dictionary
        # Assume Dict olecule supply as input
        mol_supply_json = args.inputDict
    else:
        raise Exception('Error: either smiles or input file *.smi, sdf, or Json string should be given. Exit. ')
        sys.exit(1)

    # Run calculations
    dict_output_pIChemiSt={}

    for molid_ind in mol_supply_json.keys():
        mol_name = mol_supply_json[molid_ind]['mol_name']    
        mol = mol_supply_json[molid_ind]['mol_obj']    
        mol = standardize_molecule(mol)

        # break amide bonds
        frags_smi_list = break_amide_bonds_and_cap(mol)

        # match known amino-acids with defined pKas
        unknown_frags,base_pkas_fasta,acid_pkas_fasta,diacid_pkas_fasta = get_pkas_for_known_AAs(frags_smi_list)
        #print('UNKNOWN_FRAGMENTS: '+'   '.join(unknown_frags))

        # caclulate pKas for unknown fragmets
        if len(unknown_frags) > 0: base_pkas_calc,acid_pkas_calc,diacid_pkas_calc,net_Qs = calc_pkas(unknown_frags,use_acdlabs=args.use_acdlabs,use_pkamatcher=args.use_pkamatcher)
        else: base_pkas_calc,acid_pkas_calc,diacid_pkas_calc,net_Qs = [],[],[],[]

        # loop over all pKa sets
        seq_dict={}
        pI_dict={}
        Q_dict={}
        pH_Q_dict={}
        for pKaset in pKa_sets_to_use:

            # merge fasta and calcualted pkas
            base_pkas = base_pkas_fasta[pKaset] + base_pkas_calc
            acid_pkas = acid_pkas_fasta[pKaset] + acid_pkas_calc
            diacid_pkas = diacid_pkas_fasta[pKaset] + diacid_pkas_calc

            all_base_pkas=[]
            all_acid_pkas=[]
            all_diacid_pkas=[]
            if len(base_pkas) != 0: all_base_pkas,all_base_pkas_smi = zip(*base_pkas) 
            else: all_base_pkas,all_base_pkas_smi = [],[]
            if len(acid_pkas) != 0: all_acid_pkas,all_acid_pkas_smi = zip(*acid_pkas) 
            else: all_acid_pkas,all_acid_pkas_smi = [],[]
            if len(diacid_pkas) != 0: all_diacid_pkas,all_diacid_pkas_smi = zip(*diacid_pkas) 
            else: all_diacid_pkas,all_diacid_pkas_smi = [],[]

            # calculate isoelectric point
            molecule_constant_charge = calc_molecule_constant_charge(net_Qs)
            pI = calculateIsoelectricPoint(all_base_pkas, all_acid_pkas, all_diacid_pkas, constant_q = molecule_constant_charge)
            pH_Q = CalcChargepHCurve(all_base_pkas, all_acid_pkas, all_diacid_pkas, constant_q = molecule_constant_charge)
            pH_Q = pH_Q.T
            #print( "pI ACDlabs %6.3f" % (pI) )

            # calculate net charge at pH 7.4
            Q = calculateMolCharge(all_base_pkas, all_acid_pkas, all_diacid_pkas, 7.4, constant_q = molecule_constant_charge)
            #print( "Q at pH7.4 ACDlabs %4.1f" % (Q) )



            pI_dict[pKaset] = pI
            Q_dict[pKaset] = Q
            pH_Q_dict[pKaset] = pH_Q



        # calcualte isoelectric interval and reset undefined pI values 
        int_tr = 0.2    # TODO define it elsewhere 
        
        pH_lim = define_pH_span()

        interval_low_l = []
        interval_high_l = []
        for pKaset in pKa_sets_to_use:
            pH_Q = pH_Q_dict[pKaset]
            Q=pH_Q[:,1]
            pH=pH_Q[:,0]
            pH_int = ( pH[(Q>-int_tr) & (Q<int_tr)] )
            
            
            # isoelectric interval - pH range where the charge is within the given threshold. If molecule permanently has a charge the interval is not defined and NaN are provided. 

            #print(pH_int[0],pH_int[-1],pH_lim[0],pH_lim[1])
            if len(pH_int) > 0:

                # case when pI is not defined
                if round(pH_int[0],4) == round(pH_lim[0],4) and round(pH_int[-1],4) == round(pH_lim[1],4):
                    pI_dict[pKaset] = float('NaN')    

                interval_low_l.append(pH_int[0])
                interval_high_l.append(pH_int[-1])
            
            else:
                # case when pI is not defined
                pI_dict[pKaset] = float('NaN')    
                interval_low_l.append(pH[0])
                interval_high_l.append(pH[-1])
           
        if len(interval_low_l)>0:
            interval_low = mean(interval_low_l)
        else:
            interval_low = float('NaN')
            
        if len(interval_high_l)>0:
            interval_high = mean(interval_high_l)
        else:
            interval_high = float('NaN')
            
        interval = (interval_low, interval_high) 




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

        # plot titration curve
        if args.l_plot_titration_curve:
            figFileName = "OUT_titration_curve_pIChemiSt.png"
            plot_titration_curve(pH_Q_dict,figFileName)
        else:
            figFileName = ""
        
        # output dict for given molecule 
        dict_output_pIChemiSt[molid_ind]={'mol_name':mol_name,
                                    'pI':pI_dict,
                                    'QpH7':Q_dict,
                                    'pI_interval':interval,
                                    'plot_filename':figFileName,
                                    'pI_interval_threshold':int_tr
                                    }
        
        # define pKaset for reporting pKa of individual amino acids and fragments
        pKaset='IPC2_peptide'
        dict_output_pIChemiSt[molid_ind].update({'pKa_set':pKaset })

        
        if args.l_print_fragments:
            dict_output_pIChemiSt[molid_ind].update({
                                    'base_pkas_fasta':base_pkas_fasta,
                                    'acid_pkas_fasta':acid_pkas_fasta,
                                    'base_pkas_calc':base_pkas_calc,
                                    'acid_pkas_calc':acid_pkas_calc,
                                    'constant_Qs_calc':net_Qs
                                    })

                                    # Dicacids pkas are included as single apparent ionizaitions. No need to include diacids. 
                                    #'diacid_pkas_calc':diacid_pkas_calc,
                                    #'diacid_pkas_fasta':diacid_pkas_fasta,

    return dict_output_pIChemiSt



def print_output(dict_output_pIChemiSt,args):

    for molid_ind in dict_output_pIChemiSt.keys():
    
        molid = dict_output_pIChemiSt[molid_ind]

        print_output_prop_dict(dict_output_pIChemiSt[molid_ind]['pI'],'pI',l_print_pka_set=args.l_print_pka_set)
        print_output_prop_dict(dict_output_pIChemiSt[molid_ind]['QpH7'],'Q at pH7.4',l_print_pka_set=False)

        if args.use_acdlabs: predition_tool = 'ACDlabs'
        elif args.use_pkamatcher: predition_tool = 'pKaMatcher'

        int_tr = dict_output_pIChemiSt[molid_ind]['pI_interval_threshold']
        pKaset = dict_output_pIChemiSt[molid_ind]['pKa_set']

        print(" ")
        #print("pH interval with charge between %4.1f and %4.1f for pKa set: %s and prediction tool: %s" % (-int_tr,int_tr,pKaset,predition_tool) )
        print("pH interval with charge between %4.1f and %4.1f and prediction tool: %s" % (-int_tr,int_tr,predition_tool) )
        print("%4.1f - %4.1f" % (dict_output_pIChemiSt[molid_ind]['pI_interval'][0],dict_output_pIChemiSt[molid_ind]['pI_interval'][1]))

        if args.l_print_fragments:
            base_pkas_fasta = dict_output_pIChemiSt[molid_ind]['base_pkas_fasta']
            acid_pkas_fasta = dict_output_pIChemiSt[molid_ind]['acid_pkas_fasta']
            #diacid_pkas_fasta = dict_output_pIChemiSt[molid_ind]['diacid_pkas_fasta']
            base_pkas_calc = dict_output_pIChemiSt[molid_ind]['base_pkas_calc']
            acid_pkas_calc = dict_output_pIChemiSt[molid_ind]['acid_pkas_calc']
            #diacid_pkas_calc = dict_output_pIChemiSt[molid_ind]['diacid_pkas_calc']
            constant_Qs_calc = dict_output_pIChemiSt[molid_ind]['constant_Qs_calc']

            # merge fasta and calcualted pkas
            base_pkas = base_pkas_fasta[pKaset] + base_pkas_calc
            acid_pkas = acid_pkas_fasta[pKaset] + acid_pkas_calc
            #diacid_pkas = diacid_pkas_fasta[pKaset] + diacid_pkas_calc
            all_base_pkas=[]
            all_acid_pkas=[]
            all_diacid_pkas=[]
            if len(base_pkas) != 0: all_base_pkas,all_base_pkas_smi = zip(*base_pkas) 
            else: all_base_pkas,all_base_pkas_smi = [],[]
            if len(acid_pkas) != 0: all_acid_pkas,all_acid_pkas_smi = zip(*acid_pkas) 
            else: all_acid_pkas,all_acid_pkas_smi = [],[]
            #if len(diacid_pkas) != 0: all_diacid_pkas,all_diacid_pkas_smi = zip(*diacid_pkas) 
            #else: all_diacid_pkas,all_diacid_pkas_smi = [],[]
        
            print(" ")
            print("List of calculated BASE pKa's with the corresponding fragments")
            for pkas,smi in zip(all_base_pkas,all_base_pkas_smi):
                s_pkas = ["%4.1f"%(pkas)]
                print("smiles or AA, base pKa : %-15s %s" % (smi,' '.join(s_pkas)))

            print(" ")
            print("List of calculated ACID pKa's with the corresponding fragments")
            for pkas,smi in zip(all_acid_pkas,all_acid_pkas_smi):
                s_pkas = ["%4.1f"%(pkas)]
                print("smiles or AA, acid pKa : %-15s %s" % (smi,' '.join(s_pkas)))

#            print(" ")
#            print("List of calculated DIACID pKa's with the corresponding fragments")
#            for pkas,smi in zip(all_diacid_pkas,all_diacid_pkas_smi):
#                s_pkas = ["%4.1f  %4.1f"%(pkas[0],pkas[1])]
#                print("smiles or AA, diacid pka : %-15s %s" % (smi,' '.join(s_pkas)))

            print(" ")
            print("List of constantly ionized fragments")
            for v in constant_Qs_calc:
                pkas = v[0]
                smi = v[1]
                s_pkas = ["%4.1f"%(pkas)]
                print("smiles, charge : %-15s %s" % (smi,' '.join(s_pkas)))

    return







def args_parser():
    parser = argparse.ArgumentParser(description="Script caclultes isoelectric point of a molecules by cutting all amide bonds, retreiving stored pka values for known AAs, predicting pKas of unknown fragemnts using pKaMatcher or ACDlabs, and calculating Henderson-Hasselbalch equation.")
    parser.add_argument("-j", dest="inputJSON", help="input molecule supply in JSON format",default='')
    parser.add_argument("-i", dest="inputFile", help="input file with molecule structure. smi or sdf",default='')
    parser.add_argument("-s", dest="smiles", help="input smiles. if used then single smi is assumed and fasta returned in stdout. Input filenames are ignored",default='')
    parser.add_argument("-o", dest="outputFile", help="output file with molecule structure. csv or sdf",default='')

    parser.add_argument("--plot_titration_curve",default=False, action='store_true',dest="l_plot_titration_curve", help="Plot titration curve and store it in OUT_titration_curve_pIChemiSt.png file.")
    parser.add_argument("--print_fragment_pkas",default=False, action='store_true',dest="l_print_fragments", help="Print out fragments with corresponding pKas used in pI calcution.")
    parser.add_argument("--print_pka_set",default=False, action='store_true',dest="l_print_pka_set", help="Print out stored pka sets explicitly.")
    parser.add_argument("--use_acdlabs",default=False, action='store_true',dest="use_acdlabs", help="Use acdlabs for pka prediction of unknown fragments.")
    parser.add_argument("--use_pkamatcher",default=False, action='store_true',dest="use_pkamatcher", help="Use pKaMatcher for pka prediction of unknown fragments.")
    parser.add_argument("--json",default=False, action='store_true',dest="l_json", help="Output in JSON format.")

    args = parser.parse_args()
    options = args.__dict__

    return args,options



if __name__ == "__main__":

    args,options = args_parser()

    dict_output_pIChemiSt = calc_pIChemiSt(options)

    ### ----------------------------------------------------------------------
    # Output 
    if args.outputFile == '': # output plain text
        if args.l_json:
            print(json.dumps(dict_output_pIChemiSt, indent=2))
        else:
            print_output(dict_output_pIChemiSt,args)    

    else: # output file

        known_out_file_types = ['.sdf','.csv']
        filename, out_fext = os.path.splitext(args.outputFile)
        if out_fext not in known_out_file_types:
            raise Exception('Error! Output file extention not in supported file types:'+str(known_file_types))
            sys.exit(1)

        mol_list=[]
        for mi in mol_supply_json.keys():
            mol = mol_supply_json[mi]['mol_obj']
            mol.SetProp('pI mean',"%.2f" % dict_output_pIChemiSt[mi]['pI']['pI mean'])
            mol.SetProp('pI std',"%.2f" % dict_output_pIChemiSt[mi]['pI']['std'])
            mol.SetProp('pI interval',' - '.join([ "%.2f" % x for x in dict_output_pIChemiSt[mi]['pI_interval'] ] ))
            mol.SetProp('pI interval threshold',"%.2f" % dict_output_pIChemiSt[mi]['pI_interval_threshold'])

            mol_list.append(mol)

        if out_fext == '.sdf':
            with Chem.SDWriter(args.outputFile) as sdf_w:
                for mol in mol_list:
                    sdf_w.write(mol)

        elif out_fext == '.csv':
            with open(args.outputFile,'w') as csv_f:
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
                        
        # print info 
        dict_file = {'outputFile':args.outputFile,'outputInfo':'Number of molecules processed:'+str(len(dict_output_pIChemiSt.keys()))}
        print(json.dumps(dict_file))
   

















 
