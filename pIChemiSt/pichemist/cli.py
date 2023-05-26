import sys, os, csv
import argparse
import json
import os.path

# TODO: Explicit imports 
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from rdkit import Chem
from rdkit import RDLogger

from pichemist.config import PKA_SETS_NAMES
from pichemist.api import calc_pkas
from pichemist.io import generate_input
from pichemist.molecule import MolStandardiser
from pichemist.molecule import PeptideCutter
from pichemist.fasta.matcher import get_aa_pkas_for_list
from pichemist.model import PKaMethod
from pichemist.model import InputFormat
from pichemist.model import OutputFormat
from pichemist.model import MODELS

from pichemist.stats import mean
from pichemist.stats import stddev
from pichemist.stats import stderr
from pichemist.utils import get_logger

# Configure logging
log = get_logger(__name__)
RDLogger.DisableLog("rdApp.info")


def calc_molecule_constant_charge(net_Qs):
    q = [ v[0] for v in net_Qs]
    if len(q)>0:
        constant_q = float(sum(q))/float(len(q))
    else:
        constant_q = 0.0

    return constant_q


        
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
    pH_a = np.arange(pH_lim[0],pH_lim[1]+dpH,dpH)
    Q_a=pH_a*0.0    
    for i in range(len(pH_a)):
        Q = calculateMolCharge(base_pkas, acid_pkas, diacid_pkas, pH_a[i],constant_q=constant_q)
        Q_a[i]=Q
    pH_Q = np.vstack((pH_a,Q_a))
    return pH_Q


def print_output_prop_dict(prop_dict, prop, print_pka_set=False):
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

    if print_pka_set: print_pka_set()   # ANDREY: print_pka_set() does not exist anywhere

    return


### PLOT titration curve
def _plot_titration_curve(pH_Q_dict,figFileName):
    plt.matplotlib.rcParams.update({'font.size': 16})
    lines = ["-","--","-.",":"]
    w1=4.0 ; w2=3.0 ; w3=2.0 ; w4=1.0
    linew = [w1,w1, w2,w2, w3,3, w4,w4]
    linecycler = cycle(lines)
    linewcycler = cycle(linew)

    plt.figure(figsize=(8,6))
    i=0
    for pKaset in PKA_SETS_NAMES:
        i+=1
        pH_Q = pH_Q_dict[pKaset] 
        l = plt.plot(pH_Q[:,0],pH_Q[:,1],next(linecycler),label=pKaset,linewidth=next(linewcycler)) 
        if pKaset == 'IPC2_peptide': 
            plt.setp(l,linewidth=8,linestyle='-',color='k')

        # Store data for output
        if i==1: 
            pH = pH_Q[:,0]
            Q_M = pH_Q[:,1]
        else:
            Q_M = np.column_stack([Q_M,pH_Q[:,1]])

    plt.plot(pH,pH*0,'k-')
    plt.plot([7,7],[np.min(Q_M),np.max(Q_M)],'k-')
    #xlim([np.min(pH),np.max(pH)])
    plt.xlim([2,12])
    plt.ylim([np.min(Q_M),np.max(Q_M)])

    plt.legend(loc="center right", bbox_to_anchor=[1.1, 0.5],ncol=1, shadow=True, fontsize=10).get_frame().set_alpha(1)
    plt.ylabel('peptide charge')
    plt.xlabel('pH')
	
    plt.minorticks_on()
    plt.grid(True)

    #show()
    plt.savefig(figFileName)
    return

# TODO: use input and type to determine input
# TODO: use method instead of use_acdlabs and use_pkamatcher
def calc_pichemist(input_dict, method,
                   plot_titration_curve=False,
                   print_fragments=False,
                   print_pka_set=False):

    # Run calculations
    dict_output={}
    for mol_idx in input_dict.keys():
        # Prepare molecule and break into fragments
        mol_name = input_dict[mol_idx]['mol_name']
        mol = input_dict[mol_idx]['mol_obj']    
        mol = MolStandardiser().standardise_molecule(mol)
        frags_smi_list = PeptideCutter().break_amide_bonds_and_cap(mol)
        # ANDREY: calc_net_Qs wants standardise_molecule() for each mol in the list
        # should that be done here instead or maybe not at all since we do for the
        # full peptide?

        # TODO: Match and calc pKas should be an API that contains the following
            # TODO: Match AA pKas from SMILES list should also be an API
            # TODO: Predict pKas from SMILES list should also be an API

        # Match known pKas from FASTA definitions
        unknown_frags, base_pkas_fasta, acid_pkas_fasta, diacid_pkas_fasta = get_aa_pkas_for_list(frags_smi_list)

        # caclulate pKas for unknown fragmets
        base_pkas_calc, acid_pkas_calc, diacid_pkas_calc, net_Qs = list(), list(), list(), list()
        if len(unknown_frags) > 0:
            base_pkas_calc, acid_pkas_calc, diacid_pkas_calc, net_Qs = calc_pkas(unknown_frags,
                                                                                 method=method)

        # loop over all pKa sets
        pI_dict={}
        Q_dict={}
        pH_Q_dict={}
        for pKaset in PKA_SETS_NAMES:

            # merge fasta and calcualted pkas
            base_pkas = base_pkas_fasta[pKaset] + base_pkas_calc
            acid_pkas = acid_pkas_fasta[pKaset] + acid_pkas_calc
            diacid_pkas = diacid_pkas_fasta[pKaset] + diacid_pkas_calc

            all_base_pkas=[]
            all_acid_pkas=[]
            all_diacid_pkas=[]
            if len(base_pkas) != 0:
                all_base_pkas, all_base_pkas_smi = zip(*base_pkas) 
            else:
                all_base_pkas, all_base_pkas_smi = [],[]
            if len(acid_pkas) != 0:
                all_acid_pkas, all_acid_pkas_smi = zip(*acid_pkas) 
            else:
                all_acid_pkas, all_acid_pkas_smi = [],[]
            if len(diacid_pkas) != 0:
                all_diacid_pkas, all_diacid_pkas_smi = zip(*diacid_pkas) 
            else:
                all_diacid_pkas, all_diacid_pkas_smi = [],[]

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
        for pKaset in PKA_SETS_NAMES:
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


        # print isoelectric interval
        int_tr = 0.2    # TODO define it elsewhere 

        interval_low_l = []
        interval_high_l = []
        for pKaset in PKA_SETS_NAMES:
            pH_Q = pH_Q_dict[pKaset]
            Q=pH_Q[:,1]
            pH=pH_Q[:,0]
            pH_int = ( pH[(Q>-int_tr) & (Q<int_tr)] )
            
            # isoelectric interval - pH range where the charge is within the given threshold. If molecule permanently has a charge the interval is not defined and NaN are provided. 
            if len(pH_int) > 1:
                #interval = (pH_int[0], pH_int[-1])
                interval_low_l.append(pH_int[0])
                interval_high_l.append(pH_int[-1])
            #else:
            #    interval_low.append(float('NaN'))
            #    interval_low.append(float('NaN'))
           
        if len(interval_low_l)>0:
            interval_low = mean(interval_low_l)
        else:
            interval_low = float('NaN')
            
        if len(interval_high_l)>0:
            interval_high = mean(interval_high_l)
        else:
            interval_high = float('NaN')
            
        interval = (interval_low, interval_high) 

        # plot titration curve
        if plot_titration_curve:
            figFileName = "OUT_titration_curve_pIChemiSt.png"
            _plot_titration_curve(pH_Q_dict,figFileName)
        else:
            figFileName = ""
        
        # output dict for given molecule 
        dict_output[mol_idx]={'mol_name':mol_name,
                            'pI':pI_dict,
                            'QpH7':Q_dict,
                            'pI_interval':interval,
                            'plot_filename':figFileName,
                            'pI_interval_threshold':int_tr
                            }
        
        # define pKaset for reporting pKa of individual amino acids and fragments
        pKaset='IPC2_peptide'
        dict_output[mol_idx].update({'pKa_set':pKaset })

        
        if print_fragments:
            dict_output[mol_idx].update({
                                    'base_pkas_fasta':base_pkas_fasta,
                                    'acid_pkas_fasta':acid_pkas_fasta,
                                    'base_pkas_calc':base_pkas_calc,
                                    'acid_pkas_calc':acid_pkas_calc,
                                    'constant_Qs_calc':net_Qs
                                    })

                                    # Dicacids pkas are included as single apparent ionizaitions. No need to include diacids. 
                                    #'diacid_pkas_calc':diacid_pkas_calc,
                                    #'diacid_pkas_fasta':diacid_pkas_fasta,

    return dict_output



def print_output(dict_output, args):

    for mol_idx in dict_output.keys():
    
        molid = dict_output[mol_idx]

        print_output_prop_dict(dict_output[mol_idx]['pI'],'pI', print_pka_set=args.print_pka_set)
        print_output_prop_dict(dict_output[mol_idx]['QpH7'],'Q at pH7.4', print_pka_set=False)

        if args.method == "acd":
            predition_tool = 'ACDlabs'
        if args.method == "pkamatcher":
            predition_tool = 'pKaMatcher'

        int_tr = dict_output[mol_idx]['pI_interval_threshold']
        pKaset = dict_output[mol_idx]['pKa_set']

        print(" ")
        #print("pH interval with charge between %4.1f and %4.1f for pKa set: %s and prediction tool: %s" % (-int_tr,int_tr,pKaset,predition_tool) )
        print("pH interval with charge between %4.1f and %4.1f and prediction tool: %s" % (-int_tr,int_tr,predition_tool) )
        print("%4.1f - %4.1f" % (dict_output[mol_idx]['pI_interval'][0],dict_output[mol_idx]['pI_interval'][1]))

        if args.print_fragments:
            base_pkas_fasta = dict_output[mol_idx]['base_pkas_fasta']
            acid_pkas_fasta = dict_output[mol_idx]['acid_pkas_fasta']
            #diacid_pkas_fasta = dict_output[mol_idx]['diacid_pkas_fasta']
            base_pkas_calc = dict_output[mol_idx]['base_pkas_calc']
            acid_pkas_calc = dict_output[mol_idx]['acid_pkas_calc']
            #diacid_pkas_calc = dict_output[mol_idx]['diacid_pkas_calc']
            constant_Qs_calc = dict_output[mol_idx]['constant_Qs_calc']

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



__doc__ = "Calculates the isoelectric point (pI) of a molecules by cutting its amide bonds, retreiving the pKa values for known AAs, predicting pKas of unknown fragments using pKaMatcher or ACDlabs, and finally calculating the pI using the Henderson-Hasselbalch equation."



def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", dest="input",
                        help="Input to the algorithm. Depends on the format used",
                        default=None)
    parser.add_argument("-if", dest="input_format",
                        help="Format of the input", choices=MODELS[InputFormat],
                        default=InputFormat.SMILES_FILE)
    parser.add_argument("-o", dest="output_file",
                        help="Output file for the results of the calculation",
                        default=None)
    parser.add_argument("-of", dest="output_format",
                        help="Format of the input", choices=MODELS[OutputFormat],
                        default=OutputFormat.JSON)
    parser.add_argument("--plot_titration_curve", default=False,
                        action='store_true', dest="plot_titration_curve",
                        help="TODO:")
    parser.add_argument("--print_fragment_pkas", default=False,
                        action='store_true', dest="print_fragment_pkas",
                        help="TODO: Print out fragments with corresponding pKas used in pI calcution.")
    parser.add_argument("--print_pka_set", default=False,
                        action='store_true', dest="print_pka_set",
                        help="TODO: Print out stored pka sets explicitly.")
    parser.add_argument("--method",
                        choices=MODELS[PKaMethod],
                        default=PKaMethod.PKA_MATCHER,
                        help="Method for the prediction of pKas of unknown fragments")
    parser.add_argument("--json",default=False, action='store_true',dest="l_json", help="Output in JSON format.")

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = arg_parser()
    input_dict = generate_input(args.input_format, args.input)
    dict_output = calc_pichemist(input_dict, args.method,
                                 args.plot_titration_curve,
                                 args.print_fragment_pkas)
    print(dict_output)
    exit(0)

    ### ----------------------------------------------------------------------
    # Output
    if args.output_file is None: # output plain text
        if args.output_format == OutputFormat.JSON:
            print(json.dumps(dict_output, indent=2))
        else:
            print_output(dict_output,args)    

    else: # output file

        known_out_file_types = ['.sdf','.csv']
        filename, out_fext = os.path.splitext(args.output_file)
        if out_fext not in known_out_file_types:
            raise Exception('Error! Output file extention not in supported file types:'+str(known_file_types))
            sys.exit(1)

        mol_list=[]
        for mi in input_dict.keys():
            mol = input_dict[mi]['mol_obj']
            mol.SetProp('pI mean',"%.2f" % dict_output[mi]['pI']['pI mean'])
            mol.SetProp('pI std',"%.2f" % dict_output[mi]['pI']['std'])
            mol.SetProp('pI interval',' - '.join([ "%.2f" % x for x in dict_output[mi]['pI_interval'] ] ))
            mol.SetProp('pI interval threshold',"%.2f" % dict_output[mi]['pI_interval_threshold'])

            mol_list.append(mol)

        if out_fext == '.sdf':
            with Chem.SDWriter(args.output_file) as sdf_w:
                for mol in mol_list:
                    sdf_w.write(mol)

        elif out_fext == '.csv':
            with open(args.output_file,'w') as csv_f:
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
        dict_file = {'output_file':args.output_file,'outputInfo':'Number of molecules processed:'+str(len(dict_output.keys()))}
        print(json.dumps(dict_file))
   

















 
