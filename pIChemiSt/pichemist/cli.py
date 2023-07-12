import os
import sys
import csv
import json
import argparse

from rdkit import Chem
from rdkit import RDLogger
from pichemist.io import generate_input
from pichemist.io import output_property_dict
from pichemist.model import PKaMethod
from pichemist.model import InputFormat
from pichemist.model import OutputFormat
from pichemist.model import MODELS
from pichemist.utils import get_logger

# Configure logging
log = get_logger(__name__)
RDLogger.DisableLog("rdApp.info")


def print_output(dict_output, method, print_fragments=False):

    for mol_idx in dict_output.keys():
        molid = dict_output[mol_idx]
        output_property_dict(dict_output[mol_idx]['pI'], 'pI')
        output_property_dict(dict_output[mol_idx]['QpH7'],'Q at pH7.4')

        if method == "acd":
            prediction_tool = 'ACDlabs'
        if method == "pkamatcher":
            prediction_tool = 'pKaMatcher'

        int_tr = dict_output[mol_idx]['pI_interval_threshold']
        pka_set = dict_output[mol_idx]['pKa_set']

        print(" ")
        print("pH interval with charge between %4.1f and %4.1f and prediction tool: %s" % (-int_tr, int_tr, prediction_tool) )
        print("pI interval: %4.1f - %4.1f" % (dict_output[mol_idx]['pI_interval'][0], dict_output[mol_idx]['pI_interval'][1]))

        if print_fragments:
            base_pkas_fasta = dict_output[mol_idx]['base_pkas_fasta']
            acid_pkas_fasta = dict_output[mol_idx]['acid_pkas_fasta']
            #diacid_pkas_fasta = dict_output[mol_idx]['diacid_pkas_fasta']
            base_pkas_calc = dict_output[mol_idx]['base_pkas_calc']
            acid_pkas_calc = dict_output[mol_idx]['acid_pkas_calc']
            #diacid_pkas_calc = dict_output[mol_idx]['diacid_pkas_calc']
            constant_Qs_calc = dict_output[mol_idx]['constant_Qs_calc']

            # merge fasta and calcualted pkas
            base_pkas = base_pkas_fasta[pka_set] + base_pkas_calc
            acid_pkas = acid_pkas_fasta[pka_set] + acid_pkas_calc
            #diacid_pkas = diacid_pkas_fasta[pka_set] + diacid_pkas_calc
            all_base_pkas=[]
            acid_pkas=[]
            diacid_pkas=[]
            if len(base_pkas) != 0: all_base_pkas,all_base_pkas_smi = zip(*base_pkas) 
            else: all_base_pkas,all_base_pkas_smi = [],[]
            if len(acid_pkas) != 0: acid_pkas,all_acid_pkas_smi = zip(*acid_pkas) 
            else: acid_pkas,all_acid_pkas_smi = [],[]
            #if len(diacid_pkas) != 0: diacid_pkas,all_diacid_pkas_smi = zip(*diacid_pkas) 
            #else: diacid_pkas,all_diacid_pkas_smi = [],[]
        
            print(" ")
            print("List of calculated BASE pKa's with the corresponding fragments")
            for pkas,smi in zip(all_base_pkas,all_base_pkas_smi):
                s_pkas = ["%4.1f"%(pkas)]
                print("smiles or AA, base pKa : %-15s %s" % (smi,' '.join(s_pkas)))

            print(" ")
            print("List of calculated ACID pKa's with the corresponding fragments")
            for pkas,smi in zip(acid_pkas,all_acid_pkas_smi):
                s_pkas = ["%4.1f"%(pkas)]
                print("smiles or AA, acid pKa : %-15s %s" % (smi,' '.join(s_pkas)))

#            print(" ")
#            print("List of calculated DIACID pKa's with the corresponding fragments")
#            for pkas,smi in zip(diacid_pkas,all_diacid_pkas_smi):
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



__doc__ = """Calculates the isoelectric point (pI) of a molecules by cutting its 
amide bonds, retreiving the pKa values for known AAs, predicting pKas of unknown 
fragments using pKaMatcher or ACDlabs, and finally calculating the pI using the 
Henderson-Hasselbalch equation."""

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
                        default=OutputFormat.CONSOLE)
    parser.add_argument("--plot_titration_curve", default=False,
                        action='store_true', dest="plot_titration_curve",
                        help="Generate an image of the titration curve into a file")
    parser.add_argument("--print_fragment_pkas", default=False,
                        action='store_true', dest="print_fragment_pkas",
                        help="TODO: Print out fragments with corresponding pKas used in pI calcution.")
    parser.add_argument("--method",
                        choices=MODELS[PKaMethod],
                        default=PKaMethod.PKA_MATCHER,
                        help="Method for the prediction of pKas of unknown fragments")
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = arg_parser()
    input_dict = generate_input(args.input_format, args.input)
    dict_output = pichemist_from_list(input_dict, args.method,
                                      args.plot_titration_curve,
                                      args.print_fragment_pkas)

    ### ----------------------------------------------------------------------
    # Output
    if args.output_file is None: # output plain text
        if args.output_format == OutputFormat.JSON:
            print(json.dumps(dict_output, indent=2))
        elif args.output_format == OutputFormat.CONSOLE:
            print_output(dict_output, args.method, args.print_fragment_pkas)
        else:
            raise RuntimeError("TODO")  

    else: # output file

        known_out_file_types = ['.sdf','.csv']
        filename, out_fext = os.path.splitext(args.output_file)
        if out_fext not in known_out_file_types:
            raise Exception('Error! Output file extention not in supported file types:' + str(known_out_file_types))
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
   

















 
