###
### Last  update: Andrey Frolov, AstraZeneca, Molndal. 08/01/2020
### First verion: Andrey Frolov, AstraZeneca, Molndal. 11/02/2016
###

import sys
import optparse
import json

#import numpy as np
#import pylab

#import json
#from json import encoder
#encoder.FLOAT_REPR = lambda o: format(o, '.2f')


# Conversion maps

pka_json_type_matching = {
     "acidic": "acidic",
     "basic": "basic",
     "terminus_ionizable": "terminus_ionizable"
}

# pka_json_indices = {
#      "primary": 0,
#      "n-terminal": 1,
#      "c-terminal": 2
# }

# def invert_dict(my_map):
#     inv_map = {v: k for k, v in my_map.items()}
#     return inv_map

######
###### Some global definitions
######

all_known_pKa_sets=['ProMoST',
'IPC_peptide',
'IPC2_peptide',
'Gauci',
'Bjellqvist',
'Rodwell',
'Grimsley',
'Thurlkill',
'EMBOSS',
'DTASelect',
'Solomon',
'Sillero',
'Lehninger',
'Toseland',
'Nozaki',
'Dawson']


def list_to_comma_seprated_string(l):
	s=""
	for v in l: s+=str(v)+","
	return s[:-1]

### Preselected set of pKa to display
PKA_SETS_NAMES = ['IPC2_peptide',
                 'IPC_peptide',
                 'ProMoST',
                 'Gauci',
                 'Grimsley',
                 'Thurlkill',
                 'Lehninger',
                 'Toseland']

known_basic_res=['K','R','H']
known_acidic_res=['D','E','C','Y','U']
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
                sys.exit(1)
    
    return terminus_ionizable
    


###------------------------------------------------------------------------------------------
### Sets of pKa 

PKA_SETS = {}


def create_logic_set_from_standardised_json(filepath):
    """
    Reads a pKa set from standardised JSON and converts it
    into a set that is compatible with the legacy logic
    of the software. It returns the set and its name.

    A standardise JSON follows the format:
    {
        "name": $name-of-the-set
        "values": {
            "acidic": {
                $amino-acid: [
                    {"value": $value,
                     "position": $position},
                    {"value": $value,
                     "position": $position}
                ]
            },
            "basic": {...},
            "terminus_ionizable": {...}
        }
    }

    """
    with open(filepath) as f:
        data = json.load(f)

    logic_set = dict()
    values = data["values"]
    for new_name, original_name in pka_json_type_matching.items():

        # Initialise a dict and populate it using
        # a list of values instead of objects
        logic_set[original_name] = dict()
        pka_type_set = values[new_name]
        for aa, values_list in pka_type_set.items():
            value_list = list()
            for val_obj in values_list:
                pka = val_obj["value"]
                value_list.append(pka)
            logic_set[original_name][aa] = value_list
    return logic_set, data["name"]


# Get the root dir of the software repo
import os
root_dir = os.path.dirname(
    os.path.dirname(
    os.path.dirname(os.path.realpath(__file__))))

# Get all the pKa file paths
standardised_pka_sets_dir = f"{root_dir}/data/standardised"
with os.scandir(standardised_pka_sets_dir) as pka_files:
    pka_filepaths = [f"{standardised_pka_sets_dir}/{pka_file.name}"
                        for pka_file in pka_files]
    
    for filepath in pka_filepaths:
        data_set, name = create_logic_set_from_standardised_json(filepath)
        PKA_SETS[name] = data_set





###
### Set ProMoST From http://proteomics.mcw.edu/promost_adv.html
###

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

PKA_SETS['ProMoST']={
 'acidic': acidic,
 'basic': basic,
 'terminus_ionizable': terminus_ionizable
}


# TODO: Creation of pka_set_gauci_original.json from legacy
# TODO: Creation of pka_set_gauci_standardised.json from pka_set_gauci_original.json
# TODO: Creation of ProMoST


if __name__ == '__main__':
    # print(PKA_SETS)
    # name = 'Nozaki'
    # short_set = pKa_sets_short[name]
    # print(short_set)
    # new_short_set = {
    #     "name": name,
    #     "values": short_set
    # }
    # print(new_short_set)
    # filepath = f"data/pka_set_{name.lower()}_original.json"
    # with open(filepath, "w") as f:
    #     json.dump(new_short_set, f, indent=2)
    # exit(0)

    SetName = "Gauci"
    # print(PKA_SETS[SetName])
    # exit(0)

    # Write out
    pka_set = PKA_SETS[SetName]
    print(pka_set)

    # pka_dict = dict()
    # pka_dict["name"] = SetName
    # # pka_values
    # values = {
    #      "acidic": {},
    #      "basic": {},
    #      "terminus_ionizable": {}
    # }
    # pka_dict["values"] = values
    # for new_name, original_name in pka_json_type_matching.items():
    #     pka_type_set = pka_set[original_name]
    #     for aa, value_list in pka_type_set.items():
    #         new_values = list()
    #         for i in range(len(value_list)):
    #             pka = value_list[i]
    #             position = invert_dict(pka_json_indices)[i]
    #             new_values.append({
    #                 "position": position,
    #                 "value": pka
    #             })
    #         values[new_name][aa] = new_values
    
    filepath = "data/pka_set_nozaki_refined.json"
    # with open(filepath, "w") as f:
    #     json.dump(pka_dict, f, indent=2)
    # print(pka_dict)

    # Read in
    # old_dict = read_json_set(filepath)
    # print(old_dict)

    # TODO: Ensure that pka_set_gauci_refined = FillMissingAAtoterminus_ionizable(pka_set_gauci_original)

    # old_dict["terminus_ionizable"] = FillMissingAAtoterminus_ionizable(old_dict["terminus_ionizable"])
    # assert pka_set == old_dict