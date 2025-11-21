from peptools.wrapper.descriptors import calculate_descriptors
from peptools.wrapper.ec import calculate_extinction_coefficient
from peptools.wrapper.liabilities import calculate_liabilities
from peptools.wrapper.pi import calculate_pichemist

import json
from copy import deepcopy



def modify_fasta_according_to_cys_keys_pichemist(fasta_in,no_free_cys_thiols,n_disulfide_bonds_str):

    if no_free_cys_thiols:
        n_cys_to_replace = fasta_in.count("C")
    else:

        if n_disulfide_bonds_str == "max":
            n_cys = fasta_in.count("C")
            if n_cys % 2 == 1: # odd number
                n_cys_to_replace = n_cys-1
            else:
                n_cys_to_replace = n_cys
            n_disulfide_bonds = n_cys / 2
        else:
            try: 
                n_disulfide_bonds = int(float(n_disulfide_bonds_str))
            except:
                raise ValueError(
                    f"Specified {n_disulfide_bonds} `n_disulfide_bonds` "
                    f"Can not convert it to integer number."
                )

            if n_disulfide_bonds > 0:
                n_cys = n_disulfide_bonds * 2
                n_cys_to_replace = n_disulfide_bonds * 2
                # We can't have more bonds than Cysteine residues
                n_cys_max = fasta_in.count("C")
                if n_cys > n_cys_max:
                    raise ValueError(
                            f"Specified {n_disulfide_bonds} `n_disulfide_bonds` "
                            f"but only {n_cys_max} Cysteine residues are present."
                    )
            else:
                n_cys = 0
                n_cys_to_replace = 0

    fasta_out = fasta_in.replace("C", "X", n_cys_to_replace)

    return fasta_out



def modify_fasta_according_to_cys_keys_mec(fasta_in,no_free_cys_thiols,n_disulfide_bonds_str):

    n_cys_tot = fasta_in.count("C")

    # determine how many disulfide bridged cysterins are there?
    if n_disulfide_bonds_str == "max":
            if n_cys_tot % 2 == 1: # odd number
                n_cys_in_disulfides = n_cys_tot - 1
            else:
                n_cys_in_disulfides = n_cys_tot
    else:
            try: 
                n_disulfide_bonds = int(float(n_disulfide_bonds_str))
            except:
                raise ValueError(
                    f"Specified {n_disulfide_bonds} `n_disulfide_bonds` "
                    f"Can not convert it to integer number."
                )

            n_cys_in_disulfides = n_disulfide_bonds * 2

    # We can't have more bonds than Cysteine residues
    if n_cys_in_disulfides > n_cys_tot:
                    raise ValueError(
                            f"Specified {n_disulfide_bonds} `n_disulfide_bonds` "
                            f"but only {n_cys_tot} Cysteine residues are present."
                    )

    # replace all C involed in disulfide bonds by special character
    fasta_out = fasta_in.replace("C", "𝒞", n_cys_in_disulfides)

    if no_free_cys_thiols:
        # replace all remaining C by X
        fasta_out = fasta_out.replace("C", "X")

    return fasta_out





def run_peptide_master(mol_supply_json, params):
#    print(json.dumps(params.run.__dict__, indent=2, default=str))
#    print(json.dumps(params.chem.__dict__, indent=2, default=str))

    descriptors_dict = calculate_descriptors(mol_supply_json)

    mol_supply_json_init = deepcopy(mol_supply_json)

    ### pichemist
    #print("mol_supply_json pichemist")
    #print(mol_supply_json)
    # Modify sequence accorind to keys for pIChemiSt
    for i,d in enumerate(mol_supply_json):
        if mol_supply_json[i+1]['mol_obj'] == None:
            mol_supply_json[i+1]['fasta'] = modify_fasta_according_to_cys_keys_pichemist(mol_supply_json[i+1]['fasta'],params.chem.no_free_cys_thiols,params.chem.n_disulfide_bonds)
    pichemist_dict = calculate_pichemist(mol_supply_json, params)

    ### MEC & liabilites
    # Modify sequence accorind to keys for MEC
    mol_supply_json = deepcopy(mol_supply_json_init)
    #print("mol_supply_json mec")
    #print(mol_supply_json)

    for i,d in enumerate(mol_supply_json):
        if mol_supply_json[i+1]['mol_obj'] == None:
            mol_supply_json[i+1]['fasta'] = modify_fasta_according_to_cys_keys_mec(mol_supply_json[i+1]['fasta'],params.chem.no_free_cys_thiols,params.chem.n_disulfide_bonds)
    ext_coeff_dict = calculate_extinction_coefficient(mol_supply_json, params.run)
    liabilities_dict = calculate_liabilities(mol_supply_json)

    dict_out = {
        "output_descriptors": descriptors_dict,
        "output_extn_coeff": ext_coeff_dict,
        "output_pIChemiSt": pichemist_dict,
        "output_liabilities": liabilities_dict,
    }
    return dict_out
