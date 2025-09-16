from peptools.wrapper.descriptors import calculate_descriptors
from peptools.wrapper.ec import calculate_extinction_coefficient
from peptools.wrapper.liabilities import calculate_liabilities
from peptools.wrapper.pi import calculate_pichemist


def run_peptide_master(mol_supply_json, params):
    ext_coeff_dict = calculate_extinction_coefficient(mol_supply_json, params.run)
    pichemist_dict = calculate_pichemist(mol_supply_json, params)
    liabilities_dict = calculate_liabilities(mol_supply_json)
    descriptors_dict = calculate_descriptors(mol_supply_json)
    dict_out = {
        "output_descriptors": descriptors_dict,
        "output_extn_coeff": ext_coeff_dict,
        "output_pIChemiSt": pichemist_dict,
        "output_liabilities": liabilities_dict,
    }
    return dict_out
