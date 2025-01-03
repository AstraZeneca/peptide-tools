from peptools.wrapper.ec import calculate_extinction_coefficient
from peptools.wrapper.pi import calculate_pichemist
from peptools.wrapper.pi import calculate_pifasta


def run_peptide_master(mol_supply_json, params):
    ext_coeff_dict = calculate_extinction_coefficient(mol_supply_json, params.run)
    pifasta_dict = calculate_pifasta(mol_supply_json, params.run, params.chem)
    pichemist_dict = calculate_pichemist(mol_supply_json, params.run, params.io)
    dict_out = {
        "output_extn_coeff": ext_coeff_dict,
        "output_pI_fasta": pifasta_dict,
        "output_pIChemiSt": pichemist_dict,
    }
    return dict_out
