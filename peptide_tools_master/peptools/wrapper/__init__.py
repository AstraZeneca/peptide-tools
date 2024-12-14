from peptools.wrapper.ec import calculate_extinction_coefficient
from peptools.wrapper.pi import calculate_pichemist
from peptools.wrapper.pi import calculate_pifasta


def run_peptide_master(mol_supply_json, params):
    dict_out_extn_coeff = calculate_extinction_coefficient(mol_supply_json, params.run)
    dict_out_pI_fasta = calculate_pifasta(mol_supply_json, params.run, params.chem)
    dict_out_pIChemiSt = calculate_pichemist(mol_supply_json, params.run, params.io)
    dict_out_peptide_tools_master = {
        "output_extn_coeff": dict_out_extn_coeff,
        "output_pI_fasta": dict_out_pI_fasta,
        "output_pIChemiSt": dict_out_pIChemiSt,
    }
    return dict_out_peptide_tools_master
