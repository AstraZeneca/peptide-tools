from extn_coeff_fasta import calc_extn_coeff


def calculate_extinction_coefficient(mol_supply_json, params):
    dict_out_extn_coeff = dict()
    if params.calc_extn_coeff:
        extn_coeff_options = _configure_options(mol_supply_json)
        dict_out_extn_coeff = calc_extn_coeff(extn_coeff_options)
    return dict_out_extn_coeff


def _configure_options(mol_supply_json):
    extn_coeff_options = {
        "seq": "",
        "inputDict": mol_supply_json,
        "inputJSON": "",
        "inputFile": "",
        "outputFile": "",
        "l_json": True,
    }
    return extn_coeff_options
