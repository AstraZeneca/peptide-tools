from pichemist.api import pichemist_from_dict

from pI_fasta import calc_pI_fasta


def calculate_pifasta(mol_supply_json, params, chem_params):
    dict_out_pI_fasta = dict()
    if params.calc_pI_fasta:
        pI_fasta_options = _configure_options(mol_supply_json, params, chem_params)
        dict_out_pI_fasta = calc_pI_fasta(pI_fasta_options)
    return dict_out_pI_fasta


def _configure_options(mol_supply_json, params, chem_params):
    pI_fasta_options = {
        "seq": "",
        "inputDict": mol_supply_json,
        "inputJSON": "",
        "inputFile": "",
        "outputFile": "",
        "tol": 0.001,
        "CTermRes": "_",
        "NTermRes": "_",
        "IonizableTerminiOfCTermRes": _get_ionized_residue(chem_params.ionized_Cterm),
        "IonizableTerminiOfNTermRes": _get_ionized_residue(chem_params.ionized_Nterm),
        "lCyclic": False,
        "NPhosphateGroups": chem_params.NPhosphateGroups,
        "NAlkylLysGroups": chem_params.NAlkylLysGroups,
        "NDiAlkylLysGroups": chem_params.NDiAlkylLysGroups,
        "lPrintpKa": False,
        "lPlot": params.generate_plots,
        "lIgnoreC": False,
        "plot_filename": "OUT_titration_curve.png",
        "l_json": True,
        "pka_set_list": "",
    }
    return pI_fasta_options


def _get_ionized_residue(is_residue_ionized):
    if is_residue_ionized:
        return "_"
    else:
        return ""


def calculate_pichemist(mol_supply_json, params, print_fragment_pkas):
    dict_out_pIChemiSt = dict()
    if params.calc_pIChemiSt:
        plot_filename_prefix = "temp"
        if params.filepath_prefix:
            plot_filename_prefix = params.filepath_prefix

        dict_out_pIChemiSt = pichemist_from_dict(
            mol_supply_json,
            method="pkamatcher",
            ph_q_curve_file_prefix=plot_filename_prefix,
            plot_ph_q_curve=params.generate_plots,
            print_fragments=print_fragment_pkas,
        )
    return dict_out_pIChemiSt
