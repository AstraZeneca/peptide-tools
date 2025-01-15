from pichemist.api import pichemist_from_dict


def calculate_pichemist(mol_supply_json, params):
    dict_out_pIChemiSt = dict()
    if params.run.calc_pIChemiSt:
        plot_filename_prefix = "temp"
        if params.io.filepath_prefix:
            plot_filename_prefix = params.io.filepath_prefix

        dict_out_pIChemiSt = pichemist_from_dict(
            mol_supply_json,
            method="pkamatcher",
            ph_q_curve_file_prefix=plot_filename_prefix,
            plot_ph_q_curve=params.run.generate_plots,
            print_fragments=params.run.print_fragment_pkas,
            ionizable_nterm=params.chem.ionizable_nterm,
            ionizable_cterm=params.chem.ionizable_cterm,
        )
    return dict_out_pIChemiSt
