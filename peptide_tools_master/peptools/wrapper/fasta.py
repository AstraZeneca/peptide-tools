from pI_fasta import calc_pI_fasta


def calculate_pI_from_fasta(mol_supply_json, params, chem_params):
    dict_out_pI_fasta = dict()
    if params.calc_pI_fasta:
        ionized_Cterm_residue = ""
        if chem_params.ionized_Cterm:
            ionized_Cterm_residue = "_"

        ionized_Nterm_residue = ""
        if chem_params.ionized_Nterm:
            ionized_Nterm_residue = "_"

        pI_fasta_options = {
            "seq": "",
            "inputDict": mol_supply_json,
            "inputJSON": "",
            "inputFile": "",
            "outputFile": "",
            "tol": 0.001,
            "CTermRes": "_",
            "NTermRes": "_",
            "IonizableTerminiOfCTermRes": ionized_Cterm_residue,
            "IonizableTerminiOfNTermRes": ionized_Nterm_residue,
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
        dict_out_pI_fasta = calc_pI_fasta(pI_fasta_options)
    return dict_out_pI_fasta
