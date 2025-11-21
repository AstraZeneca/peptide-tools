class ParameterSet:
    def __init__(self, io_params, run_params, chem_params):
        self.io = io_params
        self.run = run_params
        self.chem = chem_params


class IOParameters:
    def __init__(self):
        self.mol_name = "none"
        self.filepath = None
        self.filepath_prefix = None
        self.input_filepath = None
        self.input_file_extension = None
        self.output_filename = None
        self.output_file_extension = None
        self.output_dir = None
        self.delete_temp_file = False


class RuntimeParameters:
    def __init__(self):
        self.generate_plots = True
        self.print_fragment_pkas = False
        self.generate_fragment_images = False
        self.calc_extn_coeff = False
        self.calc_pIChemiSt = False


class ChemicalParameters:
    def __init__(self, ionizable_cterm, ionizable_nterm, n_disulfide_bonds, no_free_cys_thiols):
        self.ionizable_cterm = ionizable_cterm
        self.ionizable_nterm = ionizable_nterm
        self.n_disulfide_bonds = n_disulfide_bonds
        self.no_free_cys_thiols = no_free_cys_thiols
