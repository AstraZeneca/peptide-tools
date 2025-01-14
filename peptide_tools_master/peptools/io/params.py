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
        self.calc_extn_coeff = False
        self.calc_pIChemiSt = False


# TODO: Remove sic params
class ChemicalParameters:
    def __init__(
        self,
        ionizable_cterm,
        ionizable_nterm,
        NPhosphateGroups,  # sic
        NAlkylLysGroups,  # sic
        NDiAlkylLysGroups,  # sic
    ):
        self.ionizable_cterm = ionizable_cterm
        self.ionizable_nterm = ionizable_nterm
        self.NPhosphateGroups = NPhosphateGroups  # sic
        self.NAlkylLysGroups = NAlkylLysGroups  # sic
        self.NDiAlkylLysGroups = NDiAlkylLysGroups  # sic
