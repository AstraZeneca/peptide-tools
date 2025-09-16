import inspect
import sys
from enum import Enum


class BaseEnum(str, Enum):
    """Base class for models."""

    pass


class PKaType(BaseEnum):
    ACIDIC = "acid"
    BASIC = "base"


class PKaMethod(BaseEnum):
    PKA_MATCHER = "pkamatcher"
    ACD = "acd"


class ACDPKaFlag(BaseEnum):
    CLASSIC = "-MPKAAPP"
    GALAS = "-MPKAAPPGALAS"


class InputAttribute(BaseEnum):
    MOL_NAME = "mol_name"
    MOL_OBJECT = "mol_obj"
    MOL_FASTA = "fasta"


class InputFormat(BaseEnum):
    SMILES_STDIN = "smiles_stdin"
    SMILES_FILE = "smiles_file"
    FASTA_STDIN = "fasta_stdin"
    FASTA_FILE = "fasta_file"
    SD_FILE = "sdf"


class InputFileExtension(BaseEnum):
    SDF = "sdf"
    SMILES = "smi"
    FASTA = "fasta"


class OutputAttribute(BaseEnum):
    MOL_NAME = "mol_name"
    SMILES = "SMILES"
    PI = "pI"
    Q_PH7 = "QpH7"
    PI_INTERVAL = "pI_interval"
    PI_INTERVAL_THRESHOLD = "pI_interval_threshold"
    PLOT_FILENAME = "plot_filename"
    PKA_SET = "pKa_set"
    FRAG_BASE_PKA_FASTA = "frag_base_pkas_fasta"
    FRAG_ACID_PKA_FASTA = "frag_acid_pkas_fasta"
    FRAG_BASE_PKA_CALC = "frag_base_pkas_calc"
    FRAG_ACID_PKA_CALC = "frag_acid_pkas_calc"
    FRAG_CONSTANT_QS = "frag_Qs_calc"


class OutputFragAttribute(BaseEnum):
    TYPE = "type"
    COUNT = "count"
    PKA = "pka"
    FRAGMENT = "frag"
    IMAGE = "base64_image"


class OutputFormat(BaseEnum):
    JSON = "json"
    SD_FILE = "sdf"
    CSV_FILE = "csv"
    CONSOLE = "console"


class OutputConsoleFormat(BaseEnum):
    JSON = "json"
    CONSOLE = "console"


class OutputFileFormat(BaseEnum):
    SD_FILE = "sdf"
    CSV_FILE = "csv"


class Models(object):
    """Contains the definitions of all data models."""

    def __init__(self):
        self.definitions = self._get_definitions()

    @staticmethod
    def _get_subclasses_of_enum():
        """Produces a dict of names and subclasses of BaseEnum."""
        classes = dict()
        for name, cls in inspect.getmembers(sys.modules[__name__]):
            if inspect.isclass(cls) and cls is not BaseEnum:
                if issubclass(cls, BaseEnum):
                    classes[name] = cls
        return classes

    @staticmethod
    def _get_values_from_enum(BaseEnum):
        """Generates a list from all values of a class."""
        return [BaseEnum[x].value for x in BaseEnum._member_names_]

    def _get_definitions(self):
        """Builds a dictionary with classes and their values."""
        classes = self._get_subclasses_of_enum()
        return {c: self._get_values_from_enum(c) for c in classes.values()}


MODELS = Models().definitions
