import sys
import inspect
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


class InputFormat(BaseEnum):
    SMILES = "smiles"
    SMILES_FILE = "smiles_file"
    JSON = "json"
    SD_FILE = "sdf"


class OutputFormat(BaseEnum):
    JSON = "json"
    SD_FILE = "sdf"
    CSV_FILE = "csv"


class Models(object):
    """Contains the definitions of all models."""

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
        return [BaseEnum[x].value
                for x in BaseEnum._member_names_]

    def _get_definitions(self):
        """Builds a dictionary with classes and their values."""
        classes = self._get_subclasses_of_enum()
        return {c: self._get_values_from_enum(c)
                for c in classes.values()}


MODELS = Models().definitions
