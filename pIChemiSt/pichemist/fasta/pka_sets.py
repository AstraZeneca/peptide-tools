import os
import json

from pichemist.config import PKA_DATA_SETS_DIR
from pichemist.config import PKA_JSON_TYPE_MATCHING


def create_logic_set_from_standardised_json(filepath):
    """
    Reads a pKa set from standardised JSON and converts it
    into a set that is compatible with the legacy logic
    of the software. It returns the set and its name.

    A standardise JSON follows the format:
    {
        "name": $name-of-the-set
        "values": {
            "acidic": {
                $amino-acid: [
                    {"value": $value,
                     "position": $position},
                    {"value": $value,
                     "position": $position}
                ]
            },
            "basic": {...},
            "terminus_ionizable": {...}
        }
    }

    """
    with open(filepath) as f:
        data = json.load(f)

    logic_set = dict()
    values = data["values"]
    for new_name, original_name in PKA_JSON_TYPE_MATCHING.items():

        # Initialise a dict and populate it using
        # a list of values instead of objects
        logic_set[original_name] = dict()
        pka_type_set = values[new_name]
        for aa, values_list in pka_type_set.items():
            value_list = list()
            for val_obj in values_list:
                pka = val_obj["value"]
                value_list.append(pka)
            logic_set[original_name][aa] = value_list
    return logic_set, data["name"]


def _generate_pka_sets():
    """
    Produces a set of pKa sets that are compatible
    with the business logic.

    """
    # Get all the pKa file paths
    pka_sets = dict()
    standardised_pka_sets_dir = PKA_DATA_SETS_DIR
    with os.scandir(standardised_pka_sets_dir) as pka_files:
        pka_filepaths = [f"{standardised_pka_sets_dir}/{pka_file.name}"
                         for pka_file in pka_files]

        for filepath in pka_filepaths:
            data_set, name = create_logic_set_from_standardised_json(filepath)
            pka_sets[name] = data_set
    return pka_sets


PKA_SETS = _generate_pka_sets()
