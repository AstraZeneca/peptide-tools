import json

from pichemist.config import PKA_JSON_TYPE_MATCHING
from pichemist.config import PKA_JSON_INDICES


def _invert_dict(my_map):
    """Inverts keys and values of a dictionary."""
    inv_map = {v: k for k, v in my_map.items()}
    return inv_map


def convert_legacy_set_to_json(pka_set, name,
                               output_filepath):
    """
    Converts the legacy data into JSON.

    """
    pka_dict = dict()
    pka_dict["name"] = name
    values = {
            "acidic": {},
            "basic": {},
            "terminus_ionizable": {}
    }

    for new_name, original_name in PKA_JSON_TYPE_MATCHING.items():
        pka_type_set = pka_set[original_name]
        for aa, value_list in pka_type_set.items():
            new_values = list()
            for i in range(len(value_list)):
                pka = value_list[i]
                position = _invert_dict(PKA_JSON_INDICES)[i]
                new_values.append({
                    "position": position,
                    "value": pka
                })
            values[new_name][aa] = new_values

    pka_dict["values"] = values
    with open(output_filepath, "w") as f:
        json.dump(pka_dict, f, indent=2)
