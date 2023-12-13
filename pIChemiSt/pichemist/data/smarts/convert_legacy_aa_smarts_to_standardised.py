import os
import json

import legacy.smarts_naturalAAs


def _get_aa_dict_names():
    return [i for i in legacy.smarts_naturalAAs.__dir__()
            if not i.startswith("__")]


def convert_legacy_aa_smarts_to_json(output_filepath):
    """Conversion of AA SMARTS into JSON."""
    result = {"smarts": {}}
    for d in _get_aa_dict_names():
        result["smarts"][d] = getattr(legacy.smarts_naturalAAs, d)
    with open(output_filepath, "w") as f:
        json.dump(result, f, indent=2)


if __name__ == '__main__':
    # Get parent dir where the script runs
    parent_dir = os.path.dirname(os.path.realpath(__file__))

    # Converts the legacy data set into JSON
    filepath = f"{parent_dir}/standardised/" \
        f"aa_smarts_set_standardised.json"
    convert_legacy_aa_smarts_to_json(filepath)
