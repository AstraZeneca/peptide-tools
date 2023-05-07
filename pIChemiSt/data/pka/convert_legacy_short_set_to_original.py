import os
import json

from legacy.short_pka_sets import pKa_sets_short as pka_sets_short


def convert_short_set_to_original_json(pka_sets, name,
                                       output_filepath):
    """
    Converts the legacy data into JSON
    without applying any further processing.

    """
    short_set = pka_sets[name]
    new_short_set = {
        "name": name,
        "values": short_set
    }
    with open(output_filepath, "w") as f:
        json.dump(new_short_set, f, indent=2)


if __name__ == '__main__':
    # Get parent dir where the script runs
    parent_dir = os.path.dirname(os.path.realpath(__file__))

    # Converts each data set from pka_sets_short into JSON
    for name in pka_sets_short.keys():
        filepath = f"{parent_dir}/original_short/" \
            f"pka_set_{name.lower()}_original.json"
        convert_short_set_to_original_json(pka_sets_short, name,
                                           filepath)
