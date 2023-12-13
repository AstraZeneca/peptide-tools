import os
import json

from legacy import smarts_notes
from legacy.smarts_pKaMatcher import read_smarts_pKaMatcher


def _add_notes_to_smarts_dict(smarts_dict):
    smarts_dict["notes"] = smarts_notes.notes
    return smarts_dict


def convert_legacy_ss_smarts_to_json(output_filepath):
    """
    Conversion of DAT semi-structured data into JSON.
    It also adds notes on how the data set was created.

    """
    data = read_smarts_pKaMatcher()
    result = dict()
    result_smarts = dict()
    for e in data:
        res_data = list()
        for i in e:
            smarts = i["smarts"]
            name = i["name"]
            res_data.append({
                "pka": i["pka"],
                "pka_std": i["pka_std"],
                "idx": i["ind"],
                "type": i["type"]
            })

            result_smarts[smarts] = {
                "name": name,
                "data": res_data
            }
    result["smarts"] = result_smarts
    result = _add_notes_to_smarts_dict(result)
    with open(output_filepath, "w") as f:
        json.dump(result, f, indent=2)


if __name__ == '__main__':
    # Get parent dir where the script runs
    parent_dir = os.path.dirname(os.path.realpath(__file__))

    # Converts the legacy data set into JSON
    filepath = f"{parent_dir}/standardised/" \
        f"ss_smarts_set_standardised.json"
    convert_legacy_ss_smarts_to_json(filepath)
