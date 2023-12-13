import json

from pichemist.config import SS_SMARTS_PKA_SET_FILEPATH


def create_logic_set_from_standardised_json(filepath):
    """
    Reads a pKa set from standardised JSON and converts it
    into an set compatible with the legacy logic
    of the software.

    """
    with open(filepath) as f:
        data = json.load(f)

    logic_set = list()
    for smarts, v in data["smarts"].items():
        res_data = list()
        name = v["name"]
        data = v["data"]
        for d in data:
            res_data.append({
                "pka": d["pka"],
                "idx": d["idx"],
                "pka_std": d["pka_std"],
                "type": d["type"],
                "smarts": smarts,
                "name": name
            })
        logic_set.append(res_data)
    return logic_set


# Read in the data
SS_SMARTS_PKA_SET = create_logic_set_from_standardised_json(
    SS_SMARTS_PKA_SET_FILEPATH)
