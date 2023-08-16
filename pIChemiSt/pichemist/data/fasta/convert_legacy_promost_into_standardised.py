import os

from legacy import promost
from legacy_conversion import convert_legacy_set_to_json


def build_promost_set_from_module(pm: promost):
    """Builds a pKa set from a module."""
    pka_set = {
        "acidic": pm.acidic,
        "basic": pm.basic,
        "terminus_ionizable": pm.terminus_ionizable
    }
    return pka_set


if __name__ == '__main__':
    # Get parent dir where the script runs
    parent_dir = os.path.dirname(os.path.realpath(__file__))

    # Converts the legacy data set with filled missing AAs into JSON
    name = promost.SetName
    pka_set = build_promost_set_from_module(promost)
    filepath = f"{parent_dir}/standardised/" \
        f"pka_set_{name.lower()}_standardised.json"
    pka_dict = convert_legacy_set_to_json(pka_set, name, filepath)
