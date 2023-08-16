import os

from legacy import gauci
from legacy_conversion import convert_legacy_set_to_json


def build_gauci_set_from_module(gm: gauci,
                                fill_missing_aas=False):
    """Builds a pKa set from a module."""
    pka_set = {
        "acidic": gm.acidic1,
        "basic": gm.basic1,
        "terminus_ionizable": gm.terminus_ionizable1
    }
    if fill_missing_aas:
        pka_set = _fill_missing_terminus_ionizable_aas(pka_set)
    return pka_set


def _fill_missing_terminus_ionizable_aas(pka_set):
    """Adds the missing terminus ionizable AAs."""
    pka_set["terminus_ionizable"] = \
        gauci.FillMissingAAtoterminus_ionizable(
            pka_set["terminus_ionizable"]
        )
    return pka_set


if __name__ == '__main__':
    # Get parent dir where the script runs
    parent_dir = os.path.dirname(os.path.realpath(__file__))

    # Converts the original legacy data set into JSON
    name = gauci.SetName
    pka_set = build_gauci_set_from_module(gauci, fill_missing_aas=False)
    filepath = f"{parent_dir}/original_enriched/" \
        f"pka_set_{name.lower()}_original.json"
    pka_dict = convert_legacy_set_to_json(pka_set, name, filepath)

    # Converts the legacy data set with filled missing AAs into JSON
    pka_set = build_gauci_set_from_module(gauci, fill_missing_aas=True)
    filepath = f"{parent_dir}/standardised/" \
        f"pka_set_{name.lower()}_standardised.json"
    pka_dict = convert_legacy_set_to_json(pka_set, name, filepath)
