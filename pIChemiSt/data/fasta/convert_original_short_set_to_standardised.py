import os
import json

from legacy_conversion import _invert_dict
from pichemist.config import KNOWN_BASIC_RESIDUES
from pichemist.config import KNOWN_ACIDIC_RESIDUES
from pichemist.config import KNOWN_RESIDUES
from pichemist.config import PKA_JSON_TYPE_MATCHING
from pichemist.config import PKA_JSON_INDICES


def enrich_short_pka_set(pka_set):
    """
    Iterates through a short pKa set
    and enriches the data to produce
    a set that is enriched (i.e., ProMoST format).

    """
    basic = dict()
    acidic = dict()
    terminus_ionizable = dict()

    for res in KNOWN_BASIC_RESIDUES:
        if res in pka_set.keys():
            pka = pka_set[res]
            basic[res] = [pka, pka, pka]

    for res in KNOWN_ACIDIC_RESIDUES:
        if res in pka_set.keys():
            pka = pka_set[res]
            acidic[res] = [pka, pka, pka]

    for res in KNOWN_RESIDUES:
        pka_cterm = pka_set['Cterm']
        pka_nterm = pka_set['Nterm']
        terminus_ionizable[res] = [pka_nterm, pka_cterm]

    pka_set = {
     'acidic': acidic,
     'basic': basic,
     'terminus_ionizable': terminus_ionizable
    }
    return pka_set


def convert_original_set_to_standardised_json(pka_set, name,
                                              output_filepath):
    """
    Converts an original short pKa set
    into a standardised JSON.

    """
    pka_dict = dict()
    pka_dict["name"] = name
    values = {
         "acidic": dict(),
         "basic": dict(),
         "terminus_ionizable": dict()
    }
    pka_dict["values"] = values

    # NOTE: new_name and original_name are the same now
    for new_name, original_name in PKA_JSON_TYPE_MATCHING.items():
        pka_type_set = pka_set[original_name]

        # For each pKa type set (e.g. acidic)
        # Create a new set of values with value and position attributes
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
    with open(output_filepath, "w") as f:
        json.dump(pka_dict, f, indent=2)


if __name__ == '__main__':
    # Get parent dir where the script runs
    parent_dir = os.path.dirname(os.path.realpath(__file__))

    # Get all the pKa file paths
    original_pka_sets_dir = f"{parent_dir}/original_short"
    with os.scandir(original_pka_sets_dir) as pka_files:
        pka_filepaths = [f"{original_pka_sets_dir}/{pka_file.name}"
                         for pka_file in pka_files]

    # Convert all original pKa files into standardised
    standardised_pka_sets_dir = f"{parent_dir}/standardised"
    for filepath in pka_filepaths:
        with open(filepath) as f:
            data = json.load(f)

        name = data["name"]
        pka_values = data["values"]
        enriched_pka_set = enrich_short_pka_set(pka_values)
        filepath = f"{standardised_pka_sets_dir}/" \
                   f"pka_set_{name.lower()}_standardised.json"
        convert_original_set_to_standardised_json(enriched_pka_set, name,
                                                  filepath)
