from copy import deepcopy

from peptools.wrapper.adapters import modify_fasta_according_to_cys_keys_mec
from peptools.wrapper.adapters import modify_fasta_according_to_cys_keys_pichemist
from peptools.wrapper.descriptors import calculate_descriptors
from peptools.wrapper.ec import calculate_extinction_coefficient
from peptools.wrapper.liabilities import calculate_liabilities
from peptools.wrapper.pi import calculate_pichemist


def _modify_fasta_if_needed(mol_supply_json, modifier_fn, chem_params):
    """
    Return a deep-copied mol_supply_json with FASTA modified
    only for entries without a mol_obj.
    """
    data = deepcopy(mol_supply_json)

    for _, entry in data.items():
        if entry.get("mol_obj") is None:
            entry["fasta"] = modifier_fn(
                entry["fasta"],
                chem_params.no_free_cys_thiols,
                chem_params.n_disulfide_bonds,
            )

    return data


def run_peptide_master(mol_supply_json, params):
    descriptors_dict = calculate_descriptors(mol_supply_json)

    # pIChemiSt branch
    pichemist_input = _modify_fasta_if_needed(
        mol_supply_json,
        modify_fasta_according_to_cys_keys_pichemist,
        params.chem,
    )
    pichemist_dict = calculate_pichemist(pichemist_input, params)

    # MEC + liabilities branch
    mec_input = _modify_fasta_if_needed(
        mol_supply_json,
        modify_fasta_according_to_cys_keys_mec,
        params.chem,
    )
    ext_coeff_dict = calculate_extinction_coefficient(mec_input, params.run)
    liabilities_dict = calculate_liabilities(mec_input)

    return {
        "output_descriptors": descriptors_dict,
        "output_extn_coeff": ext_coeff_dict,
        "output_pIChemiSt": pichemist_dict,
        "output_liabilities": liabilities_dict,
    }
