from pichemist.core import merge_pkas_lists
from pichemist.charges import SmartsChargeCalculator
from pichemist.fasta.matcher import FastaPKaMatcher
from pichemist.pka.acd import calc_pkas_acdlabs
from pichemist.pka.pkamatcher import PKaMatcher
from pichemist.model import PKaMethod
from pichemist.model import MODELS


class ApiExcepton(Exception):
    pass


def get_fasta_pkas_from_list(smiles_list):
    """Match pKas from FASTA definitions for a SMILES list."""
    return FastaPKaMatcher().get_aa_pkas_from_list(smiles_list)


def calculate_pkas_from_list(smiles_list, method):
    """Calculates pKa values using ACD or pKa Matcher."""
    if method not in MODELS[PKaMethod]:
        raise ApiExcepton("Invalid method. Only the formats "
                          f"'{MODELS[PKaMethod]} are accepted")
    if method == PKaMethod.ACD.value:
        base_pkas, acid_pkas, diacid_pkas = calc_pkas_acdlabs(smiles_list)
    if method == PKaMethod.PKA_MATCHER.value:
        base_pkas, acid_pkas, diacid_pkas = \
            PKaMatcher().calculate_pka_from_list(smiles_list)
    return base_pkas, acid_pkas, diacid_pkas


def match_and_calculate_pkas_and_charges_from_list(smiles_list, method):
    """
    Produces the pKa values for a list of SMILES by matching them against
    some FASTA definitions and calculating the unmatched ones.

    """
    # FASTA match
    unknown_frags, base_pkas_fasta, acid_pkas_fasta, diacid_pkas_fasta = \
        get_fasta_pkas_from_list(smiles_list)

    # Unknown fragment calculation
    if len(unknown_frags) > 0:
        base_pkas_calc, acid_pkas_calc, diacid_pkas_calc = \
            calculate_pkas_from_list(unknown_frags,
                                     method=method)

    # Merge FASTA-matched and calculated pKas
    base_pkas_list = [base_pkas_fasta, base_pkas_calc]
    acid_pkas_list = [acid_pkas_fasta, acid_pkas_calc]
    diacid_pkas_list = [diacid_pkas_fasta, diacid_pkas_calc]
    base_pkas_dict, acid_pkas_dict, diacid_pkas_dict = \
        merge_pkas_lists(base_pkas_list, acid_pkas_list, diacid_pkas_list)

    # Calculate the charges and return all results
    net_qs = SmartsChargeCalculator().calculate_net_qs_from_list(unknown_frags)
    return base_pkas_dict, acid_pkas_dict, diacid_pkas_dict, net_qs
