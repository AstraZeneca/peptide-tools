from pichemist.core import merge_pkas_lists
from pichemist.charges import SmartsChargeCalculator
from pichemist.charges import PKaChargeCalculator
from pichemist.fasta.matcher import FastaPKaMatcher
from pichemist.isoelectric import CurveCalculator
from pichemist.isoelectric import IsoelectricPointCalculator
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
    base_pkas_calc, acid_pkas_calc, diacid_pkas_calc = \
        calculate_pkas_from_list(unknown_frags,
                                 method=method)

    # Merge FASTA-matched and calculated pKas
    base_pkas_list = [base_pkas_fasta, base_pkas_calc]
    acid_pkas_list = [acid_pkas_fasta, acid_pkas_calc]
    diacid_pkas_list = [diacid_pkas_fasta, diacid_pkas_calc]
    base_pkas_dict, acid_pkas_dict, diacid_pkas_dict = \
        merge_pkas_lists(base_pkas_list, acid_pkas_list, diacid_pkas_list)

    # Calculate charges and their matching fragments, and return results
    net_qs_and_frags = SmartsChargeCalculator().calculate_net_qs_from_list(
        unknown_frags)
    return base_pkas_fasta, acid_pkas_fasta, diacid_pkas_fasta, \
        base_pkas_dict, acid_pkas_dict, diacid_pkas_dict, net_qs_and_frags


def _get_net_qs_from_qs_and_frags(net_qs_and_frags):
    """TODO:..."""
    return [v[0] for v in net_qs_and_frags]


def calculate_pI_pH_and_charge_dicts(base_pkas_dict, acid_pkas_dict,
                                     diacid_pkas_dict, net_qs_and_frags):
    """TODO:..."""
    pI_dict = dict()
    q_dict = dict()
    pH_q_dict = dict()
    pka_sets_names = FastaPKaMatcher().get_pka_sets_names()
    for pka_set in pka_sets_names:

        # Merge fasta and calculated pkas
        base_pkas = base_pkas_dict[pka_set]
        acid_pkas = acid_pkas_dict[pka_set]
        diacid_pkas = diacid_pkas_dict[pka_set]

        # Calculate isoelectric point
        net_qs = _get_net_qs_from_qs_and_frags(net_qs_and_frags)
        constant_q = PKaChargeCalculator().calculate_constant_charge(net_qs)
        Q = PKaChargeCalculator().calculate_charge(base_pkas, acid_pkas,
                                                   diacid_pkas, pH=7.4,
                                                   constant_q=constant_q)
        pI = IsoelectricPointCalculator().calculate_pI(base_pkas, acid_pkas,
                                                       diacid_pkas,
                                                       constant_q=constant_q)
        pH_Q = CurveCalculator().calculate_charged_curve(base_pkas, acid_pkas,
                                                         diacid_pkas,
                                                         constant_q=constant_q)
        pI_dict[pka_set] = pI
        q_dict[pka_set] = Q
        pH_q_dict[pka_set] = pH_Q
    return pI_dict, q_dict, pH_q_dict
