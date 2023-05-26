from pichemist.charges import ChargeCalculator
from pichemist.pka.acd import calc_pkas_acdlabs
from pichemist.pka.pkamatcher import PKaMatcher
from pichemist.model import PKaMethod
from pichemist.model import MODELS


class ApiException(Exception):
    pass


def calc_pkas(smiles_list, method):
    """TODO: """
    if method not in MODELS[PKaMethod]:
        raise ApiException("Invalid method. Only the formats "
                           f"'{MODELS[PKaMethod]} are accepted")
    if method == PKaMethod.ACD.value:
        base_pkas, acid_pkas, diacid_pkas = calc_pkas_acdlabs(smiles_list)
    if method == PKaMethod.PKA_MATCHER.value:
        base_pkas, acid_pkas, diacid_pkas = \
            PKaMatcher().calculate_pka_from_list(smiles_list)
    # TODO: Remove?
    # if use_dimorphite:
    #    base_pkas,acid_pkas,diacid_pkas = calc_pkas_dimorphite_dl(smiles_list)
    net_Qs = ChargeCalculator().calculate_net_qs_from_list(smiles_list)
    return (base_pkas, acid_pkas, diacid_pkas, net_Qs)
