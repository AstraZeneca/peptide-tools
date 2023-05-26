from pichemist.charges import calc_net_Qs
from pichemist.pka.acd import calc_pkas_acdlabs
from pichemist.pka.pkamatcher import PKaMatcher
from pichemist.model import PKaMethod
from pichemist.model import MODELS


class ApiException(Exception):
    pass


# calcualtes pKa for the list of smiles. 
def calc_pkas(smi_list, method):

    if method not in MODELS[PKaMethod]:
        raise Exception("TODO...")

    if method == PKaMethod.ACD.value:
        base_pkas,acid_pkas,diacid_pkas = calc_pkas_acdlabs(smi_list)

    if method == PKaMethod.PKA_MATCHER.value:
        base_pkas,acid_pkas,diacid_pkas = PKaMatcher().calculate_pka_from_list(smi_list)

    #if use_dimorphite:
    #    base_pkas,acid_pkas,diacid_pkas = calc_pkas_dimorphite_dl(smi_list)

    net_Qs = calc_net_Qs(smi_list) 
    
    return (base_pkas,acid_pkas,diacid_pkas,net_Qs)