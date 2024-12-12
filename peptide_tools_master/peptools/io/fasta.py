def _is_input_fasta(input_list):
    # https://en.wikipedia.org/wiki/FASTA_format
    if input_list[0] == ">" or input_list[0] == ";":
        return True
    return False


def configure_fasta_input(fasta_str, params):
    params.calc_extn_coeff = True
    params.calc_pI_fasta = True
    params.calc_pIChemiSt = False
    return {1: {"mol_name": params.mol_name, "mol_obj": None, "fasta": fasta_str}}
