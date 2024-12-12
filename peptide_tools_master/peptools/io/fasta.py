import os

from Bio import SeqIO


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


def read_fasta_file(input_filepath):
    filename, ext = os.path.splitext(input_filepath)
    if not ext == ".fasta":
        raise FileFormatException()

    biosuppl = SeqIO.parse(open(input_filepath), "fasta")
    mol_supply_json = {}
    mol_id = 0
    for biofasta in biosuppl:
        mol_id += 1
        mol_supply_json[mol_id] = {
            "mol_name": biofasta.id,
            "mol_obj": None,
            "fasta": str(biofasta.seq),
        }
    return mol_supply_json
