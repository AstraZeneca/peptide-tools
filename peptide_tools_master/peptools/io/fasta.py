import os

from Bio import SeqIO
from peptools.io.file import FileFormatException
from peptools.io.file import InputFileExtension
from peptools.io.model import InputAttribute
from peptools.io.model import InputFactory


def _is_input_fasta(input_data):
    # https://en.wikipedia.org/wiki/FASTA_format
    if input_data[0] == ">" or input_data[0] == ";":
        return True
    return False


def configure_fasta_input(fasta_str, params):
    params.calc_extn_coeff = True
    params.calc_pI_fasta = True
    params.calc_pIChemiSt = False
    return {1: InputFactory.new(None, params.mol_name, fasta_str)}


def read_fasta_file(input_filepath):
    filename, ext = os.path.splitext(input_filepath)
    if not ext == InputFileExtension.FASTA:
        raise FileFormatException()

    biosuppl = SeqIO.parse(open(input_filepath), "fasta")
    mol_supply_json = {}
    mol_id = 0
    for biofasta in biosuppl:
        mol_id += 1
        mol_supply_json[mol_id] = InputFactory.new(None, biofasta.id, str(biofasta.seq))
    return mol_supply_json
