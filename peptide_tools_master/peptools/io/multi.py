import tempfile

from peptools.io.fasta import _is_input_fasta
from peptools.io.file import FileExtension
from peptools.io.structure import _is_input_sdf
from peptools.io.structure import _is_input_smi


def is_input_multiline(input_data):
    return "\n" in input_data


def multiline_input_to_filepath(input_data, params):
    input_list = input_data.split("\n")
    suffix = recognize_input_suffix(input_list)
    tf = tempfile.NamedTemporaryFile(prefix="tmp_peptools", suffix=suffix, delete=False)
    input_data = tf.name
    with open(input_data, "w") as f:
        for line in input_list:
            f.write(line + "\n")
    params.delete_temp_file = True
    return input_data


def recognize_input_suffix(input_list):
    if _is_input_fasta(input_list[0]):
        suffix = FileExtension.FASTA
    elif _is_input_sdf(input_list):
        suffix = FileExtension.SDF
    elif _is_input_smi(input_list[0]):
        suffix = FileExtension.SMI
    return suffix
