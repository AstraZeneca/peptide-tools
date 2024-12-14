class FileFormatException(Exception):
    pass


class InputFileExtension:
    SDF = ".sdf"
    SMI = ".smi"
    FASTA = ".fasta"


ACCEPTED_FILE_FORMATS = [
    InputFileExtension.SDF,
    InputFileExtension.SMI,
    InputFileExtension.FASTA,
]
