from smi2scrambledfasta import get_scrambled_fasta_from_mol


class InputFileExtension:
    SDF = ".sdf"
    SMI = ".smi"
    FASTA = ".fasta"


ACCEPTED_FILE_FORMATS = [
    InputFileExtension.SDF,
    InputFileExtension.SMI,
    InputFileExtension.FASTA,
]


class InputAttribute:
    MOL_NAME = "mol_name"
    MOL_OBJECT = "mol_obj"
    FASTA = "fasta"


class InputFactory:
    INPUT_TEMPLATE = {
        InputAttribute.MOL_NAME: None,
        InputAttribute.MOL_OBJECT: None,
        InputAttribute.FASTA: None,
    }

    def new(mol, mol_name, fasta=None):
        input_data = InputFactory.INPUT_TEMPLATE.copy()
        input_data[InputAttribute.MOL_NAME] = mol_name
        input_data[InputAttribute.MOL_OBJECT] = mol
        input_data = InputFactory._set_fasta(input_data, mol, fasta)
        return input_data

    def _set_fasta(input_data, mol, fasta):
        if fasta is not None:
            input_data[InputAttribute.FASTA] = fasta
        else:
            input_data[InputAttribute.FASTA] = get_scrambled_fasta_from_mol(mol)
        return input_data
