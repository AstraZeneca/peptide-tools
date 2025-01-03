import csv
import json
import os

from peptools.io.fasta import _is_input_fasta_sequence
from peptools.io.fasta import configure_fasta_input
from peptools.io.fasta import read_fasta_file
from peptools.io.file import FileFormatException
from peptools.io.input import ACCEPTED_FILE_FORMATS
from peptools.io.input import InputFileExtension
from peptools.io.multi import is_input_multiline
from peptools.io.multi import multiline_input_to_filepath
from peptools.io.output import OutputFileExtension
from peptools.io.params import ChemicalParameters
from peptools.io.params import IOParameters
from peptools.io.params import ParameterSet
from peptools.io.params import RuntimeParameters
from peptools.io.structure import _is_input_smi
from peptools.io.structure import configure_smi_input
from peptools.io.structure import read_structure_file
from rdkit import Chem


class IOException(Exception):
    pass


def generate_input(input_data):
    input_data = input_data.encode("utf-8").decode("unicode_escape")
    params = IOParameters()
    mol_supply_json = dict()
    input_data = _polish_input(input_data, params)

    # Validate input
    if not input_data:
        raise IOException("Input data is empty.")

    # Multiline input
    if is_input_multiline(input_data):
        input_data = multiline_input_to_filepath(input_data, params)

    # Input is a file path
    if os.path.exists(input_data):
        mol_supply_json = read_file(input_data, params)

    # Input is FASTA sequence
    elif _is_input_fasta_sequence(input_data):
        mol_supply_json = configure_fasta_input(input_data, params)

    # Input is SMILES
    elif _is_input_smi(input_data):
        mol_supply_json = configure_smi_input(input_data, params)
    else:
        raise FileFormatException()
    return mol_supply_json, params


def _polish_input(input_data, params):
    input_data = input_data.strip()
    input_data = input_data.replace("ENDOFLINE", "\n")
    return input_data


def read_file(input_data, params):
    params.input_filepath = input_data
    params.workdir = os.path.dirname(params.input_filepath)
    params.filename = os.path.splitext(os.path.basename(params.input_filepath))[0]

    # Validation
    params.filepath_prefix, params.input_file_extension = os.path.splitext(
        params.input_filepath
    )
    if params.input_file_extension not in ACCEPTED_FILE_FORMATS:
        raise FileFormatException(
            "Extension not supported: " + params.input_file_extension
        )

    # Configure output
    _configure_output_file(params)

    # Read file
    if params.input_file_extension in [InputFileExtension.SDF, InputFileExtension.SMI]:
        mol_supply_json = read_structure_file(input_data)
    elif params.input_file_extension == InputFileExtension.FASTA:
        mol_supply_json = read_fasta_file(input_data)
    else:
        raise FileFormatException()

    # Delete if temporary file
    if params.delete_temp_file:
        os.remove(params.input_filepath)
    return mol_supply_json


def _configure_output_file(params):
    OUTPUT_SUFFIX = "_OUTPUT"
    params.output_file_extension = OutputFileExtension.CSV
    if params.input_file_extension == InputFileExtension.SDF:
        params.output_file_extension = OutputFileExtension.SDF
    params.output_filename = (
        f"{params.filepath_prefix}{OUTPUT_SUFFIX}{params.output_file_extension}"
    )


def configure_runtime_parameters(args, input_file_extension):
    params = RuntimeParameters()
    params.generate_plot = False
    params.print_fragment_pkas = bool(args.print_fragment_pkas)
    # TODO: Merge logic since pIfasta has been removed
    if input_file_extension in [
        InputFileExtension.SMI,
        InputFileExtension.SDF,
    ]:
        params.calc_extn_coeff = True
        params.calc_pIChemiSt = True
    elif input_file_extension == InputFileExtension.FASTA:
        params.calc_extn_coeff = True
        params.calc_pIChemiSt = True
    return params


def configure_chemical_parameters(args):
    return ChemicalParameters(
        args.ionizable_cterm,
        args.ionizable_nterm,
        args.NPhosphateGroups,
        args.NAlkylLysGroups,
        args.NDiAlkylLysGroups,
    )


def generate_parameter_set(args, io_params):
    input_file_extension = io_params.input_file_extension
    run_params = configure_runtime_parameters(args, input_file_extension)
    chem_params = configure_chemical_parameters(args)
    return ParameterSet(io_params, run_params, chem_params)


def generate_output(mol_supply_json, dict_out, params):
    pichemist_dict = dict_out["output_pIChemiSt"]
    ext_coeff_dict = dict_out["output_extn_coeff"]
    if not params.io.output_filename:
        _print_to_console_and_exit(dict_out)

    if params.io.input_file_extension != ".fasta":
        mol_list = []
        for mi in mol_supply_json.keys():
            mol = mol_supply_json[mi]["mol_obj"]

            if params.run.calc_pIChemiSt:
                mol.SetProp("pI mean", "%.2f" % pichemist_dict[mi]["pI"]["pI mean"])
                mol.SetProp("pI std", "%.2f" % pichemist_dict[mi]["pI"]["std"])
                mol.SetProp(
                    "pI interval",
                    " - ".join(["%.2f" % x for x in pichemist_dict[mi]["pI_interval"]]),
                )
                mol.SetProp(
                    "pI interval threshold",
                    "%.2f" % pichemist_dict[mi]["pI_interval_threshold"],
                )

            if params.run.calc_extn_coeff:
                mol.SetProp("mol_name", ext_coeff_dict[mi]["mol_name"])
                mol.SetProp("Sequence(FASTA)", ext_coeff_dict[mi]["fasta"])
                mol.SetProp("e205(nm)", "%i" % ext_coeff_dict[mi]["e205"])
                mol.SetProp("e214(nm)", "%i" % ext_coeff_dict[mi]["e214"])
                mol.SetProp("e280(nm)", "%i" % ext_coeff_dict[mi]["e280"])

            mol_list.append(mol)

        if params.io.output_file_extension == ".sdf":
            with Chem.SDWriter(params.io.output_filename) as sdf_w:
                for mol in mol_list:
                    sdf_w.write(mol)

        elif params.io.output_file_extension == ".csv":
            with open(params.io.output_filename, "w") as csv_f:
                csv_w = csv.writer(csv_f)
                count = 0
                for mol in mol_list:
                    props = mol.GetPropsAsDict()

                    count += 1
                    if count == 1:
                        header = ["SMILES"] + list(props.keys())
                        csv_w.writerow(header)

                    row = [Chem.MolToSmiles(mol)]
                    for p in header[1:]:
                        row += [props[p]]
                    csv_w.writerow(row)

    if params.io.input_file_extension == ".fasta":
        if not params.io.output_file_extension == ".csv":
            raise FileFormatException("Only CSV output is supported for FASTA input.")
        dict_list = _generate_fasta_output(mol_supply_json, dict_out, params)
        dict_list_to_file(dict_list, params.io.output_filename)

    dict_file = {
        "outputFile": params.io.output_filename,
        "outputInfo": "Number of molecules processed:"
        + str(len(mol_supply_json.keys())),
    }
    print(json.dumps(dict_file))


def _print_to_console_and_exit(dict_out):
    print(json.dumps(dict_out, indent=2))
    exit()


def _generate_fasta_output(mol_supply_json, dict_out, params):
    ext_coeff_dict = dict_out["output_extn_coeff"]
    pichemist_dict = dict_out["output_pIChemiSt"]

    dict_list = list()
    for mi in mol_supply_json.keys():
        res = dict()
        fasta = mol_supply_json[mi]["fasta"]

        if params.run.calc_pIChemiSt:
            res["pI mean"] = "%.2f" % pichemist_dict[mi]["pI"]["pI mean"]
            res["pI std"] = "%.2f" % pichemist_dict[mi]["pI"]["std"]

        if params.run.calc_extn_coeff:
            res["mol_name"] = ext_coeff_dict[mi]["mol_name"]
            res["Sequence(FASTA)"] = ext_coeff_dict[mi]["fasta"]
            res["e205(nm)"] = "%i" % ext_coeff_dict[mi]["e205"]
            res["e214(nm)"] = "%i" % ext_coeff_dict[mi]["e214"]
            res["e280(nm)"] = "%i" % ext_coeff_dict[mi]["e280"]
        dict_list.append(res)
    return dict_list


def dict_list_to_file(dict_list, filename):
    with open(filename, "w") as file:
        count = 0
        writer = csv.writer(file)
        for attributes in dict_list:
            count += 1
            if count == 1:
                header = list(attributes.keys())
                writer.writerow(header)

            row = []
            for h in header:
                row += [attributes[h]]
            writer.writerow(row)
