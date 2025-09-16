import csv
import json
import os

from pichemist.model import InputAttribute
from pichemist.model import InputFileExtension
from pichemist.model import InputFormat
from pichemist.model import MODELS
from pichemist.model import OutputAttribute
from pichemist.model import OutputConsoleFormat
from pichemist.model import OutputFileFormat
from rdkit import Chem


class FileFormatError(Exception):
    pass


class IOException(Exception):
    pass


def generate_input(input_format, input_data):
    """Produces an input dictionary compatible with the API."""
    input_dict = dict()
    if input_format == InputFormat.SMILES_FILE or input_format == InputFormat.SD_FILE:
        input_dict = read_structure_file(input_data)
    if input_format == InputFormat.SMILES_STDIN:
        input_dict[1] = {
            InputAttribute.MOL_NAME.value: input_data,
            InputAttribute.MOL_OBJECT.value: Chem.MolFromSmiles(input_data),
            InputAttribute.MOL_FASTA.value: None,
        }
    if input_format == InputFormat.FASTA_STDIN:
        input_dict[1] = {
            InputAttribute.MOL_NAME.value: input_data,
            InputAttribute.MOL_OBJECT.value: None,
            InputAttribute.MOL_FASTA.value: input_data,
        }
    if input_format == InputFormat.FASTA_FILE:
        input_dict = read_fasta_file(input_data)
    return input_dict


def read_structure_file(input_filepath):
    """
    Reads a file containing molecule structures.
    It guesses the input type according to the extension
    of the file.

    """
    # Initialize file reader
    ext = _get_input_extension(input_filepath)
    try:
        if ext[1:] == InputFileExtension.SDF.value:
            suppl = Chem.SDMolSupplier(input_filepath)
        elif ext[1:] == InputFileExtension.SMILES.value:
            suppl = Chem.SmilesMolSupplier(input_filepath, titleLine=False)
        else:
            raise FileFormatError(
                "Invalid format. Only the formats "
                f"{MODELS[InputFileExtension]} are accepted"
            )
    except OSError:
        raise OSError(f"File error: Does the file exist? {input_filepath}")

    # Populate input and assign properties
    dict_input = dict()
    uuid = 1
    for mol in suppl:
        if not mol.HasProp("_Name"):
            mol.SetProp("_Name", "tmp_" + str(uuid))
        dict_input[uuid] = {
            InputAttribute.MOL_NAME.value: mol.GetProp("_Name"),
            InputAttribute.MOL_OBJECT.value: mol,
            InputAttribute.MOL_FASTA.value: None,
        }
        uuid += 1
    return dict_input


def read_fasta_file(input_filepath):
    """
    Reads a file containing FASTA entries using BioPython.

    """
    from Bio import SeqIO  # conditional import

    # Initialize file reader
    ext = _get_input_extension(input_filepath)
    try:
        if ext[1:] == InputFileExtension.FASTA.value:
            suppl = SeqIO.parse(open(input_filepath), "fasta")
        else:
            raise FileFormatError(
                "Invalid format. Only the formats "
                f"{MODELS[InputFileExtension]} are accepted"
            )
    except OSError:
        raise OSError(f"File error: Does the file exist? {input_filepath}")

    # Populate input and assign properties
    dict_input = dict()
    uuid = 1
    for fasta in suppl:
        dict_input[uuid] = {
            InputAttribute.MOL_NAME.value: fasta.id,
            InputAttribute.MOL_OBJECT.value: None,
            InputAttribute.MOL_FASTA.value: str(fasta.seq),
        }
        uuid += 1
    return dict_input


def _get_input_extension(input_filepath):
    _, ext = os.path.splitext(input_filepath)
    if not ext:
        raise FileFormatError("Something is wrong with the file " f"{input_filepath}")
    return ext


def _format_results_for_console_output(prop_dict, prop):
    """Prints a formatted output for a dictionary of results."""
    lj = 12
    keys = list(prop_dict.keys())
    keys.remove("std")
    keys.insert(0, "std")
    keys.remove("err")
    keys.insert(0, "err")
    keys.remove(prop + " mean")
    keys.insert(0, prop + " mean")
    print("")
    print("=" * 150)
    print(prop)
    print("-" * 33)
    for k in keys:
        p = prop_dict[k]
        print(k.rjust(lj) + "  " + str(round(p, 2)).ljust(lj))
    print("")
    return


def _output_json_to_console(dict_output):
    """Dumps a JSON string to console."""
    print(json.dumps(dict_output, indent=2))


def _output_text_to_console(dict_output, method, print_fragments=False):
    """Generates a formatted text output to console."""
    for mol_idx in dict_output.keys():
        dict_mol = dict_output[mol_idx]
        _format_results_for_console_output(dict_mol[OutputAttribute.PI.value], "pI")
        _format_results_for_console_output(
            dict_mol[OutputAttribute.Q_PH7.value], "Q at pH7.4"
        )

        int_tr = dict_mol[OutputAttribute.PI_INTERVAL_THRESHOLD.value]

        print(
            "\npH interval with charge between %4.1f and %4.1f and "
            "prediction tool: %s" % (-int_tr, int_tr, method)
        )
        print(
            "pI interval: %4.1f - %4.1f"
            % (
                dict_mol["pI_interval"][0],
                dict_mol[OutputAttribute.PI_INTERVAL.value][1],
            )
        )

        if print_fragments:
            frag_base_pkas_fasta = dict_mol[OutputAttribute.FRAG_BASE_PKA_FASTA.value]
            frag_acid_pkas_fasta = dict_mol[OutputAttribute.FRAG_ACID_PKA_FASTA.value]
            frag_base_pkas_calc = dict_mol[OutputAttribute.FRAG_BASE_PKA_CALC.value]
            frag_acid_pkas_calc = dict_mol[OutputAttribute.FRAG_ACID_PKA_CALC.value]
            frag_constant_Qs_calc = dict_mol[OutputAttribute.FRAG_CONSTANT_QS.value]

            # Print the results
            print("\nList of calculated pKa values or constant charges")
            format_header = "{:15s} {:5s}  {:17s}  {:s}"
            format_results = "{:15s} {:5d}  {:17.1f}  {:s}"
            print(
                format_header.format("Type", "Count", "pKa_or_constant_Q", "Fragment")
            )

            for _, frag in frag_base_pkas_fasta.items():
                print(
                    format_results.format(
                        frag["type"], frag["count"], frag["pka"], frag["frag"]
                    )
                )

            for _, frag in frag_base_pkas_calc.items():
                print(
                    format_results.format(
                        frag["type"], frag["count"], frag["pka"], frag["frag"]
                    )
                )

            for _, frag in frag_acid_pkas_fasta.items():
                print(
                    format_results.format(
                        frag["type"], frag["count"], frag["pka"], frag["frag"]
                    )
                )

            for _, frag in frag_acid_pkas_calc.items():
                print(
                    format_results.format(
                        frag["type"], frag["count"], frag["pka"], frag["frag"]
                    )
                )

            for _, frag in frag_constant_Qs_calc.items():
                print(
                    format_results.format(
                        frag["type"], frag["count"], frag["pka"], frag["frag"]
                    )
                )


def output_results(
    input_dict, output_dict, output_file, output_format, method, print_fragment_pkas
):
    """Interface to dispatch different outputs."""
    if not output_file:
        _check_supported_output_console_format(output_format)
        _output_to_console(output_dict, output_format, method, print_fragment_pkas)
    else:
        _check_supported_output_file_format(output_format)
        _output_to_file(input_dict, output_dict, output_file, output_format)
        _print_summary_to_console(output_dict, output_file)


def _check_supported_output_console_format(output_format):
    """Checks the supported output format."""
    if output_format not in MODELS[OutputConsoleFormat]:
        raise IOException(
            "Only the following output formats "
            f"are supported: {MODELS[OutputConsoleFormat]}"
        )


def _output_to_console(output_dict, output_format, method, print_fragment_pkas):
    """Deals with console outputs."""
    if output_format == OutputConsoleFormat.JSON:
        _output_json_to_console(output_dict)
    elif output_format == OutputConsoleFormat.CONSOLE:
        _output_text_to_console(output_dict, method, print_fragment_pkas)


def _check_supported_output_file_format(output_format):
    """Checks the supported output format."""
    if output_format not in MODELS[OutputFileFormat]:
        raise IOException(
            "Only the following file output formats "
            f"are supported: {MODELS[OutputFileFormat]} "
            f"(got {output_format})"
        )


def _prepare_output_list(input_dict, output_dict):
    """Merges input and output to prepare a list of results."""
    mol_list = list()
    for mi in input_dict.keys():
        mol = input_dict[mi][InputAttribute.MOL_OBJECT.value]
        mol.SetProp(
            "pI mean", "%.2f" % output_dict[mi][OutputAttribute.PI.value]["pI mean"]
        )
        mol.SetProp("pI std", "%.2f" % output_dict[mi][OutputAttribute.PI.value]["std"])
        mol.SetProp(
            "pI interval lower bound",
            "%.2f" % output_dict[mi][OutputAttribute.PI_INTERVAL.value][0],
        )
        mol.SetProp(
            "pI interval upper bound",
            "%.2f" % output_dict[mi][OutputAttribute.PI_INTERVAL.value][1],
        )
        mol.SetProp(
            "pI interval threshold",
            "%.2f" % output_dict[mi][OutputAttribute.PI_INTERVAL_THRESHOLD.value],
        )
        mol_list.append(mol)
    return mol_list


def _output_csv_to_file(mol_list, output_file):
    """Dumps the results to a CSV file."""
    with open(output_file, "w") as csv_f:
        csv_w = csv.writer(csv_f)
        count = 0
        for mol in mol_list:
            props = mol.GetPropsAsDict()
            count += 1
            if count == 1:
                header = [OutputAttribute.SMILES.value] + list(props.keys())
                csv_w.writerow(header)
            row = [Chem.MolToSmiles(mol)]
            for p in header[1:]:
                row += [props[p]]
            csv_w.writerow(row)


def _output_to_file(input_dict, output_dict, output_file, output_format):
    """Deals with file outputs."""
    mol_list = _prepare_output_list(input_dict, output_dict)
    if output_format == OutputFileFormat.SD_FILE:
        with Chem.SDWriter(output_file) as sdf_w:
            for mol in mol_list:
                sdf_w.write(mol)

    elif output_format == OutputFileFormat.CSV_FILE:
        _output_csv_to_file(mol_list, output_file)


def _print_summary_to_console(output_dict, output_file):
    """Print summary to console."""
    dict_file = {
        "output_file": output_file,
        "output_info": "Number of molecules processed: " f"{len(output_dict.keys())}",
    }
    print(json.dumps(dict_file))
