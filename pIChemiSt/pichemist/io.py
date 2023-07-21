import os
import csv
import json

from rdkit import Chem
from pichemist.model import InputFormat
from pichemist.model import FileExtension
from pichemist.model import OutputFileFormat
from pichemist.model import OutputConsoleFormat
from pichemist.model import MODELS


class FileFormatError(Exception):
    pass


class IOException(Exception):
    pass


def generate_input(input_format, input):
    """Produces an input dictionary compatible with the logic."""
    input_dict = dict()
    if input_format == InputFormat.SMILES_FILE:
        input_dict = read_structure_file(input)
    if input_format == InputFormat.SMILES:
        input_dict[1] = {"mol_name": input,
                         "mol_obj": Chem.MolFromSmiles(input),
                         "fasta": None}
    if input_format == InputFormat.JSON:
        input_dict = json.loads(input)
    return input_dict


def read_structure_file(input_filepath):
    """
    Reads a file containing molecule structures.
    It guesses the input type according to the extension
    of the file.

    """
    _, ext = os.path.splitext(input_filepath)
    if not ext:
        raise FileFormatError("Something wrong with the file "
                              f"{input_filepath}")

    # Initialize file reader
    try:
        if ext[1:] == FileExtension.SDF.value:
            suppl = Chem.SDMolSupplier(input_filepath)
        elif ext[1:] == FileExtension.SMILES.value:
            suppl = Chem.SmilesMolSupplier(input_filepath, titleLine=False)
        else:
            raise FileFormatError("Invalid format. Only the formats "
                                  f"{MODELS[FileExtension]} are accepted")
    except OSError:
        raise OSError(f"File error: Does the file exist? {input_filepath}")

    # Populate input and assign properties
    dict_input = dict()
    uuid = 1
    for mol in suppl:
        if not mol.HasProp("_Name"):
            mol.SetProp("_Name", "tmp_" + str(uuid))
        dict_input[uuid] = {"mol_name": mol.GetProp("_Name"),
                            "mol_obj": mol,
                            "fasta": None}
        uuid += 1
    return dict_input


def _format_results_for_console_output(prop_dict, prop):
    """Prints a formatted output for a dictionary of results."""
    lj = 12
    keys = list(prop_dict.keys())
    keys.remove("std")
    keys.insert(0, "std")
    keys.remove("err")
    keys.insert(0, "err")
    keys.remove(prop + " mean")
    keys.insert(0, prop+" mean")
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
        _format_results_for_console_output(
            dict_output[mol_idx]["pI"], "pI")
        _format_results_for_console_output(
            dict_output[mol_idx]["QpH7"], "Q at pH7.4")

        int_tr = dict_output[mol_idx]["pI_interval_threshold"]
        pka_set = dict_output[mol_idx]["pKa_set"]

        print("\npH interval with charge between %4.1f and %4.1f and "
              "prediction tool: %s" % (-int_tr, int_tr, method))
        print("pI interval: %4.1f - %4.1f" % (dict_output[mol_idx][
              "pI_interval"][0], dict_output[mol_idx]["pI_interval"][1]))

        if print_fragments:
            base_pkas_fasta = dict_output[mol_idx]["base_pkas_fasta"]
            acid_pkas_fasta = dict_output[mol_idx]["acid_pkas_fasta"]
            base_pkas_calc = dict_output[mol_idx]["base_pkas_calc"]
            acid_pkas_calc = dict_output[mol_idx]["acid_pkas_calc"]
            constant_Qs_calc = dict_output[mol_idx]["constant_Qs_calc"]

            # Merge fasta and calculated pKas
            base_pkas = base_pkas_fasta[pka_set] + base_pkas_calc
            acid_pkas = acid_pkas_fasta[pka_set] + acid_pkas_calc
            all_base_pkas = list()
            acid_pkas = list()

            # NOTE: Diacids prints disabled
            # diacid_pkas_fasta = dict_output[mol_idx]["diacid_pkas_fasta"]
            # diacid_pkas_calc = dict_output[mol_idx]["diacid_pkas_calc"]
            # diacid_pkas = diacid_pkas_fasta[pka_set] + diacid_pkas_calc
            # diacid_pkas = list()

            # Zip values and structures
            all_base_pkas, all_base_pkas_smi = list(), list()
            acid_pkas, all_acid_pkas_smi = list(), list()
            if len(base_pkas) != 0:
                all_base_pkas, all_base_pkas_smi = zip(*base_pkas)
            if len(acid_pkas) != 0:
                acid_pkas, all_acid_pkas_smi = zip(*acid_pkas)

            # Print the results
            print("\nList of calculated BASE pKa values with their fragments")
            for pkas, smi in zip(all_base_pkas, all_base_pkas_smi):
                s_pkas = ["%4.1f" % (pkas)]
                print("smiles or AA, base pKa : %-15s %s"
                      % (smi, " ".join(s_pkas)))
            print("\nList of calculated ACID pKa values with their fragments")
            for pkas, smi in zip(acid_pkas, all_acid_pkas_smi):
                s_pkas = ["%4.1f" % (pkas)]
                print("smiles or AA, acid pKa : %-15s %s"
                      % (smi, " ".join(s_pkas)))
            print("\nList of constantly ionized fragments")
            for v in constant_Qs_calc:
                pkas = v[0]
                smi = v[1]
                s_pkas = ["%4.1f" % (pkas)]
                print("smiles, charge : %-15s %s"
                      % (smi, " ".join(s_pkas)))


def _check_supported_output_console_format(output_format):
    """Checks the supported output format."""
    if output_format not in MODELS[OutputConsoleFormat]:
        raise IOException("Only the following output formats "
                          f"are supported: {MODELS[OutputConsoleFormat]}")


def _output_to_console(output_dict, output_format,
                       method, print_fragment_pkas):
    """Deals with console outputs."""
    if output_format == OutputConsoleFormat.JSON:
        _output_json_to_console(output_dict)
    elif output_format == OutputConsoleFormat.CONSOLE:
        _output_text_to_console(output_dict, method, print_fragment_pkas)
    else:
        raise RuntimeError("TODO")


def _check_supported_output_file_format(output_format):
    """Checks the supported output format."""
    if output_format not in MODELS[OutputFileFormat]:
        raise IOException("Only the following file output formats "
                          f"are supported: {MODELS[OutputFileFormat]} "
                          f"(got {output_format})")


def _prepare_output_list(input_dict, output_dict):
    """Merges input and output to prepare a list of results."""
    mol_list = list()
    for mi in input_dict.keys():
        mol = input_dict[mi]["mol_obj"]
        mol.SetProp("pI mean", "%.2f"
                    % output_dict[mi]["pI"]["pI mean"])
        mol.SetProp("pI std", "%.2f"
                    % output_dict[mi]["pI"]["std"])
        # NOTE: String interval is disabled
        # print(output_dict[mi]["pI_interval"])
        # mol.SetProp("pI interval", " - ".join(
        #     ["%.2f" % x for x in output_dict[mi]["pI_interval"]]))
        mol.SetProp("pI interval lower bound", "%.2f"
                    % output_dict[mi]["pI_interval"][0])
        mol.SetProp("pI interval upper bound", "%.2f"
                    % output_dict[mi]["pI_interval"][1])
        mol.SetProp("pI interval threshold", "%.2f"
                    % output_dict[mi]["pI_interval_threshold"])
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
                header = ["SMILES"] + list(props.keys())
                csv_w.writerow(header)
            row = [Chem.MolToSmiles(mol)]
            for p in header[1:]:
                row += [props[p]]
            csv_w.writerow(row)


def _output_to_file(input_dict, output_dict,
                    output_file, output_format):
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
    dict_file = {"output_file": output_file,
                 "output_info": "Number of molecules processed: "
                 f"{len(output_dict.keys())}"}
    print(json.dumps(dict_file))


def output_results(input_dict, output_dict, output_file, output_format,
                   method, print_fragment_pkas):
    """Interface to dispatch different outputs."""
    if not output_file:
        _check_supported_output_console_format(output_format)
        _output_to_console(output_dict, output_format,
                           method, print_fragment_pkas)
    else:
        _check_supported_output_file_format(output_format)
        _output_to_file(input_dict, output_dict,
                        output_file, output_format)
        _print_summary_to_console(output_dict, output_file)
