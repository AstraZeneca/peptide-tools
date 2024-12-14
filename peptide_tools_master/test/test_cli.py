import json
import os
import subprocess
import textwrap

import pytest
from helpers import raise_if_file_exists_list
from helpers import raise_if_file_not_exists_list
from helpers import remove_file_list
from helpers import stringify_list
from peptools.io import IOException


test_dir = os.path.dirname(os.path.abspath(__file__))
examples_dir = os.path.join(test_dir, "examples")
peptide_tools_dir = os.path.dirname(test_dir)


cli_base_args = ["python", f"{peptide_tools_dir}/peptide_tools_master.py"]


def test_smiles_stdin_input_1():
    """Validity of console JSON output."""
    smiles = "C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N"  # noqa: E501
    test_args = cli_base_args + ["--input", smiles]
    subprocess_output = subprocess.run(
        stringify_list(test_args), capture_output=True, text=True
    )
    # print(" ".join(test_args))
    result = json.loads(subprocess_output.stdout)
    with open(f"{examples_dir}/payload_1_out.json", "r") as file:
        expected = json.load(file)
    assert result == expected


def test_smiles_stdin_input_1_fragments():
    """Validity of console JSON output."""
    smiles = "C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N"  # noqa: E501
    test_args = cli_base_args + ["--input", smiles, "--print_fragment_pkas"]
    subprocess_output = subprocess.run(
        stringify_list(test_args), capture_output=True, text=True
    )
    # print(" ".join(test_args))
    result = json.loads(subprocess_output.stdout)
    with open(f"{examples_dir}/payload_1_out_fragments.json", "r") as file:
        expected = json.load(file)
    assert result == expected


def test_smiles_file_input_1():
    """Validity of CSV file output."""
    # Prevalidation
    smiles_file = os.path.join(examples_dir, "payload_1.smi")
    temporary_result_file = os.path.join(examples_dir, "payload_1_OUTPUT.csv")
    temporary_plot_file = os.path.join(examples_dir, "payload_1_1.png")
    temporary_file_list = [temporary_result_file, temporary_plot_file]
    try:
        raise_if_file_exists_list(temporary_file_list)

        # Validation
        test_args = cli_base_args + ["--input", smiles_file]
        _ = subprocess.run(stringify_list(test_args), capture_output=True, text=True)
        # print(" ".join(test_args))
        raise_if_file_not_exists_list(temporary_file_list)
        with open(temporary_result_file, "r") as file:
            temp_content = file.read()
        with open(f"{examples_dir}/payload_1_out.csv", "r") as file:
            expected_content = file.read()
        assert (
            temp_content == expected_content
        ), "Expected output file content does not match"
    finally:
        remove_file_list(temporary_file_list)


def test_fasta_file_input_1():
    """Validity of CSV file output."""
    fasta_file = os.path.join(examples_dir, "payload_2.fasta")
    temporary_result_file = os.path.join(examples_dir, "payload_2_OUTPUT.csv")
    try:
        raise_if_file_exists_list([temporary_result_file])
        test_args = cli_base_args + ["--input", fasta_file]
        _ = subprocess.run(stringify_list(test_args), capture_output=True, text=True)
        # print(" ".join(test_args))

        assert os.path.exists(temporary_result_file)
        with open(temporary_result_file, "r") as file:
            temp_content = file.read()
        with open(f"{examples_dir}/payload_2_out.csv", "r") as file:
            expected_content = file.read()
        assert (
            temp_content == expected_content
        ), "Expected output file content does not match"
    finally:
        os.remove(temporary_result_file)


def test_smiles_stdin_input_3():
    """Validity of console text output."""
    smiles = "C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)NCCCCN"  # noqa: E501
    test_args = cli_base_args + [
        "--input",
        smiles,
        "--print_fragment_pkas",
        "--ionized_Cterm",
        "--ionized_Nterm",
        "-p",
        0,
        "-l",
        0,
        "-l",
        0,
    ]
    subprocess_output = subprocess.run(
        stringify_list(test_args), capture_output=True, text=True
    )
    # print(" ".join(stringify_list(test_args)))
    result = json.loads(subprocess_output.stdout)
    with open(f"{examples_dir}/payload_3_out.json", "r") as file:
        expected = json.load(file)
    assert result == expected


def test_smiles_stdin_input_multiline_1():
    """Validity of CSV file output for multiple SMILES."""
    multi_smiles = textwrap.dedent(
        """\
        C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N
        C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N
    """
    )
    test_args = cli_base_args + ["--input", multi_smiles]
    subprocess_output = subprocess.run(
        stringify_list(test_args), capture_output=True, text=True
    )
    # print(" ".join(stringify_list(test_args)))
    result = json.loads(subprocess_output.stdout)
    assert "outputFile" in result
    results_filepath = result["outputFile"]
    raise_if_file_not_exists_list([results_filepath])
    with open(results_filepath, "r") as file:
        result = file.read()

    with open(f"{examples_dir}/payload_4_out.csv", "r") as file:
        expected = file.read()
    assert result == expected
    # TODO: Check deletion of temporary files


def test_fasta_stdin_input_1():
    """Validity of CSV file output for multiple SMILES."""
    fasta = textwrap.dedent(
        """\
        >sp|P43220|GLP1R_HUMAN Glucagon-like peptide 1 receptor OS=Homo sapiens OX=9606 GN=GLP1R PE=1 SV=2
        MAGAPGPLRLALLLLGMVGRAGPRPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDL
        FCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPW
        RDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIH
        LNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAAN
        YYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRN
        SNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLL
        GTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERW
        RLEHLHIQRDSSMKPLKCPTSSLSSGATAGSSMYTATCQASCS
    """
    )
    test_args = cli_base_args + ["--input", fasta]
    subprocess_output = subprocess.run(
        stringify_list(test_args), capture_output=True, text=True
    )
    # print(" ".join(stringify_list(test_args)))
    result = json.loads(subprocess_output.stdout)
    assert "outputFile" in result
    results_filepath = result["outputFile"]
    raise_if_file_not_exists_list([results_filepath])
    with open(results_filepath, "r") as file:
        result = file.read()

    with open(f"{examples_dir}/payload_5_out.csv", "r") as file:
        expected = file.read()
    assert result == expected


def test_sdf_file_input_1():
    """Validity of CSV file output."""
    sdf_file = os.path.join(examples_dir, "payload_5.sdf")
    temporary_result_file = os.path.join(examples_dir, "payload_5_OUTPUT.sdf")
    temporary_plot_file = os.path.join(examples_dir, "payload_5_1.png")
    temporary_file_list = [temporary_result_file, temporary_plot_file]
    try:
        raise_if_file_exists_list(temporary_file_list)

        test_args = cli_base_args + ["--input", sdf_file]
        _ = subprocess.run(stringify_list(test_args), capture_output=True, text=True)
        # print(" ".join(test_args))

        raise_if_file_not_exists_list(temporary_file_list)
        with open(temporary_result_file, "r") as file:
            temp_content = file.read()
        with open(f"{examples_dir}/payload_5_out.sdf", "r") as file:
            expected_content = file.read()
        assert (
            temp_content == expected_content
        ), "Expected output file content does not match"
    finally:
        remove_file_list(temporary_file_list)


def test_empty_input():
    test_args = cli_base_args + ["--input", " "]
    subprocess_output = subprocess.run(
        stringify_list(test_args), capture_output=True, text=True
    )
    assert subprocess_output.returncode == 1
