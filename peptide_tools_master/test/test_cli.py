import json
import os
import subprocess


test_dir = os.path.dirname(os.path.abspath(__file__))
examples_dir = os.path.join(test_dir, "examples")
peptide_tools_dir = os.path.dirname(test_dir)


cli_base_args = ["python", f"{peptide_tools_dir}/peptide_tools_master.py"]


def test_smiles_stdin_input_1():
    """Validity of console JSON output."""
    smiles = "C[C@@H](C(=O)N[C@@H](CCC(=O)O)C(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc1ccc(cc1)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc3ccccc3)N"  # noqa: E501
    test_args = cli_base_args + ["--input", smiles]
    subprocess_output = subprocess.run(test_args, capture_output=True, text=True)
    result = json.loads(subprocess_output.stdout)
    with open(f"{examples_dir}/payload_1_out.json", "r") as file:
        expected = json.load(file)
    assert result == expected


def test_smiles_file_input_1():
    """Validity of CSV file output."""
    smiles_file = os.path.join(examples_dir, "payload_1.smi")
    temporary_result_file = os.path.join(examples_dir, "payload_1_OUTPUT.csv")
    assert not os.path.exists(
        temporary_result_file
    ), "Expected output file exist prior creation"

    test_args = cli_base_args + ["--input", smiles_file]
    _ = subprocess.run(test_args, capture_output=True, text=True)

    assert os.path.exists(temporary_result_file)
    with open(temporary_result_file, "r") as file:
        temp_content = file.read()
    with open(f"{examples_dir}/payload_1_out.csv", "r") as file:
        expected_content = file.read()
    assert (
        temp_content == expected_content
    ), "Expected output file content does not match"
    os.remove(temporary_result_file)


def test_fasta_file_input_1():
    """Validity of CSV file output."""
    fasta_file = os.path.join(examples_dir, "payload_2.fasta")
    temporary_result_file = os.path.join(examples_dir, "payload_2_OUTPUT.csv")
    assert not os.path.exists(
        temporary_result_file
    ), "Expected output file exist prior creation"

    test_args = cli_base_args + ["--input", fasta_file]
    subprocess_output = subprocess.run(test_args, capture_output=True, text=True)

    assert os.path.exists(temporary_result_file)
    with open(temporary_result_file, "r") as file:
        temp_content = file.read()
    with open(f"{examples_dir}/payload_2_out.csv", "r") as file:
        expected_content = file.read()
    assert (
        temp_content == expected_content
    ), "Expected output file content does not match"
    os.remove(temporary_result_file)
