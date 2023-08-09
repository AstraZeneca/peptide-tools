import os
import base64
import json
import hashlib
import pytest
import tempfile

from helpers import TestError
from helpers import examples_dir
from helpers import stdout_to_variable
from pichemist.cli import arg_parser
from pichemist.cli import run_cli


def test_parser_creation():
    """Parser creation."""
    args = arg_parser(["-i", f"{examples_dir}/payload_1.smi",
                       "--print_fragment_pkas", "--method", "pkamatcher"])
    assert args is not None


def test_console_text_output_1():
    """Validity of console text output using pKaMatcher."""
    args = arg_parser(["-i", f"{examples_dir}/payload_1.smi",
                       "--print_fragment_pkas", "--method", "pkamatcher"])
    result = stdout_to_variable(run_cli, args)
    with open(f"{examples_dir}/payload_1_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_console_text_output_2():
    """Validity of console text output using pKaMatcher."""
    args = arg_parser(["-i", f"{examples_dir}/payload_3.smi",
                       "--print_fragment_pkas", "--method", "pkamatcher"])
    result = stdout_to_variable(run_cli, args)
    with open(f"{examples_dir}/payload_3_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


@pytest.mark.acd
def test_console_text_output_3():
    """Validity of console text output using ACD."""
    args = arg_parser(["-i", f"{examples_dir}/payload_5.smi",
                       "--print_fragment_pkas", "--method", "acd"])
    result = stdout_to_variable(run_cli, args)
    with open(f"{examples_dir}/payload_5_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_console_json_output():
    """Validity of console JSON output."""
    args = arg_parser(["-i", f"{examples_dir}/payload_1.smi",
                       "-of", "json", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    result = stdout_to_variable(run_cli, args)
    result = json.loads(result)
    with open(f"{examples_dir}/payload_1_out.json", "r") as f:
        expected = json.load(f)
    assert result == expected


def test_file_csv_output():
    """Validity of CSV file output."""
    tmp_filepath = tempfile.NamedTemporaryFile(suffix='.csv').name
    args = arg_parser(["-i", f"{examples_dir}/payload_1.smi",
                       "-o", tmp_filepath,
                       "-of", "csv", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    run_cli(args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(f"{examples_dir}/payload_1_out.csv", "r") as f:
        expected = f.read()
        with open(tmp_filepath, "r") as f:
            results = f.read()
            assert results == expected
    os.remove(tmp_filepath)


def test_file_sdf_output_1():
    """Validity of SDF file output."""
    tmp_filepath = tempfile.NamedTemporaryFile(suffix='.sdf').name
    args = arg_parser(["-i", f"{examples_dir}/payload_1.smi",
                       "-o", tmp_filepath,
                       "-of", "sdf", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    run_cli(args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(f"{examples_dir}/payload_1_out.sdf", "r") as f:
        expected = f.read()
        with open(tmp_filepath, "r") as f:
            results = f.read()
            assert results == expected
    os.remove(tmp_filepath)


def test_file_sdf_output_2():
    """Validity of SDF file output."""
    tmp_filepath = tempfile.NamedTemporaryFile(suffix='.sdf').name
    args = arg_parser(["-i", f"{examples_dir}/payload_3.smi",
                       "-o", tmp_filepath,
                       "-of", "sdf", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    run_cli(args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(f"{examples_dir}/payload_3_out.sdf", "r") as f:
        expected = f.read()
        with open(tmp_filepath, "r") as f:
            results = f.read()
            assert results == expected
    os.remove(tmp_filepath)


def test_file_sdf_input_output():
    """Validity of SDF file input and output."""
    tmp_filepath = tempfile.NamedTemporaryFile(suffix='.sdf').name
    args = arg_parser(["-i", f"{examples_dir}/payload_4.sdf",
                       "-if", "sdf",
                       "-o", tmp_filepath,
                       "-of", "sdf", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    run_cli(args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(f"{examples_dir}/payload_4_out.sdf", "r") as f:
        expected = f.read()
        with open(tmp_filepath, "r") as f:
            results = f.read()
            assert results == expected
    os.remove(tmp_filepath)


def test_smiles_stdin_input():
    """Validity of SMILES stdin input and text output."""
    args = arg_parser(["-i", "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)"
                             "(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H]"
                             "(N)Cc1ccccc1)C(=O)O",
                       "-if", "smiles_stdin",
                       "--print_fragment_pkas", "--method", "pkamatcher"])
    result = stdout_to_variable(run_cli, args)
    with open(f"{examples_dir}/payload_1_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_file_ph_q_plot_1():
    """Validity of the pH/Q plot file."""
    tmp_file_prefix = tempfile.NamedTemporaryFile().name
    tmp_filepath = f"{tmp_file_prefix}_1.png"
    args = arg_parser(["-i", f"{examples_dir}/payload_1.smi",
                       "-of", "json", "--print_fragment_pkas",
                       "-pp", tmp_file_prefix,
                       "--plot_ph_q_curve",
                       "--method", "pkamatcher"])
    _ = stdout_to_variable(run_cli, args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(tmp_filepath, "rb") as f:
        img_txt = f.read()
    os.remove(tmp_filepath)
    base64_hash = base64.b64encode(img_txt)
    hash_object = hashlib.sha224(base64_hash)
    hex_dig = hash_object.hexdigest()
    expected = "6e612de5b63edc90da7b090bdf209f74f096a6feb9a0167385d425b5"
    assert hex_dig == expected


def test_file_ph_q_plot_2():
    """Validity of the pH/Q plot file."""
    tmp_file_prefix = tempfile.NamedTemporaryFile().name
    tmp_filepath = f"{tmp_file_prefix}_1.png"
    args = arg_parser(["-i", f"{examples_dir}/payload_2.smi",
                       "-of", "json", "--print_fragment_pkas",
                       "-pp", tmp_file_prefix,
                       "--plot_ph_q_curve",
                       "--method", "pkamatcher"])
    _ = stdout_to_variable(run_cli, args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(tmp_filepath, "rb") as f:
        img_txt = f.read()
    os.remove(tmp_filepath)
    base64_hash = base64.b64encode(img_txt)
    hash_object = hashlib.sha224(base64_hash)
    hex_dig = hash_object.hexdigest()
    expected = "33ae664e79247014266e751b2ac9a013d2f974c21559dd1ab8daeb7e"
    assert hex_dig == expected
