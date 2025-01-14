import json
import os
import tempfile

import pytest
from helpers import examples_dir
from helpers import stdout_to_variable
from helpers import TestError
from pichemist.cli import arg_parser
from pichemist.cli import run_pichemist


def test_parser_creation():
    """Parser creation."""
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_1.smi",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    assert args is not None


def test_console_text_output_1():
    """Validity of console text output using pKaMatcher."""
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_1.smi",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_1_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_console_text_output_2():
    """Validity of console text output using pKaMatcher."""
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_3.smi",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_3_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


@pytest.mark.acd
def test_console_text_output_3():
    """Validity of console text output using ACD."""
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_5.smi",
            "--print_fragment_pkas",
            "--method",
            "acd",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_5_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_console_json_output():
    """Validity of console JSON output."""
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_1.smi",
            "-of",
            "json",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    result = json.loads(result)
    with open(f"{examples_dir}/payload_1_out.json", "r") as f:
        expected = json.load(f)
    assert result == expected


def test_console_json_output_images():
    """Validity of console JSON output."""
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_1.smi",
            "-of",
            "json",
            "--print_fragment_pkas",
            "--generate_fragment_images",
            "--method",
            "pkamatcher",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    result = json.loads(result)
    with open(f"{examples_dir}/payload_1_out_images.json", "r") as f:
        expected = json.load(f)
    assert result == expected


def test_file_csv_output():
    """Validity of CSV file output."""
    tmp_filepath = tempfile.NamedTemporaryFile(suffix=".csv").name
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_1.smi",
            "-o",
            tmp_filepath,
            "-of",
            "csv",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    run_pichemist(args)
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
    tmp_filepath = tempfile.NamedTemporaryFile(suffix=".sdf").name
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_1.smi",
            "-o",
            tmp_filepath,
            "-of",
            "sdf",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    run_pichemist(args)
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
    tmp_filepath = tempfile.NamedTemporaryFile(suffix=".sdf").name
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_3.smi",
            "-o",
            tmp_filepath,
            "-of",
            "sdf",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    run_pichemist(args)
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
    tmp_filepath = tempfile.NamedTemporaryFile(suffix=".sdf").name
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_4.sdf",
            "-if",
            "sdf",
            "-o",
            tmp_filepath,
            "-of",
            "sdf",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    run_pichemist(args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(f"{examples_dir}/payload_4_out.sdf", "r") as f:
        expected = f.read()
        with open(tmp_filepath, "r") as f:
            results = f.read()
            assert results == expected
    os.remove(tmp_filepath)


def test_smiles_stdin_input_1():
    """Validity of SMILES stdin input and text output."""
    args = arg_parser(
        [
            "-i",
            "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)"
            "(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H]"
            "(N)Cc1ccccc1)C(=O)O",
            "-if",
            "smiles_stdin",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_1_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_smiles_stdin_input_2():
    """Validity of SMILES stdin input and text output."""
    args = arg_parser(
        [
            "-i",
            "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=O)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=O)N)C(=O)O",  # noqa: E501
            "-if",
            "smiles_stdin",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_2_out_smiles.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_fasta_stdin_input():
    """Validity of FASTA stdin input and text output."""
    args = arg_parser(
        [
            "-i",
            "CNCN",
            "-if",
            "fasta_stdin",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_2_out_fasta.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_fasta_stdin_capped_cterm_input():
    """Validity of FASTA stdin input and text output."""
    args = arg_parser(
        [
            "-i",
            "CNCN",
            "-if",
            "fasta_stdin",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
            "--ionizable_cterm",
            "false",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_2_out_fasta_capped_cterm.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_fasta_stdin_capped_nterm_input():
    """Validity of FASTA stdin input and text output."""
    args = arg_parser(
        [
            "-i",
            "CNCN",
            "-if",
            "fasta_stdin",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
            "--ionizable_nterm",
            "false",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_2_out_fasta_capped_nterm.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_fasta_stdin_all_capped_input():
    """Validity of FASTA stdin input and text output."""
    args = arg_parser(
        [
            "-i",
            "CNCN",
            "-if",
            "fasta_stdin",
            "--print_fragment_pkas",
            "--method",
            "pkamatcher",
            "--ionizable_nterm",
            "false",
            "--ionizable_cterm",
            "false",
        ]
    )
    result = stdout_to_variable(run_pichemist, args)
    with open(f"{examples_dir}/payload_2_out_fasta_all_capped.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_file_ph_q_plot_1():
    """Existence of the pH/Q plot file."""
    tmp_file_prefix = tempfile.NamedTemporaryFile().name
    tmp_filepath = f"{tmp_file_prefix}_1.png"
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_1.smi",
            "-of",
            "json",
            "--print_fragment_pkas",
            "-pp",
            tmp_file_prefix,
            "--plot_ph_q_curve",
            "--method",
            "pkamatcher",
        ]
    )
    _ = stdout_to_variable(run_pichemist, args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    os.remove(tmp_filepath)


def test_file_ph_q_plot_2():
    """Existence of the pH/Q plot file."""
    tmp_file_prefix = tempfile.NamedTemporaryFile().name
    tmp_filepath = f"{tmp_file_prefix}_1.png"
    args = arg_parser(
        [
            "-i",
            f"{examples_dir}/payload_2.smi",
            "-of",
            "json",
            "--print_fragment_pkas",
            "-pp",
            tmp_file_prefix,
            "--plot_ph_q_curve",
            "--method",
            "pkamatcher",
        ]
    )
    _ = stdout_to_variable(run_pichemist, args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    os.remove(tmp_filepath)
