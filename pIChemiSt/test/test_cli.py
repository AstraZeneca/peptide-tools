import json
import os
import tempfile

from helpers import TestError
from helpers import examples_dir
from helpers import stdout_to_variable
from pichemist.cli import arg_parser
from pichemist.cli import run_cli


def test_parser_creation():
    """Parser creation."""
    args = arg_parser(["-i", f"{examples_dir}/pka_matcher_1.smi",
                       "--print_fragment_pkas", "--method", "pkamatcher"])
    assert args is not None


def test_console_text_output():
    """Validity of console text output."""
    args = arg_parser(["-i", f"{examples_dir}/pka_matcher_1.smi",
                       "--print_fragment_pkas", "--method", "pkamatcher"])
    result = stdout_to_variable(run_cli, args)
    with open(f"{examples_dir}/pka_matcher_1_out.txt", "r") as f:
        expected = f.read()
    assert result == expected


def test_console_json_output():
    """Validity of console JSON output."""
    args = arg_parser(["-i", f"{examples_dir}/pka_matcher_1.smi",
                       "-of", "json", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    result = stdout_to_variable(run_cli, args)
    result = json.loads(result)
    with open(f"{examples_dir}/pka_matcher_1_out.json", "r") as f:
        expected = json.load(f)
    assert result == expected


def test_file_csv_output():
    """Validity of CSV file output."""
    tmp_filepath = tempfile.NamedTemporaryFile(suffix='.csv').name
    args = arg_parser(["-i", f"{examples_dir}/pka_matcher_1.smi",
                       "-o", tmp_filepath,
                       "-of", "csv", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    run_cli(args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(f"{examples_dir}/pka_matcher_1_out.csv", "r") as f:
        expected = f.read()
        with open(tmp_filepath, "r") as f:
            results = f.read()
            assert results == expected
    os.remove(tmp_filepath)


def test_file_sdf_output():
    """Validity of SDF file output."""
    tmp_filepath = tempfile.NamedTemporaryFile(suffix='.sdf').name
    args = arg_parser(["-i", f"{examples_dir}/pka_matcher_1.smi",
                       "-o", tmp_filepath,
                       "-of", "sdf", "--print_fragment_pkas",
                       "--method", "pkamatcher"])
    run_cli(args)
    if not os.path.exists(tmp_filepath):
        raise TestError("File was not created.")
    with open(f"{examples_dir}/pka_matcher_1_out.sdf", "r") as f:
        expected = f.read()
        with open(tmp_filepath, "r") as f:
            results = f.read()
            assert results == expected
