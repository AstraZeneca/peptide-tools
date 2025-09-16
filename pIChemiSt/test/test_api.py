import tempfile

import pytest
from helpers import examples_dir
from helpers import jsonify_output
from pichemist.api import pichemist_from_dict
from pichemist.io import generate_input


def test_pka_matcher_json_1_file():
    """Example with mixed amino acids using pKaMatcher."""
    args = {
        "input_data": f"{examples_dir}/payload_1.smi",
        "input_format": "smiles_file",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "pkamatcher",
    }

    expected = {
        "1": {
            "mol_name": "Phe-Ornithine-aMeAsp-Lys-dAla",
            "pI": {
                "IPC2_peptide": 8.046875,
                "IPC_peptide": 9.8125,
                "ProMoST": 8.375,
                "Gauci": 8.6875,
                "Grimsley": 8.9375,
                "Thurlkill": 9.0625,
                "Lehninger": 9.859375,
                "Toseland": 9.40625,
                "pI mean": 9.0234,
                "std": 1.7216,
                "err": 0.6087,
            },
            "QpH7": {
                "IPC2_peptide": 0.6315,
                "IPC_peptide": 0.9916,
                "ProMoST": 0.2617,
                "Gauci": 0.5541,
                "Grimsley": 0.6645,
                "Thurlkill": 0.7975,
                "Lehninger": 0.9932,
                "Toseland": 0.9516,
                "Q at pH7.4 mean": 0.7307,
                "std": 0.675,
                "err": 0.2386,
            },
            "pI_interval": [8.625, 9.3625],
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_dict(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == jsonify_output(output)


def test_pka_matcher_json_1_stdin():
    """Example with mixed amino acids using pKaMatcher."""
    args = {
        "input_data": "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O",  # noqa
        "input_format": "smiles_stdin",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "pkamatcher",
    }

    expected = {
        "1": {
            "mol_name": "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O",  # noqa
            "pI": {
                "IPC2_peptide": 8.046875,
                "IPC_peptide": 9.8125,
                "ProMoST": 8.375,
                "Gauci": 8.6875,
                "Grimsley": 8.9375,
                "Thurlkill": 9.0625,
                "Lehninger": 9.859375,
                "Toseland": 9.40625,
                "pI mean": 9.0234,
                "std": 1.7216,
                "err": 0.6087,
            },
            "QpH7": {
                "IPC2_peptide": 0.6315,
                "IPC_peptide": 0.9916,
                "ProMoST": 0.2617,
                "Gauci": 0.5541,
                "Grimsley": 0.6645,
                "Thurlkill": 0.7975,
                "Lehninger": 0.9932,
                "Toseland": 0.9516,
                "Q at pH7.4 mean": 0.7307,
                "std": 0.675,
                "err": 0.2386,
            },
            "pI_interval": [8.625, 9.3625],
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_dict(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == jsonify_output(output)


def test_natural_aa_json_1():
    """Example with only natural amino acids (only FASTA matching)."""
    args = {
        "input_data": f"{examples_dir}/payload_2.smi",
        "input_format": "smiles_file",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "pkamatcher",
    }

    expected = {
        "1": {
            "mol_name": "Cys-Asn-Cys-Asn",
            "pI": {
                "IPC2_peptide": 5.0,
                "IPC_peptide": 5.0,
                "ProMoST": 5.5,
                "Gauci": 5.0,
                "Grimsley": 4.875,
                "Thurlkill": 5.5,
                "Lehninger": 5.0,
                "Toseland": 4.875,
                "pI mean": 5.0938,
                "std": 0.6789,
                "err": 0.24,
            },
            "QpH7": {
                "IPC2_peptide": -0.2391,
                "IPC_peptide": -0.2318,
                "ProMoST": -0.4771,
                "Gauci": -0.937,
                "Grimsley": -1.9323,
                "Thurlkill": -0.3328,
                "Lehninger": -0.2154,
                "Toseland": -1.5909,
                "Q at pH7.4 mean": -0.7446,
                "std": 1.7898,
                "err": 0.6328,
            },
            "pI_interval": [3.7625, 6.6875],
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_dict(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == jsonify_output(output)


@pytest.mark.acd
def test_acd_json_1():
    """Example with mixed amino acids using ACD."""
    args = {
        "input_data": f"{examples_dir}/payload_1.smi",
        "input_format": "smiles_file",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "acd",
    }

    expected = {
        "1": {
            "mol_name": "Phe-Ornithine-aMeAsp-Lys-dAla",
            "pI": {
                "IPC2_peptide": 8.046875,
                "IPC_peptide": 9.75,
                "ProMoST": 8.3125,
                "Gauci": 8.625,
                "Grimsley": 8.875,
                "Thurlkill": 9.0,
                "Lehninger": 9.796875,
                "Toseland": 9.34375,
                "pI mean": 8.9688,
                "std": 1.6868,
                "err": 0.5964,
            },
            "QpH7": {
                "IPC2_peptide": 0.6313,
                "IPC_peptide": 0.9914,
                "ProMoST": 0.2616,
                "Gauci": 0.5539,
                "Grimsley": 0.6644,
                "Thurlkill": 0.7974,
                "Lehninger": 0.9931,
                "Toseland": 0.9514,
                "Q at pH7.4 mean": 0.7306,
                "std": 0.675,
                "err": 0.2386,
            },
            "pI_interval": [8.6125, 9.3125],
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_dict(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == jsonify_output(output)


@pytest.mark.acd
def test_acd_json_2():
    """Example with mixed amino acids using ACD (GALAS)."""
    args = {
        "input_data": f"{examples_dir}/payload_6.smi",
        "input_format": "smiles_file",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "acd",
    }

    expected = {
        "1": {
            "mol_name": "Galas789",
            "pI": {
                "IPC2_peptide": 7.875,
                "IPC_peptide": 7.75,
                "ProMoST": 8.125,
                "Gauci": 7.75,
                "Grimsley": 8.0,
                "Thurlkill": 8.0,
                "Lehninger": 7.75,
                "Toseland": 7.875,
                "pI mean": 7.8906,
                "std": 0.3724,
                "err": 0.1317,
            },
            "QpH7": {
                "IPC2_peptide": 0.0887,
                "IPC_peptide": 0.03,
                "ProMoST": 0.2262,
                "Gauci": 0.0269,
                "Grimsley": 0.127,
                "Thurlkill": 0.1116,
                "Lehninger": 0.0284,
                "Toseland": 0.0686,
                "Q at pH7.4 mean": 0.0884,
                "std": 0.1794,
                "err": 0.0634,
            },
            "pI_interval": [6.9875, 8.8],
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_dict(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == jsonify_output(output)


@pytest.mark.acd
def test_natural_aa_json_2():
    """Example with only natural amino acids (only FASTA matching)."""
    args = {
        "input_data": f"{examples_dir}/payload_2.smi",
        "input_format": "smiles_file",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "acd",
    }

    expected = {
        "1": {
            "mol_name": "Cys-Asn-Cys-Asn",
            "pI": {
                "IPC2_peptide": 5.0,
                "IPC_peptide": 5.0,
                "ProMoST": 5.5,
                "Gauci": 5.0,
                "Grimsley": 4.875,
                "Thurlkill": 5.5,
                "Lehninger": 5.0,
                "Toseland": 4.875,
                "pI mean": 5.0938,
                "std": 0.6789,
                "err": 0.24,
            },
            "QpH7": {
                "IPC2_peptide": -0.2391,
                "IPC_peptide": -0.2318,
                "ProMoST": -0.4771,
                "Gauci": -0.937,
                "Grimsley": -1.9323,
                "Thurlkill": -0.3328,
                "Lehninger": -0.2154,
                "Toseland": -1.5909,
                "Q at pH7.4 mean": -0.7446,
                "std": 1.7898,
                "err": 0.6328,
            },
            "pI_interval": [3.7625, 6.6875],
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_dict(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == jsonify_output(output)


def test_natural_aa_json_3():
    """
    Example with only natural amino acids (only FASTA matching).
    It also checks for the presence of "plot_filename" in the output.

    """
    tmp_file_prefix = tempfile.NamedTemporaryFile().name
    tmp_filepath = f"{tmp_file_prefix}_1.png"
    args = {
        "input_data": f"{examples_dir}/payload_2.smi",
        "input_format": "smiles_file",
        "ph_q_curve_file_prefix": tmp_file_prefix,
        "plot_ph_q_curve": True,
        "print_fragments": False,
        "method": "pkamatcher",
    }

    expected = {
        "1": {
            "mol_name": "Cys-Asn-Cys-Asn",
            "pI": {
                "IPC2_peptide": 5.0,
                "IPC_peptide": 5.0,
                "ProMoST": 5.5,
                "Gauci": 5.0,
                "Grimsley": 4.875,
                "Thurlkill": 5.5,
                "Lehninger": 5.0,
                "Toseland": 4.875,
                "pI mean": 5.0938,
                "std": 0.6789,
                "err": 0.24,
            },
            "QpH7": {
                "IPC2_peptide": -0.2391,
                "IPC_peptide": -0.2318,
                "ProMoST": -0.4771,
                "Gauci": -0.937,
                "Grimsley": -1.9323,
                "Thurlkill": -0.3328,
                "Lehninger": -0.2154,
                "Toseland": -1.5909,
                "Q at pH7.4 mean": -0.7446,
                "std": 1.7898,
                "err": 0.6328,
            },
            "pI_interval": [3.7625, 6.6875],
            "pI_interval_threshold": 0.2,
            "plot_filename": tmp_filepath,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_dict(
        input_dict,
        args["method"],
        args["ph_q_curve_file_prefix"],
        args["plot_ph_q_curve"],
        args["print_fragments"],
    )
    assert expected == jsonify_output(output)
