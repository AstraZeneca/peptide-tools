import tempfile

import pytest
from helpers import examples_dir
from pichemist.api import pichemist_from_list
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
        1: {
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
                "pI mean": 9.0234375,
                "std": 1.721588565104915,
                "err": 0.6086734743994516,
            },
            "QpH7": {
                "IPC2_peptide": 0.6314906212267486,
                "IPC_peptide": 0.9915539516610472,
                "ProMoST": 0.26174063515548607,
                "Gauci": 0.5540630760817584,
                "Grimsley": 0.6645409545014482,
                "Thurlkill": 0.797542965316429,
                "Lehninger": 0.9932283675959863,
                "Toseland": 0.9515959465104951,
                "Q at pH7.4 mean": 0.7307195647561748,
                "std": 0.6749606913955383,
                "err": 0.23863464096007284,
            },
            "pI_interval": (8.624999999999998, 9.362499999999997),
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_list(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == output


def test_pka_matcher_json_1_stdin():
    """Example with mixed amino acids using pKaMatcher."""
    args = {
        "input_data": "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O",  # ignore
        "input_format": "smiles_stdin",
        "plot_ph_q_curve": False,
        "print_fragments": False,
        "method": "pkamatcher",
    }

    expected = {
        1: {
            "mol_name": "C[C@@H](NC(=O)[C@H](CCCCN)NC(=O)[C@](C)(CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O",  # ignore
            "pI": {
                "IPC2_peptide": 8.046875,
                "IPC_peptide": 9.8125,
                "ProMoST": 8.375,
                "Gauci": 8.6875,
                "Grimsley": 8.9375,
                "Thurlkill": 9.0625,
                "Lehninger": 9.859375,
                "Toseland": 9.40625,
                "pI mean": 9.0234375,
                "std": 1.721588565104915,
                "err": 0.6086734743994516,
            },
            "QpH7": {
                "IPC2_peptide": 0.6314906212267486,
                "IPC_peptide": 0.9915539516610472,
                "ProMoST": 0.26174063515548607,
                "Gauci": 0.5540630760817584,
                "Grimsley": 0.6645409545014482,
                "Thurlkill": 0.797542965316429,
                "Lehninger": 0.9932283675959863,
                "Toseland": 0.9515959465104951,
                "Q at pH7.4 mean": 0.7307195647561748,
                "std": 0.6749606913955383,
                "err": 0.23863464096007284,
            },
            "pI_interval": (8.624999999999998, 9.362499999999997),
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_list(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == output


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
        1: {
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
                "pI mean": 5.09375,
                "std": 0.6789237807000135,
                "err": 0.240035804620894,
            },
            "QpH7": {
                "IPC2_peptide": -0.23913646012640216,
                "IPC_peptide": -0.23180587970079417,
                "ProMoST": -0.47710727633330485,
                "Gauci": -0.9370497315178979,
                "Grimsley": -1.9322611310762747,
                "Thurlkill": -0.33280195754105635,
                "Lehninger": -0.21536815651977892,
                "Toseland": -1.5908867442915533,
                "Q at pH7.4 mean": -0.7445521671383828,
                "std": 1.7898169252151255,
                "err": 0.6327958424510355,
            },
            "pI_interval": (3.7624999999999984, 6.687499999999999),
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_list(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == output


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
        1: {
            "mol_name": "Phe-Ornithine-aMeAsp-Lys-dAla",
            "pI": {
                "IPC2_peptide": 8.0625,
                "IPC_peptide": 10.03125,
                "ProMoST": 8.375,
                "Gauci": 8.75,
                "Grimsley": 9.125,
                "Thurlkill": 9.1875,
                "Lehninger": 10.09375,
                "Toseland": 9.5625,
                "pI mean": 9.1484375,
                "std": 1.9449202666819019,
                "err": 0.6876331547189606,
            },
            "QpH7": {
                "IPC2_peptide": 0.6323748200446536,
                "IPC_peptide": 0.9924381504789521,
                "ProMoST": 0.2626248339733912,
                "Gauci": 0.5549472748996638,
                "Grimsley": 0.665425153319353,
                "Thurlkill": 0.7984271641343339,
                "Lehninger": 0.9941125664138912,
                "Toseland": 0.9524801453284001,
                "Q at pH7.4 mean": 0.7316037635740799,
                "std": 0.6749606913955382,
                "err": 0.2386346409600728,
            },
            "pI_interval": (8.687499999999998, 9.612499999999997),
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_list(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == output


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
        1: {
            "mol_name": "Galas789",
            "pI": {
                "Gauci": 4.75,
                "Grimsley": 5.0,
                "IPC2_peptide": 4.75,
                "IPC_peptide": 4.25,
                "Lehninger": 4.25,
                "ProMoST": 5.0,
                "Thurlkill": 5.125,
                "Toseland": 4.75,
                "err": 0.30896437395758103,
                "pI mean": 4.734375,
                "std": 0.8738832158818477,
            },
            "QpH7": {
                "Gauci": -0.9632323250192959,
                "Grimsley": -0.8631136848823902,
                "IPC2_peptide": -0.9013539271954132,
                "IPC_peptide": -0.9601482453547936,
                "Lehninger": -0.9617047865577028,
                "ProMoST": -0.7638658182464676,
                "Q at pH7.4 mean": -0.9016797792224913,
                "Thurlkill": -0.8785187771654149,
                "Toseland": -0.9215006693584518,
                "err": 0.06346230799231138,
                "std": 0.17949851332445046,
            },
            "pI_interval": (3.6999999999999984, 5.724999999999998),
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_list(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == output


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
        1: {
            "mol_name": "Cys-Asn-Cys-Asn",
            "pI": {
                "IPC2_peptide": 5,
                "IPC_peptide": 5,
                "ProMoST": 5.5,
                "Gauci": 5,
                "Grimsley": 4.875,
                "Thurlkill": 5.5,
                "Lehninger": 5,
                "Toseland": 4.875,
                "pI mean": 5.09375,
                "std": 0.6789237807000135,
                "err": 0.240035804620894,
            },
            "QpH7": {
                "IPC2_peptide": -0.23913646012640216,
                "IPC_peptide": -0.23180587970079417,
                "ProMoST": -0.47710727633330485,
                "Gauci": -0.9370497315178979,
                "Grimsley": -1.9322611310762747,
                "Thurlkill": -0.33280195754105635,
                "Lehninger": -0.21536815651977892,
                "Toseland": -1.5908867442915533,
                "Q at pH7.4 mean": -0.7445521671383828,
                "std": 1.7898169252151255,
                "err": 0.6327958424510355,
            },
            "pI_interval": (3.7624999999999984, 6.687499999999999),
            "pI_interval_threshold": 0.2,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_list(
        input_dict, args["method"], args["plot_ph_q_curve"], args["print_fragments"]
    )
    assert expected == output


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
        1: {
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
                "pI mean": 5.09375,
                "std": 0.6789237807000135,
                "err": 0.240035804620894,
            },
            "QpH7": {
                "IPC2_peptide": -0.23913646012640216,
                "IPC_peptide": -0.23180587970079417,
                "ProMoST": -0.47710727633330485,
                "Gauci": -0.9370497315178979,
                "Grimsley": -1.9322611310762747,
                "Thurlkill": -0.33280195754105635,
                "Lehninger": -0.21536815651977892,
                "Toseland": -1.5908867442915533,
                "Q at pH7.4 mean": -0.7445521671383828,
                "std": 1.7898169252151255,
                "err": 0.6327958424510355,
            },
            "pI_interval": (3.7624999999999984, 6.687499999999999),
            "pI_interval_threshold": 0.2,
            "plot_filename": tmp_filepath,
            "pKa_set": "IPC2_peptide",
        }
    }

    input_dict = generate_input(args["input_format"], args["input_data"])
    output = pichemist_from_list(
        input_dict,
        args["method"],
        args["ph_q_curve_file_prefix"],
        args["plot_ph_q_curve"],
        args["print_fragments"],
    )
    assert expected == output
