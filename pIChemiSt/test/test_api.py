from helpers import examples_dir
from pichemist.io import generate_input
from pichemist.api import pichemist_from_list


def test_pka_matcher_json_1():
    """Example with mixed amino acids using pKaMatcher."""
    args = {"input": f"{examples_dir}/pka_matcher_1.smi",
            "input_format": "smiles_file",
            "plot_titration_curve": False,
            "print_fragments": False,
            "method": "pkamatcher"}

    expected = {1:
                {"mol_name": "Phe-Ornithine-aMeAsp-Lys-dAla",
                 "pI":
                    {
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
                        "err": 0.6086734743994516},
                 "QpH7":
                    {
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
                        "err": 0.23863464096007284
                    },
                 "pI_interval": (8.624999999999998, 9.362499999999997),
                 "plot_filename": "",
                 "pI_interval_threshold": 0.2,
                 "pKa_set": "IPC2_peptide"}
                }

    input_dict = generate_input(args["input_format"], args["input"])
    output = pichemist_from_list(input_dict, args["method"],
                                 args["plot_titration_curve"],
                                 args["print_fragments"])
    assert expected == output


def test_pka_matcher_json_2():
    """Example with only natural amino acids (only FASTA matching)."""
    args = {"input": f"{examples_dir}/pka_matcher_2.smi",
            "input_format": "smiles_file",
            "plot_titration_curve": False,
            "print_fragments": False,
            "method": "pkamatcher"}

    expected = {1:
                {"mol_name": "Cys-Asn-Cys-Asn",
                 "pI":
                    {
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
                        "err": 0.240035804620894},
                 "QpH7":
                    {
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
                        "err": 0.6327958424510355
                    },
                 "pI_interval": (3.7624999999999984, 6.687499999999999),
                 "plot_filename": "",
                 "pI_interval_threshold": 0.2,
                 "pKa_set": "IPC2_peptide"}
                }

    input_dict = generate_input(args["input_format"], args["input"])
    output = pichemist_from_list(input_dict, args["method"],
                                 args["plot_titration_curve"],
                                 args["print_fragments"])
    assert expected == output