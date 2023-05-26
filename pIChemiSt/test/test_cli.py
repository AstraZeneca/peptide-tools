import os

from pichemist.io import generate_input
from pichemist.cli import calc_pichemist

script_dir = os.path.dirname(os.path.realpath(__file__))


def test_pka_matcher_json():
    args = {"input": f"{script_dir}/example_pka_matcher_1.smi",
            "input_format": "smiles_file",
            "plot_titration_curve": False,
            "print_fragments": False,
            "print_pka_set": False,
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
    output = calc_pichemist(input_dict, args["method"],
                            args["plot_titration_curve"],
                            args["print_fragments"],
                            args["print_pka_set"])
    assert expected == output