import os
import json

from pichemist import cli
from pichemist.fasta.smarts import AA_SMARTS_SET
from pichemist.fasta.pka_sets import FASTA_PKA_SETS
from pichemist.smarts.pka_set import SS_SMARTS_PKA_SET

script_dir = os.path.dirname(os.path.realpath(__file__))
examples_dir = os.path.join(script_dir, "examples")


def test_pka_matcher_json():
    options = {"inputJSON": "",
               "inputFile": f"{script_dir}/example_pka_matcher_1.smi",
               "smiles": "",
               "outputFile": "",
               "l_plot_titration_curve": False,
               "l_print_fragments": False,
               "l_print_pka_set": False,
               "use_acdlabs": False,
               "use_pkamatcher": True,
               "l_json": True}

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

    output = cli.calc_pIChemiSt(options)
    # print(output)
    assert expected == output


def test_fasta_pka_sets():
    with open(f"{examples_dir}/fasta_pka_sets_expected.json") as f:
        expected = json.load(f)
    assert expected == FASTA_PKA_SETS


def test_ss_smarts_pka_set():
    with open(f"{examples_dir}/ss_smarts_pka_set_expected.json") as f:
        expected = json.load(f)
    assert expected == SS_SMARTS_PKA_SET


def test_aa_smarts_pka_set():
    with open(f"{examples_dir}/aa_smarts_pka_set_expected.json") as f:
        expected = json.load(f)
    assert expected == AA_SMARTS_SET
