import json

from helpers import examples_dir
from pichemist.fasta.smarts import AA_SMARTS_SET
from pichemist.fasta.pka_sets import FASTA_PKA_SETS
from pichemist.smarts.pka_set import SS_SMARTS_PKA_SET


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
