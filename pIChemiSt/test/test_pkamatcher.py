from pichemist.pka.pkamatcher import PKaMatcher

matcher = PKaMatcher()


def test_pka_matcher_list():
    smiles = ["CC(=O)N[C@@H](CCCN)C(C)=O", "CC(=O)N[C@@](C)(CC(=O)O)C(C)=O"]
    res = matcher.calculate_pka_from_list(smiles)
    expected = (
        [(10.4, "CC(=O)N[C@@H](CCCN)C(C)=O")],
        [(3.46, "CC(=O)N[C@@](C)(CC(=O)O)C(C)=O")],
    )
    assert res == expected, f"got {res}"


def test_pka_matcher_internal_smiles():
    fragment = "CC(=O)N[C@@H](CCCN)C(C)=O"
    res = matcher._pka_dict_from_smiles(fragment)
    expected = {"acid": [], "base": [(10.4, "CC(=O)N[C@@H](CCCN)C(C)=O")]}
    assert res == expected, f"got {res}"


def test_pka_matcher_smiles():
    fragment = "CC(=O)N[C@@H](CCCN)C(C)=O"
    res = matcher.calculate_pka_from_smiles(fragment)
    expected = {"acid": [], "base": [10.4]}
    assert res == expected, f"got {res}"
