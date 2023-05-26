from pichemist.charges import ChargeCalculator

qcalc = ChargeCalculator()


def test_charges_smiles_1():
    res = qcalc.calculate_net_qs_from_smiles("CCCC[N+](C)(C)[N+]C")
    expected = 2
    assert res == expected, f"got {res}"


def test_charges_smiles_2():
    res = qcalc.calculate_net_qs_from_smiles("CCCC[N+](C)(C)[N-]C")
    expected = 0
    assert res == expected, f"got {res}"

# TODO: ANDREY: Need examples to test this
# TODO: Verify with Andrey
# def test_charges_smiles_3():
#     res = qcalc.calculate_net_qs_from_smiles("CCCC[N+](C)(C)C[B-]")
#     expected = 0    # Should not this be 0?
#     assert res == expected, f"got {res}"


def test_charges_smiles_4():
    res = qcalc.calculate_net_qs_from_list(["CCCC[N+](C)(C)[N+]C"])
    expected = [(1, 'CCCC[N+](C)(C)[N+]C'), (1, 'CCCC[N+](C)(C)[N+]C')]
    assert res == expected, f"got {res}"
