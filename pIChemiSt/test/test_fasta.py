from pichemist.fasta.scrambler import FastaScrambler

scrambler = FastaScrambler()


def test_fasta_scrambler_list():
    """Produces a FASTA sequence from a set of amino acids."""
    smiles = [
        "CC(=O)N[C@@H](C)C(C)=O",
        "CC(=O)N[C@@H](CCC(=O)O)C(=O)O",
        "CC(=O)N[C@@H](Cc1ccc(O)cc1)C(C)=O",
        "CC(=O)N[C@H](C(C)=O)C(C)C",
        "CC(=O)[C@@H](N)Cc1ccccc1",
        "CC(=O)[C@@H]1CCCN1C(C)=O"
    ]
    res = scrambler.get_scrambled_fasta_from_list(smiles)
    expected = "AEYVFP"
    assert res == expected, f"got {res}"
