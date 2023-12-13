from rdkit import Chem


def _pattern_match_rdkit(smiles, smarts):
    """RDKit SMARTS matching counter."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    pat = Chem.MolFromSmarts(smarts)
    n = 0
    for _ in mol.GetSubstructMatches(pat):
        n += 1
    return n


def pattern_match(smiles, smarts):
    """Interface for SMARTS matching."""
    return _pattern_match_rdkit(smiles, smarts)
