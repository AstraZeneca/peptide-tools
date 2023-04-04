from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# https://www.rdkit.org/docs/Cookbook.html
# Replace all ionized centers by the corresponsing neutral form. Used to track constanty ionized fragments, like quaternary amines

def _neutralize_atoms(mol):
    """
    Replace all ionized centers by the corresponsing
    neutral form. Used to track constanty ionized
    fragments, like quaternary amines.
    https://www.rdkit.org/docs/Cookbook.html

    """
    ION_SMARTS = "[+1!h0!$([*]~[-1,-2,-3,-4])" \
                 ",-1!$([*]~[+1,+2,+3,+4])]"

    pattern = Chem.MolFromSmarts(ION_SMARTS)
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def standardize_molecule(mol):
    """
    Runs molecule standardisation using a built-in
    method and some SMARTS-based heuristics.

    """
    clean_mol = rdMolStandardize.Cleanup(mol)
    m = _neutralize_atoms(clean_mol)
    return m


def fragmentMol(mol, breakableBonds, atomNum):
    # Function from breaking and capping
    readwritemol = Chem.RWMol(mol)
    for bond in breakableBonds: 
        atm1 = bond[0]
        atm2 = bond[1]
        readwritemol.RemoveBond(atm1, atm2)
        amineC = readwritemol.AddAtom(Chem.Atom(6))
        amineCO = readwritemol.AddAtom(Chem.Atom(8))
        amineCOC = readwritemol.AddAtom(Chem.Atom(6))
        readwritemol.AddBond(atm1, amineC, Chem.BondType.SINGLE) 
        readwritemol.AddBond(amineC, amineCO, Chem.BondType.DOUBLE) 
        readwritemol.AddBond(amineC, amineCOC, Chem.BondType.SINGLE) 
        newatom2 = readwritemol.AddAtom(Chem.Atom(atomNum))
        readwritemol.AddBond(atm2, newatom2, Chem.BondType.SINGLE) 
        fragmentedSmiles = Chem.MolToSmiles(readwritemol)
    return fragmentedSmiles

def break_amide_bonds_and_cap(mol):
    ### SMARTS pattern
    amideSMARTS = '[NX3,NX4;H0,H1][CX3](=[OX1])' # secondary and tertiary amide bonds
    amidePattern = Chem.MolFromSmarts(amideSMARTS)
    atomNum = 6

    breakableBonds = []
    #molname=mol.GetProp("_Name")
    smiles = Chem.MolToSmiles(mol)
    if mol.HasSubstructMatch(amidePattern):
        #smiles = Chem.MolToSmiles(mol)
        #results.write("%s," % molname)
        atomIDs = mol.GetSubstructMatches(amidePattern)
        for bond in atomIDs:
            breakableBonds.append((bond[0],bond[1]))
        fragmentedSmiles = fragmentMol(mol, breakableBonds, atomNum)
    else:
        fragmentedSmiles = smiles

    fragmentedSmilesList = fragmentedSmiles.split(".")
    numFrags = len(fragmentedSmilesList)
    #print("%s fragmented into %s pieces" % (molname,numFrags))
    return fragmentedSmilesList