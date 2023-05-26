from rdkit import Chem
from pichemist.molecule import MolStandardiser

def calc_net_Qs(smi_list):
    net_Qs = []
    for smi in smi_list:
        mol = Chem.MolFromSmiles(smi)
        mol = MolStandardiser().standardise_molecule(mol)

        # positively charged N, but not connected to any negatively charged atom. To avoid azido and nitro groups being counted.
        pattern = Chem.MolFromSmarts("[#7+;!H1;!H2;!H3;!H4!$([#7+]~[*-])]")
        at_matches = mol.GetSubstructMatches(pattern)

        at_matches_list = [y[0] for y in at_matches]

        for _ in at_matches_list:
            net_Qs.append( (1,smi) )
    return net_Qs