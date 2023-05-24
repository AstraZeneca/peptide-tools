from enum import Enum
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from pichemist.utils import get_logger

log = get_logger(__name__)


class Standardiser(object):
    """Deals with molecule standardisation."""

    def __init__(self):
        pass

    def neutralise_molecule(self, mol):
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

    def standardise_molecule(self, mol):
        """
        Runs molecule standardisation using a built-in
        method and some SMARTS-based heuristics.

        """
        clean_mol = rdMolStandardize.Cleanup(mol)
        m = self.neutralise_molecule(clean_mol)
        return m


class Atoms(Enum):
    CARBON = Chem.Atom(6)
    OXYGEN = Chem.Atom(8)


class PeptideCapper(object):
    """Deals with splitting and capping of peptides."""

    def __init__(self):
        pass

    @staticmethod
    def _add_acetyl_to_mol(rw_mol, atom_attach):
        """Adds a C(=O)(C) to an atom in a RW molecule."""
        c_carbonyl_idx = rw_mol.AddAtom(Atoms.CARBON.value)
        o_carbonyl_idx = rw_mol.AddAtom(Atoms.OXYGEN.value)
        methyl_idx = rw_mol.AddAtom(Atoms.CARBON.value)
        rw_mol.AddBond(atom_attach, c_carbonyl_idx, Chem.BondType.SINGLE)
        rw_mol.AddBond(c_carbonyl_idx, o_carbonyl_idx, Chem.BondType.DOUBLE)
        rw_mol.AddBond(c_carbonyl_idx, methyl_idx, Chem.BondType.SINGLE)
        return rw_mol

    def _fragment_and_cap(self, mol, breakable_points):
        """
        Function for breaking and capping.
        Creates a RW mol; iterates through a set of breakable
        atom indices; for each pair of indices, splits the
        molecule in two parts and adds an acetyl capping
        to the first atom (amine NH), and a carbon to the
        second atom (carbonyl C=O). This ensures that
        both fragments are capped with acetyles.

        """
        rw_mol = Chem.RWMol(mol)
        for atom in breakable_points:
            atom_1_idx = atom[0]
            atom_2_idx = atom[1]
            # Break bond
            rw_mol.RemoveBond(atom_1_idx, atom_2_idx)
            # Add acetyl to the amine
            rw_mol = PeptideCapper._add_acetyl_to_mol(rw_mol, atom_1_idx)
            # Add carbon to the carbonyl
            c_idx = rw_mol.AddAtom(Atoms.CARBON.value)
            rw_mol.AddBond(atom_2_idx, c_idx, Chem.BondType.SINGLE)
        return Chem.MolToSmiles(rw_mol)

    def break_amide_bonds_and_cap(self, mol):
        """
        Breaks the amine bonds in a molecule and
        caps them with acetyl groups.

        """
        # SMARTS - Secondary and tertiary amide bonds
        AMIDE_SMARTS = '[NX3,NX4;H0,H1][CX3](=[OX1])'
        amide_pattern = Chem.MolFromSmarts(AMIDE_SMARTS)

        breakable_points = list()
        smiles = Chem.MolToSmiles(mol)
        if mol.HasSubstructMatch(amide_pattern):
            atom_idxs = mol.GetSubstructMatches(amide_pattern)
            for atom in atom_idxs:
                breakable_points.append((atom[0], atom[1]))
            smiles = self._fragment_and_cap(mol, breakable_points)

        smiles_list = smiles.split(".")
        log.debug(f"Obtained {len(smiles_list)} fragments")
        return smiles_list
