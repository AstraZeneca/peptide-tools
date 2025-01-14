import base64
from enum import Enum
from io import BytesIO

from pichemist.utils import get_logger
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.MolStandardize import rdMolStandardize

log = get_logger(__name__)


class Atoms(Enum):
    CARBON = Chem.Atom(6)
    OXYGEN = Chem.Atom(8)


class MolStandardiser(object):
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
        ION_SMARTS = (
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"  # noqa: E501
        )

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


class PeptideCutter(object):
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

    def _fragment_and_cap(self, mol, amide_atoms):
        """
        Function for breaking and capping.
        Creates a RW mol; iterates through a set of breakable
        atom indices; for each pair of indices, splits the
        molecule in two parts and adds an acetyl capping
        to the first atom (nitrogen of the amine NH), and a carbon
        to the second atom (carbon of the carbonyl C=O). This
        ensures that both fragments are capped with acetyles.

        """
        rw_mol = Chem.RWMol(mol)
        for atoms in amide_atoms:
            amine_idx = atoms[0]
            carbonyl_idx = atoms[1]
            # Break bond
            rw_mol.RemoveBond(amine_idx, carbonyl_idx)
            # Add acetyl to the amine
            rw_mol = PeptideCutter._add_acetyl_to_mol(rw_mol, amine_idx)
            # Add carbon to the carbonyl
            c_idx = rw_mol.AddAtom(Atoms.CARBON.value)
            rw_mol.AddBond(carbonyl_idx, c_idx, Chem.BondType.SINGLE)
        return Chem.MolToSmiles(rw_mol)

    def break_amide_bonds_and_cap(self, mol):
        """
        Breaks the amine bonds in a molecule and
        caps them with acetyl groups.

        """
        # Secondary and tertiary amide bonds
        AMIDE_SMARTS = "[NX3,NX4;H0,H1][CX3](=[OX1])"
        amide_pattern = Chem.MolFromSmarts(AMIDE_SMARTS)

        amide_atoms = list()
        smiles = Chem.MolToSmiles(mol)
        if mol.HasSubstructMatch(amide_pattern):
            amide_matches = mol.GetSubstructMatches(amide_pattern)
            for atoms in amide_matches:
                amine_idx = atoms[0]
                carbonyl_idx = atoms[1]
                amide_atoms.append((amine_idx, carbonyl_idx))
            smiles = self._fragment_and_cap(mol, amide_atoms)

        smiles_list = smiles.split(".")
        log.debug(f"Obtained {len(smiles_list)} fragments")
        return smiles_list


def smiles_to_image(smiles, b64encode=True):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, kekulize=True)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    img_byte_data = buffered.getvalue()
    if b64encode:
        base64_image = base64.b64encode(img_byte_data).decode("utf-8")
    return base64_image
