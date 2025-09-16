from pichemist.config import ROUNDING_DIGITS
from pichemist.molecule import MolStandardiser
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.info")


class SmartsChargeCalculator(object):
    """Calculates the nitrogen cation charges for an input using SMARTS."""

    def __init__(self):
        self.smarts = {
            # Positively charged N but not connected to any negatively
            # charged atom. To avoid azido and nitro groups being counted
            "nitrogen_cation": "[#7+;!H1;!H2;!H3;!H4!$([#7+]~[*-])]"
        }
        self.patterns = self._prepare_patterns()

    def _prepare_patterns(self):
        """Converts SMARTS strings into objects."""
        return {name: Chem.MolFromSmarts(s) for name, s in self.smarts.items()}

    def _get_mol_from_smiles(self, smiles):
        """Get mol object from SMILES."""
        return Chem.MolFromSmiles(smiles)

    def calculate_net_qs_from_list(self, smiles_list):
        """Produces a list of charges and matching SMILES."""
        net_qs = list()
        for smiles in smiles_list:
            mol = self._get_mol_from_smiles(smiles)
            mol = MolStandardiser().standardise_molecule(mol)
            # sic - Appends charge 1 and the SMILES for each match
            for _ in self._get_net_qs_matches_from_mol(mol):
                net_qs.append((1, smiles))
        return net_qs

    def _get_net_qs_matches_from_mol(self, mol):
        """
        Matches a nitrogen cation pattern against a
        molecule and returns a list of matches.

        """
        pattern = self.patterns["nitrogen_cation"]
        matches = mol.GetSubstructMatches(pattern)
        return [y[0] for y in matches]

    def calculate_net_qs_from_smiles(self, smiles):
        """Returns the number of charges for a given SMILES."""
        mol = self._get_mol_from_smiles(smiles)
        mol = MolStandardiser().standardise_molecule(mol)
        # sic - The charges corresponds to the number of matches
        number_of_charges = len(self._get_net_qs_matches_from_mol(mol))
        return number_of_charges


class PKaChargeCalculator(object):
    """Calculate the molecule charge given its pKas."""

    def __init__(self):
        pass

    def _calculate_basic_charge(self, pH, pKa):
        """Calculate charge contribution by basic pKas."""
        return 1 / (1 + 10 ** (pH - pKa))

    def _calculate_acidic_charge(self, pH, pKa):
        """Calculate charge contribution by acidic pKas."""
        return -1 / (1 + 10 ** (pKa - pH))

    def _calculate_diacidic_charge(self, pH, pKa1, pKa2):
        """Calculate charge contribution by diacidic pKas."""
        Ka1 = 10 ** (-pKa1)
        Ka2 = 10 ** (-pKa2)
        H = 10 ** (-pH)
        f1 = (H * Ka1) / (H**2 + H * Ka1 + Ka1 * Ka2)  # fraction of [AH-]
        f2 = f1 * Ka2 / H  # fraction of [A2-]
        return -2 * f2 + (-1) * f1  # average charge of phosphate

    def calculate_charge(self, base_pkas, acid_pkas, pH, constant_q=0):
        """Calculate the molecule charge from all contributions."""
        charge = constant_q
        for pka in base_pkas:
            charge += self._calculate_basic_charge(pH, pka)
        for pka in acid_pkas:
            charge += self._calculate_acidic_charge(pH, pka)
        return round(charge, ROUNDING_DIGITS)

    def calculate_constant_charge(self, net_qs):
        """Calculates the constant charge from the net charges."""
        constant_q = 0.0
        if len(net_qs) > 0:
            constant_q = float(sum(net_qs)) / float(len(net_qs))
        return constant_q
