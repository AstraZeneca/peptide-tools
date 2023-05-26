from rdkit import Chem
from rdkit import RDLogger
from pichemist.molecule import MolStandardiser

RDLogger.DisableLog("rdApp.info")


class ChargeCalculator(object):
    """Calculates the nitrogen cation charges for an input."""

    def __init__(self):
        self.smarts = {
            # Positively charged N but not connected to any negatively
            # charged atom. To avoid azido and nitro groups being counted
            "nitrogen_cation": "[#7+;!H1;!H2;!H3;!H4!$([#7+]~[*-])]"
        }
        self.patterns = self._prepare_patterns()

    def _prepare_patterns(self):
        """Converts SMARTS strings into objects."""
        return {name: Chem.MolFromSmarts(s)
                for name, s in self.smarts.items()}

    def calculate_net_qs_from_list(self, smiles_list):
        """Produces a list of charges and matching SMILES."""
        net_qs = list()
        for smiles in smiles_list:
            # sic - Appends one charge and the SMILES for each match
            for _ in self._get_net_qs_matches_from_smiles(smiles):
                net_qs.append((1, smiles))
        return net_qs

    def _get_net_qs_matches_from_smiles(self, smiles):
        """
        Matches a nitrogen cation pattern against a
        molecule and returns a list of matches.

        """
        mol = Chem.MolFromSmiles(smiles)
        # ANDREY: Do we need to standardise again???
        mol = MolStandardiser().standardise_molecule(mol)
        pattern = self.patterns["nitrogen_cation"]
        matches = mol.GetSubstructMatches(pattern)
        return [y[0] for y in matches]

    def calculate_net_qs_from_smiles(self, smiles):
        """Returns the number of charges for a given SMILES."""
        # sic - The charges corresponds to the number of matches
        return len(self._get_net_qs_matches_from_smiles(smiles))
