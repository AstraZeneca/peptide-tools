from pichemist.config import PKA_LIMITS
from pichemist.config import SKIP_SMARTS_NAMES
from pichemist.model import MODELS
from pichemist.model import PKaType
from pichemist.smarts.pka_set import SS_SMARTS_PKA_SET
from rdkit import Chem


class PKaMatcher(object):
    """
    Uses a set of SMARTS definitions and their pKa annotations
    against input SMILES to match their pKa(s).

    """

    def __init__(self):
        self.smarts_set = SS_SMARTS_PKA_SET
        self.skip_names = SKIP_SMARTS_NAMES

    def _pka_dict_from_smiles(self, smiles: str):
        """
        Iterates through a list of dicts
        where each dict represents a group with
        its SMARTS and pKa value. If the SMARTS matches
        the input SMILES argument, then its value, and
        the SMILES are added to the results.

        """
        results = {PKaType.ACIDIC.value: list(), PKaType.BASIC.value: list()}
        mol = Chem.MolFromSmiles(smiles)
        all_used_idxs_set = set()
        # Each group is described as a list of dicts
        for pka_group_list in self.smarts_set:
            for pka_dict in pka_group_list:

                # Skip if name is included in skip config
                if pka_dict["name"] in self.skip_names:
                    continue

                # Match SMARTS against the input and store pKa
                used_idx_list = list()  # Used within larger scope
                used_idx_local = list()  # Used within smaller scope
                pat = Chem.MolFromSmarts(pka_dict["smarts"])
                for match in mol.GetSubstructMatches(pat):
                    match_idx = match[pka_dict["idx"] - 1]
                    used_idx_list += match
                    match_idxs = set(match)
                    available_idxs = match_idxs.difference(all_used_idxs_set)
                    pka = pka_dict["pka"]
                    if pka_dict["type"] not in MODELS[PKaType]:
                        raise PKaMatcherException(
                            f"Unknown site type (got {pka_dict['type']})"
                        )

                    # Basic
                    if pka_dict["type"] == PKaType.BASIC.value:
                        if (
                            match_idx in available_idxs
                            and match_idx not in used_idx_local
                        ):
                            if (
                                pka > PKA_LIMITS["base_1"]
                                and pka < PKA_LIMITS["base_2"]
                            ):
                                results["base"].append((pka, smiles))
                    # Acidic
                    if pka_dict["type"] == PKaType.ACIDIC.value:
                        if (
                            match_idx in available_idxs
                            and match_idx not in used_idx_local
                        ):
                            if (
                                pka > PKA_LIMITS["acid_1"]
                                and pka < PKA_LIMITS["acid_2"]
                            ):
                                results["acid"].append((pka, smiles))
                    used_idx_local.append(match_idx)
            all_used_idxs_set = all_used_idxs_set.union(set(used_idx_list))
        return results

    def calculate_pka_from_smiles(self, smiles: str):
        """Calculates the pKa values for a SMILES."""
        results = self._pka_dict_from_smiles(smiles)
        return {k: [t[0] for t in v] for k, v in results.items()}

    def calculate_pka_from_list(self, smiles_list: list):
        """Calculates the pKa values for a list of SMILES."""
        base_pkas = list()
        acid_pkas = list()

        for smiles in smiles_list:
            results = self._pka_dict_from_smiles(smiles)
            base_pkas.extend(results[PKaType.BASIC.value])
            acid_pkas.extend(results[PKaType.ACIDIC.value])
        return (base_pkas, acid_pkas)


class PKaMatcherException(Exception):
    def __init__(self):
        pass
