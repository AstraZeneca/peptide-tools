from rdkit import Chem
from pichemist.config import SKIP_SMARTS_NAMES
from pichemist.config import PKA_LIMITS
from pichemist.config import PKA_TYPES
from pichemist.smarts.pka_set import SS_SMARTS_PKA_SET
from pichemist.model import PKaType


class PKaMatcher(object):
    """TODO:..."""
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
        results = {PKaType.ACIDIC.value: list(),
                   PKaType.BASIC.value: list()}
        mol = Chem.MolFromSmiles(smiles)
        all_used_idxs_set = set()
        # Each group is described as a list of dicts
        for pka_group_list in self.smarts_set:
            for pka_dict in pka_group_list:
                # Skip if name is included in skip config
                if pka_dict['name'] in self.skip_names:
                    continue
                used_idx_list = list()   # Used within larger scope
                used_idx_local = list()  # Used within smaller scope
                pat = Chem.MolFromSmarts(pka_dict['smarts'])
                for match in mol.GetSubstructMatches(pat):      # ANDREY: There is no break here. Are we expecting multi matching?
                    match_idx = match[pka_dict['idx']-1]
                    used_idx_list += match
                    match_idxs = set(match)
                    available_idxs = match_idxs.difference(all_used_idxs_set)
                    pka = pka_dict['pka']
                    if pka_dict['type'] not in PKA_TYPES:
                        raise PKaMatcherException(
                            f"Unknown site type (got {pka_dict['type']})")

                    # Basic and acid pKa matching
                    # TODO: See below
                    if pka_dict['type'] == PKaType.BASIC.value:
                        if match_idx in available_idxs \
                                and match_idx not in used_idx_local:
                            if pka > PKA_LIMITS["base_1"] \
                                    and pka < PKA_LIMITS["base_2"]:
                                results["base"].append((pka, smiles))
                    if pka_dict['type'] == PKaType.ACIDIC.value:
                        if match_idx in available_idxs \
                                and match_idx not in used_idx_local:
                            if pka > PKA_LIMITS["acid_1"] \
                                    and pka < PKA_LIMITS["acid_2"]:
                                results["acid"].append((pka, smiles))
                    used_idx_local.append(match_idx)
            all_used_idxs_set = all_used_idxs_set.union(set(used_idx_list))
        return results

    def calculate_pka_from_list(self, smiles_list: list):
        base_pkas = list()
        acid_pkas = list()
        diacid_pkas = list()    # ANDREY: not used

        for smiles in smiles_list:
            results = self._pka_dict_from_smiles(smiles)
            base_pkas.extend(results[PKaType.BASIC.value])
            acid_pkas.extend(results[PKaType.ACIDIC.value])
        return (base_pkas, acid_pkas, diacid_pkas)


class PKaMatcherException(Exception):
    def __init__(self):
        pass


if __name__ == "__main__":
    matcher = PKaMatcher()
    unknown_fragments = ['CC(=O)N[C@@H](CCCN)C(C)=O',
                         'CC(=O)N[C@@](C)(CC(=O)O)C(C)=O']
    res = matcher.calculate_pka_from_list(unknown_fragments)
    expected = ([(10.4, 'CC(=O)N[C@@H](CCCN)C(C)=O')],
                [(3.46, 'CC(=O)N[C@@](C)(CC(=O)O)C(C)=O')], [])
    assert res == expected, f"got {res}"

    fragment = 'CC(=O)N[C@@H](CCCN)C(C)=O'
    res = matcher._pka_dict_from_smiles(fragment)
    # ANDREY: How possibly the same fragment can
    # return both acid and basic or multiple values?
    # This should be slimmed down into {"type": "acid/base", "value": 10.4}
    expected = {'acid': [], 'base': [(10.4, 'CC(=O)N[C@@H](CCCN)C(C)=O')]}
    assert res == expected, f"got {res}"
