from pichemist.config import AA_TYPE_KEYS
from pichemist.config import PKA_SETS_NAMES
from pichemist.fasta.helpers import pattern_match
from pichemist.fasta.pka_sets import FASTA_PKA_SETS
from pichemist.fasta.smarts import AA_SMARTS_SET


class FastaPKaMatcher(object):
    def __init__(self):
        self.pka_sets = PKA_SETS_NAMES
        self.aa_smarts_set = AA_SMARTS_SET

    def _initialise_pka_sets(self, all_pka_sets=FASTA_PKA_SETS,
                             pka_names=PKA_SETS_NAMES):
        """Generates a filtered set of pKa sets."""
        pka_sets = dict()
        for n in pka_names:
            pka_sets[n] = all_pka_sets[n]
        return pka_sets

    def get_pka_sets_names(self):
        """Getter for the names of the sets."""
        return self.pka_sets

    def _initialise_pka_dicts(self, pka_sets):
        """Initialises empty pKa dictionaries."""
        base_pka_dict = dict()
        acid_pka_dict = dict()
        diacid_pka_dict = dict()

        for n in pka_sets:
            base_pka_dict[n] = list()
            acid_pka_dict[n] = list()
            diacid_pka_dict[n] = list()
        return base_pka_dict, acid_pka_dict, diacid_pka_dict

    def _get_pka_from_set(self, aa, pka_set, pka_idx):
        """Returns a pKa value from a set."""
        return pka_set[aa][pka_idx]

    def _add_pka_to_dict(self, aa, name, pka, pka_dict):
        """Assigns a (pKa, AA) tuple to a dictionary key."""
        pka_dict[name].append((pka, aa))
        return pka_dict

    def _add_pka_to_dict_if_present(self, aa, name, pka_dict,
                                    pka_set, pka_idx):
        """Adds a pKa value to a dict if present in a query set."""
        if aa in pka_set.keys():
            pka = self._get_pka_from_set(aa, pka_set, pka_idx)
            pka_dict = self._add_pka_to_dict(aa, name, pka, pka_dict)

    def _add_pka_to_acidic_and_basic_dicts(self, aa, pka_set_name,
                                           base_pka_dict,
                                           acid_pka_dict,
                                           pka_set,
                                           aa_idx):
        """
        Wraps the logic to add a pKa value to
        acidic and basic dicts.

        """
        self._add_pka_to_dict_if_present(aa, pka_set_name,
                                         base_pka_dict,
                                         pka_set["basic"],
                                         aa_idx)
        self._add_pka_to_dict_if_present(aa, pka_set_name,
                                         acid_pka_dict,
                                         pka_set["acidic"],
                                         aa_idx)

    def _build_name(self, aa, aa_type):
        """Builds a name with AA and its type."""
        return f"{aa}_{aa_type}"

    def _add_terminus_ionizable_to_dict(self, aa, pka_set_name,
                                        pka_dict,
                                        pka_set,
                                        aa_type,
                                        aa_idx):
        """Adds a terminus pka value to a dict."""
        self._add_pka_to_dict(self._build_name(aa, aa_type),
                              pka_set_name,
                              self._get_pka_from_set(
            aa,
            pka_set['terminus_ionizable'],
            aa_idx),
                              pka_dict)

    def _get_aa_pkas(self, smiles,
                     unknown_fragments,
                     base_pka_dict,
                     acid_pka_dict,
                     diacid_pka_dict,  # sic - not used
                     pka_sets):
        """
        For a given SMILES, it matches its pKa
        value against a set of pKa sets.
        NOTE: Tried to break this function down
        further but integration tests break.

        """
        match = False

        # Middle
        CAPPED_IDX = 0
        CAPPED_KEY = AA_TYPE_KEYS["capped"]
        CAPPED_SMARTS = self.aa_smarts_set[CAPPED_KEY]
        for smarts, aa in CAPPED_SMARTS.items():
            nhits = pattern_match(smiles, smarts)
            if nhits > 0:
                for n, pka_set in pka_sets.items():
                    # Add pKa values for each match
                    for _ in range(nhits):
                        self._add_pka_to_acidic_and_basic_dicts(aa, n,
                                                                base_pka_dict,
                                                                acid_pka_dict,
                                                                pka_set,
                                                                CAPPED_IDX)
                match = True
                break

        # N-term
        NTERM_FREE_IDX = 1
        NTERM_ION_IDX = 0
        NTERM_KEY = AA_TYPE_KEYS["nterm"]
        NTERM_SMARTS = self.aa_smarts_set[NTERM_KEY]
        for smarts, aa in NTERM_SMARTS.items():
            nhits = pattern_match(smiles, smarts)
            if nhits > 0:
                for n, pka_set in pka_sets.items():
                    # Add pKa values for each match
                    for _ in range(nhits):
                        self._add_pka_to_acidic_and_basic_dicts(aa, n,
                                                                base_pka_dict,
                                                                acid_pka_dict,
                                                                pka_set,
                                                                NTERM_FREE_IDX)
                        self._add_terminus_ionizable_to_dict(aa, n,
                                                             base_pka_dict,
                                                             pka_set,
                                                             "N-term",
                                                             NTERM_ION_IDX)
                match = True
                break

        # C-term
        CTERM_FREE_IDX = 2
        CTERM_ION_IDX = 1
        CTERM_KEY = AA_TYPE_KEYS["cterm"]
        CTERM_SMARTS = self.aa_smarts_set[CTERM_KEY]
        for smarts, aa in CTERM_SMARTS.items():
            nhits = pattern_match(smiles, smarts)
            if nhits > 0:
                for n, pka_set in pka_sets.items():
                    # Add pKa values for each match
                    for _ in range(nhits):
                        self._add_pka_to_acidic_and_basic_dicts(aa, n,
                                                                base_pka_dict,
                                                                acid_pka_dict,
                                                                pka_set,
                                                                CTERM_FREE_IDX)
                        self._add_terminus_ionizable_to_dict(aa, n,
                                                             acid_pka_dict,
                                                             pka_set,
                                                             "C-term",
                                                             CTERM_ION_IDX)
                match = True
                break

        # Append as unknown if not matched
        if not match:
            unknown_fragments.append(smiles)
        return unknown_fragments, base_pka_dict, acid_pka_dict, diacid_pka_dict

    def get_aa_pkas_from_list(self, smiles_list):
        """
        Matches a list of SMILES against
        the pKa values in a set of pKa sets.

        """
        # Initialise pKa sets
        pka_sets = self._initialise_pka_sets()

        # Initialise results
        unknown_fragments = list()
        base_pka_dict, acid_pka_dict, diacid_pka_dict = \
            self._initialise_pka_dicts(PKA_SETS_NAMES)

        for s in smiles_list:
            unknown_fragments, base_pka_dict, \
                acid_pka_dict, diacid_pka_dict = self._get_aa_pkas(
                    s, unknown_fragments, base_pka_dict,
                    acid_pka_dict, diacid_pka_dict, pka_sets)
        return unknown_fragments, base_pka_dict,\
            acid_pka_dict, diacid_pka_dict
