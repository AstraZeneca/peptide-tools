from pichemist.config import AA_TYPE_KEYS
from pichemist.fasta.helpers import pattern_match
from pichemist.fasta.smarts import AA_SMARTS_SET


class FastaScrambler(object):
    def __init__(self):
        self.aa_smarts_set = AA_SMARTS_SET

    def get_scrambled_fasta_from_list(self, smiles_list):
        """
        Retrieves amino acid letters (FASTA) for a list of
        SMILES and concatenates them into a string.

        """
        aa_list = list()
        for smiles in smiles_list:
            match = False

            # Middle
            CAPPED_KEY = AA_TYPE_KEYS["capped"]
            CAPPED_SMARTS = self.aa_smarts_set[CAPPED_KEY]
            for smarts, aa in CAPPED_SMARTS.items():
                nhits = pattern_match(smiles, smarts)
                if nhits > 0:
                    aa_list += [aa]*nhits
                    match = True
                    break

            # N-term
            NTERM_KEY = AA_TYPE_KEYS["nterm"]
            NTERM_SMARTS = self.aa_smarts_set[NTERM_KEY]
            for smarts, aa in NTERM_SMARTS.items():
                nhits = pattern_match(smiles, smarts)
                if nhits > 0:
                    aa_list += [aa]*nhits
                    match = True
                    break

            # C-term
            CTERM_KEY = AA_TYPE_KEYS["cterm"]
            CTERM_SMARTS = self.aa_smarts_set[CTERM_KEY]
            for smarts, aa in CTERM_SMARTS.items():
                nhits = pattern_match(smiles, smarts)
                if nhits > 0:
                    aa_list += [aa]*nhits
                    match = True
                    break

            if not match:
                match = False
                aa_list += ["X"]
        return "".join(aa_list)
