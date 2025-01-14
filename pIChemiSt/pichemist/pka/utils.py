from pichemist.fasta.matcher import FastaPKaMatcher


def _merge_pka_list(list_of_lists):
    """Merges a list of lists of pKa values."""
    PKA_VALUE_IDX = 0
    merged_pkas = list()
    for e in list_of_lists:
        for tup in e:
            merged_pkas.append(tup[PKA_VALUE_IDX])
    return merged_pkas


def _unpack_pka_list(pka_list):
    """Unpacks a list of pKa lists."""
    return pka_list[0], pka_list[1]


def merge_pkas_lists(base_pkas_list, acid_pkas_list):
    # Initialise the merged pKa lists
    base_pkas = dict()
    acid_pkas = dict()

    # Unpack the pKa values
    base_pkas_fasta, base_pkas_calc = _unpack_pka_list(base_pkas_list)
    acid_pkas_fasta, acid_pkas_calc = _unpack_pka_list(acid_pkas_list)

    # Merge FASTA and calculated values
    pka_sets_names = FastaPKaMatcher().get_pka_sets_names()
    for pka_set in pka_sets_names:
        base_pkas[pka_set] = _merge_pka_list([base_pkas_fasta[pka_set], base_pkas_calc])
        acid_pkas[pka_set] = _merge_pka_list([acid_pkas_fasta[pka_set], acid_pkas_calc])
    return base_pkas, acid_pkas
