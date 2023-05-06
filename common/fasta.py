def get_scrambled_fasta_from_frags(smi_list):

    fasta_list = []

    for smiles in smi_list:
        match = False

        for smarts,aa in d_capped_aa_smarts.items():
            # Middle 
            nhits = pattern_match(smiles,smarts)
            if nhits > 0: 
                fasta_list += [aa]*nhits
                match=True
                break

        for smarts,aa in d_nterm_free_aa_smarts.items():
            # N-term 
            nhits = pattern_match(smiles,smarts)
            if nhits > 0: 
                fasta_list += [aa]*nhits
                match=True
                break
 
        for smarts,aa in d_cterm_free_aa_smarts.items():
            # C-term 
            nhits = pattern_match(smiles,smarts)
            if nhits > 0: 
                fasta_list += [aa]*nhits
                match=True
                break

        if not match:
            match=False
            fasta_list += ['X']

    return ''.join(fasta_list)
