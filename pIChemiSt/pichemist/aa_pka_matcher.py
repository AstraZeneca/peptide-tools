import sys
from rdkit import Chem
from pka_sets_fasta import PKA_SETS_NAMES, PKA_SETS
from aa_pka_config import D_cappedAA_smarts, D_NtermfreeAA_smarts, D_CtermfreeAA_smarts


def _patt_match_rdkit(smiles, smarts):
    """RDKit SMARTS matching counter."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    pat = Chem.MolFromSmarts(smarts)
    n = 0
    for _ in mol.GetSubstructMatches(pat):
        n += 1
    return n


def patt_match(smiles, smarts):
    """Interface for SMARTS matching."""
    return _patt_match_rdkit(smiles, smarts)


def _match_d_capped_aa(smiles,
                       base_pka_fast_dict,
                       acid_pka_fast_dict):
    """TODO: """

    for smarts, aa in D_cappedAA_smarts.items():
        nhits = patt_match(smiles, smarts)
        if nhits > 0: 
            for n in PKA_SETS_NAMES:
                pka_set = PKA_SETS[n]

                # ANDREY: the iterator is not used
                # for _ in range(nhits):
                if aa in pka_set['basic'].keys():
                    basic = pka_set['basic'][aa][0]
                    base_pka_fast_dict[n].append([basic, aa])

                if aa in pka_set['acidic'].keys():
                    pka_acid = pka_set['acidic'][aa][0]
                    acid_pka_fast_dict[n].append([pka_acid, aa]) 
            return True
    return False


def get_aa_pkas(smi_list):

    # Initialise result dicts
    unknown_fragments = list()
    base_pka_fast_dict = dict()
    acid_pka_fast_dict = dict()
    diacid_pka_fast_dict = dict()
    for name in PKA_SETS_NAMES:
        base_pka_fast_dict[name] = list()
        acid_pka_fast_dict[name] = list()
        diacid_pka_fast_dict[name] = list()

    for smiles in smi_list:
        match = False

        # Middle
        match = _match_d_capped_aa(smiles,
                                   base_pka_fast_dict,
                                   acid_pka_fast_dict)
        
        for smarts, aa in D_NtermfreeAA_smarts.items():
            # N-term 
            nhits = patt_match(smiles,smarts)
            if nhits > 0: 
                for name in PKA_SETS_NAMES:
                    pka_set = PKA_SETS[name]
                    for i in range(nhits):
                        if aa in pka_set['basic'].keys(): base_pka_fast_dict[name].append( [ pka_set['basic'][aa][1] , aa ] )
                        if aa in pka_set['acidic'].keys(): acid_pka_fast_dict[name].append( [ pka_set['acidic'][aa][1] , aa ] ) 
                        base_pka_fast_dict[name].append( [ pka_set['terminus_ionizable'][aa][0], aa+'_N-term' ])  # N-terminus
                match=True
                #print(aa)
                break
 
        for smarts,aa in D_CtermfreeAA_smarts.items():
            # C-term 
            nhits = patt_match(smiles,smarts)
            if nhits > 0: 
                for name in PKA_SETS_NAMES:
                    pka_set = PKA_SETS[name]
                    for i in range(nhits):
                        if aa in pka_set['basic'].keys(): base_pka_fast_dict[name].append( [ pka_set['basic'][aa][2] , aa ] )
                        if aa in pka_set['acidic'].keys(): acid_pka_fast_dict[name].append( [ pka_set['acidic'][aa][2] , aa ] ) 
                        acid_pka_fast_dict[name].append( [ pka_set['terminus_ionizable'][aa][1], aa+'_C-term' ])  # C-terminus
                match=True
                #print(aa)
                break

        if not match:
            unknown_fragments.append(smiles)
            match=False

    return (unknown_fragments, base_pka_fast_dict, acid_pka_fast_dict, diacid_pka_fast_dict)



def get_scrambled_fasta_from_frags(smi_list):

    fasta_list=[]

    for smiles in smi_list:
        match=False

        for smarts,aa in D_cappedAA_smarts.items():
            # Middle 
            nhits = patt_match(smiles,smarts)
            if nhits > 0: 
                fasta_list += [aa]*nhits
                match=True
                break

        for smarts,aa in D_NtermfreeAA_smarts.items():
            # N-term 
            nhits = patt_match(smiles,smarts)
            if nhits > 0: 
                fasta_list += [aa]*nhits
                match=True
                break
 
        for smarts,aa in D_CtermfreeAA_smarts.items():
            # C-term 
            nhits = patt_match(smiles,smarts)
            if nhits > 0: 
                fasta_list += [aa]*nhits
                match=True
                break

        if not match:
            match=False
            fasta_list += ['X']

    return ''.join(fasta_list)


    




if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("!!!ERROR: no cmd line arg")
        sys.exit(1)

    smilesfile = sys.argv[1]

    smi_list=[]
    smif=open(smilesfile,"r")
    for l in smif.readlines():
        entries=l.split()
        smi_list.append(entries[0])
    
    print(get_pkas_for_known_AAs(smi_list))


