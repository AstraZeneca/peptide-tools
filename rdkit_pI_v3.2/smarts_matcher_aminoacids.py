#!/usr/bin/python
import sys
from pka_sets_fasta import *
from dimorphite_dl_site_substructures_smarts import data_txt

# Get Dimorphite DL pKa sets (https://durrantlab.pitt.edu/dimorphite-dl/)
D_dimorphite_dl_type_pka = {}
#with open('site_substructures.smarts','r') as f:
for line in data_txt.splitlines():
        ln=line.split()
        if len(ln) == 0: continue
        if line[0] == '#': continue
        substr_name = ln[0]
        substr_smarts = ln[1]
        substr_type = ln[-1]
        if substr_type in ['acid','base']: 
            substr_ind = ln[2]
            substr_pka = ln[3]
            substr_pka_std = ln[4]
            D_dimorphite_dl_type_pka[substr_name] = {'pka': float(substr_pka),'name': substr_pka, 'ind':substr_ind,'pka_std':substr_pka_std,'type':substr_type,'smarts':substr_smarts,'name':substr_name}
            #Amines_primary_secondary_tertiary	[C:1]-[NX3+0:2]	1	8.159107682388349	2.5183597445318147   base
        elif substr_type == 'diacid':
            substr_ind1 = ln[2]
            substr_pka1 = ln[3]
            substr_pka_std1 = ln[4]
            substr_ind2 = ln[2]
            substr_pka2 = ln[3]
            substr_pka_std2 = ln[4]
            D_dimorphite_dl_type_pka[substr_name] = {'pka1': float(substr_pka1),'name': substr_pka, 'ind1':substr_ind1,'pka_std1':substr_pka_std1,'type':substr_type, 'pka2': float(substr_pka2),'ind2':substr_ind2,'pka_std2':substr_pka_std2,'smarts':substr_smarts,'name':substr_name}
        else:
            print('ERROR: substructure type in Dimorphine DL set not known '+substr_type)
            sys.exit(1)
            


D_cappedAA_smarts = {
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]OC(=O)C([H])([H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[O-]C(=O)C([H])([H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])O[H])N([H])C(=O)C([H])([H])[H]':   'S',
'[H]C([H])([H])C(=O)C([H])(C([H])(C([H])([H])[H])O[H])N([H])C(=O)C([H])([H])[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'N',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])S[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'A',
'[H]C([H])([H])C(=O)C([H])(C([H])(C([H])([H])[H])C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)C([H])(C([S-])([H])[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])([H])N([H])C(=O)C([H])([H])[H]':   'G',
'[H]C([H])([H])C(=O)C1([H])C([H])([H])C([H])([H])C([H])([H])N1C(=O)C([H])([H])[H]':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'I',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]OC(=O)C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[O-]C(=O)C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[H]N([H])C(=O)C([H])(C([H])([H])O[H])N([H])C(=O)C([H])([H])[H]':   'S',
'[H]C([H])([H])C([H])(O[H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'N',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]N([H])C(=O)C([H])(C([H])([H])S[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]N([H])C(=O)C([H])(C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'A',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]N([H])C(=O)C([H])(C([S-])([H])[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]N([H])C(=O)C([H])([H])N([H])C(=O)C([H])([H])[H]':   'G',
'[H]C([H])([H])C(=O)N1C([H])([H])C([H])([H])C([H])([H])C1([H])C(=O)N([H])[H]':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'I'
}


D_NtermfreeAA_smarts = {
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])(C([H])([H])[H])O[H]':   'T',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C(=O)N([H])[H]':   'N',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])[H])N([H])[H]':   'A',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'V',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)C([H])(C([S-])([H])[H])N([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])([H])N([H])[H]':   'G',
'[H]C([H])([H])C(=O)C1([H])C([H])([H])C([H])([H])C([H])([H])N1[H]':   'P',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])(C([H])([H])[H])C([H])([H])C([H])([H])[H]':   'I',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])(C([H])([H])[H])O[H]':   'T',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C(=O)N([H])[H]':   'N',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])[H])[N+]([H])([H])[H]':   'A',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'V',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)C([H])(C([S-])([H])[H])[N+]([H])([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])([H])[N+]([H])([H])[H]':   'G',
'[H]C([H])([H])C(=O)C1([H])C([H])([H])C([H])([H])C([H])([H])[N+]1([H])[H]':   'P',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])(C([H])([H])[H])C([H])([H])C([H])([H])[H]':   'I'
}

D_CtermfreeAA_smarts = {
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])(C([H])([H])[H])O[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'N',
'[H]N([H])C(=O)C([H])([H])C([H])([H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'Q',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])(C([H])([H])[H])C(=O)O[H]':   'A',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)N([H])C([H])(C([S-])([H])[H])C(=O)O[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])([H])C(=O)O[H]':   'G',
'[H]C([H])([H])C(=O)N1C([H])([H])C([H])([H])C([H])([H])C1([H])C(=O)O[H]':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'I',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])(C([H])([H])[H])O[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'N',
'[H]N([H])C(=O)C([H])([H])C([H])([H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'Q',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])[H]':   'A',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([S-])([H])[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])([H])C([O-])=O':   'G',
'[H]C([H])([H])C(=O)N1C([H])([H])C([H])([H])C([H])([H])C1([H])C([O-])=O':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'I'
}


def patmatch_OE(smiles,smarts):
    from openeye.oechem import OEmol, OEParseSmiles, OEAddExplicitHydrogens, OESubSearch
    mol = OEMol()
    mol.Clear()
    OEParseSmiles(mol, smiles)
    OEAddExplicitHydrogens(mol)
    pat = OESubSearch()
    pat.Init(smarts)
    n = 0
    for match in pat.Match(mol,True):
        n += 1
    return n

def patmatch_RDKit(smiles,smarts):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    pat = Chem.MolFromSmarts(smarts)
    #m.HasSubstructMatch(patt)
    n = 0
    for match in mol.GetSubstructMatches(pat):
        n += 1
    return n

def patmatch(smiles,smarts):
    # smarts pattern matching
    # manual switch between RDkit or OpenEye. Just for historical reason.
    lRDKit=True
    lOE=False
    if lRDKit:   n = patmatch_RDKit(smiles,smarts)
    elif lOE:    n = patmatch_OE(smiles,smarts)
    return n



def get_pkas_for_known_AAs(smilist):

    unknown_fragments = []

    base_pka_fasta_list = {}
    acid_pka_fasta_list = {}
    diacid_pka_fasta_list = {}
    for pKaset in pKa_sets_to_use:
        base_pka_fasta_list[pKaset] = []
        acid_pka_fasta_list[pKaset] = []
        diacid_pka_fasta_list[pKaset] = []

    for smiles in smilist:
        lmatch=False

        for smarts,AA in D_cappedAA_smarts.items():
            # Middle 
            nhits = patmatch(smiles,smarts)
            if nhits > 0: 
                for pKaset in pKa_sets_to_use:
                    pKa_set = pKa_sets[pKaset]
                    for i in range(nhits):
                        if AA in pKa_set['pKa_basic'].keys(): base_pka_fasta_list[pKaset].append( [ pKa_set['pKa_basic'][AA][0] , AA ] )
                        if AA in pKa_set['pKa_acidic'].keys(): acid_pka_fasta_list[pKaset].append( [ pKa_set['pKa_acidic'][AA][0] , AA ] ) 
                lmatch=True
                #print(AA)
                break

        for smarts,AA in D_NtermfreeAA_smarts.items():
            # N-term 
            nhits = patmatch(smiles,smarts)
            if nhits > 0: 
                for pKaset in pKa_sets_to_use:
                    pKa_set = pKa_sets[pKaset]
                    for i in range(nhits):
                        if AA in pKa_set['pKa_basic'].keys(): base_pka_fasta_list[pKaset].append( [ pKa_set['pKa_basic'][AA][1] , AA ] )
                        if AA in pKa_set['pKa_acidic'].keys(): acid_pka_fasta_list[pKaset].append( [ pKa_set['pKa_acidic'][AA][1] , AA ] ) 
                        base_pka_fasta_list[pKaset].append( [ pKa_set['pKa_TerminusIonizableGroup'][AA][0], AA+'_N-term' ])  # N-terminus
                lmatch=True
                #print(AA)
                break
 
        for smarts,AA in D_CtermfreeAA_smarts.items():
            # C-term 
            nhits = patmatch(smiles,smarts)
            if nhits > 0: 
                for pKaset in pKa_sets_to_use:
                    pKa_set = pKa_sets[pKaset]
                    for i in range(nhits):
                        if AA in pKa_set['pKa_basic'].keys(): base_pka_fasta_list[pKaset].append( [ pKa_set['pKa_basic'][AA][2] , AA ] )
                        if AA in pKa_set['pKa_acidic'].keys(): acid_pka_fasta_list[pKaset].append( [ pKa_set['pKa_acidic'][AA][2] , AA ] ) 
                        acid_pka_fasta_list[pKaset].append( [ pKa_set['pKa_TerminusIonizableGroup'][AA][1], AA+'_C-term' ])  # C-terminus
                lmatch=True
                #print(AA)
                break

        if not lmatch:
            unknown_fragments.append(smiles)
            lmatch=False
            #print(smiles)
        #else:
            #base_pka_fasta_list.append([base_pka_fasta,AA])
            #acid_pka_fasta_list.append([acid_pka_fasta,AA])

    return (unknown_fragments, base_pka_fasta_list, acid_pka_fasta_list, diacid_pka_fasta_list)



def get_scrambled_fasta_from_frags(smilist):

    fasta_list=[]

    for smiles in smilist:
        lmatch=False

        for smarts,AA in D_cappedAA_smarts.items():
            # Middle 
            nhits = patmatch(smiles,smarts)
            if nhits > 0: 
                fasta_list += [AA]*nhits
                lmatch=True
                break

        for smarts,AA in D_NtermfreeAA_smarts.items():
            # N-term 
            nhits = patmatch(smiles,smarts)
            if nhits > 0: 
                fasta_list += [AA]*nhits
                lmatch=True
                break
 
        for smarts,AA in D_CtermfreeAA_smarts.items():
            # C-term 
            nhits = patmatch(smiles,smarts)
            if nhits > 0: 
                fasta_list += [AA]*nhits
                lmatch=True
                break

        if not lmatch:
            lmatch=False
            fasta_list += ['X']

    return ''.join(fasta_list)


    




if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("!!!ERROR: no cmd line arg")
        sys.exit(1)

    smilesfile = sys.argv[1]

    smilist=[]
    smif=open(smilesfile,"r")
    for l in smif.readlines():
        entries=l.split()
        smilist.append(entries[0])
    
    print(get_pkas_for_known_AAs(smilist))


