#!/usr/bin/python
import sys
#from pka_sets_fasta import *
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
            


