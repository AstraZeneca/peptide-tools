#!/usr/bin/python
import sys, os
list_smarts_pka = []

pwd = os.path.dirname(os.path.realpath(__file__))
substructures_smarts_file = "{}/{}".format(pwd, "smarts_pKaMatcher.dat")
with open(substructures_smarts_file,'r') as f:
    for line in f.readlines():
        ln=line.split()
        if len(ln) == 0: continue
        if line[0] == '#': continue
        substr_name = ln[0]
        substr_smarts = ln[1]

        pka_l=[]
        for i in range(int(len(ln[2:])/4)):
           pl = ln[2+i*4:2+i*4+4]
           pd = {'pka': float(pl[1]), 'ind':int(pl[0]),'pka_std':float(pl[2]),'type':pl[3],'smarts':substr_smarts,'name':substr_name} 
           pka_l.append(pd)

        list_smarts_pka.append(pka_l)        



        








#       if substr_type in ['acid','base']: 
#           substr_ind = int(ln[2])
#           substr_pka = ln[3]
#           substr_pka_std = ln[4]
#           #D_dimorphite_dl_type_pka[substr_name] = {'pka': float(substr_pka),'name': substr_pka, 'ind':substr_ind,'pka_std':substr_pka_std,'type':substr_type,'smarts':substr_smarts,'name':substr_name}
#           list_smarts_pka.append( {'pka': float(substr_pka), 'ind':substr_ind,'pka_std':substr_pka_std,'type':substr_type,'smarts':substr_smarts,'name':substr_name} )
#           #Amines_primary_secondary_tertiary	[C:1]-[NX3+0:2]	1	8.159107682388349	2.5183597445318147   base
#       elif substr_type == 'diacid':
#           substr_ind1 = int(ln[2])
#           substr_pka1 = ln[3]
#           substr_pka_std1 = ln[4]
#           substr_ind2 = int(ln[5])
#           substr_pka2 = ln[6]
#           substr_pka_std2 = ln[7]

#           #D_dimorphite_dl_type_pka[substr_name] = {'pka1': float(substr_pka1),'name': substr_pka, 'ind1':substr_ind1,'pka_std1':substr_pka_std1,'type':substr_type, 'pka2': float(substr_pka2),'ind2':substr_ind2,'pka_std2':substr_pka_std2,'smarts':substr_smarts,'name':substr_name}

#           list_smarts_pka.append( {'pka1': float(substr_pka1),'ind1':substr_ind1,'pka_std1':substr_pka_std1,'type':substr_type, 'pka2': float(substr_pka2),'ind2':substr_ind2,'pka_std2':substr_pka_std2,'smarts':substr_smarts,'name':substr_name} )
#       else:
#           raise Exception('ERROR: substructure type in pKaMatcher set not known '+substr_type)
#           sys.exit(1)
            
