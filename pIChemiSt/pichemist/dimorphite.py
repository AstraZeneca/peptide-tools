#def calc_pkas_dimorphite_dl(unknown_fragments):
#   from rdkit import Chem
#   from dimorphite_dl import Protonate
#   from smarts_matcher_nonnaturals_dimorphite import D_dimorphite_dl_type_pka

#   skip_site_names = ['TATA','*Amide']

#   base_pkas=[]
#   acid_pkas=[]
#   diacid_pkas=[]

#   for smiles in unknown_fragments:
#       protonation_sites = Protonate({'smiles':smiles}).get_protonation_sites()

#       for sites in protonation_sites:
#           #(3, 'BOTH', '*Amide') 
#           site_name = sites[2]

#           if site_name in skip_site_names: continue

#           if site_name in D_dimorphite_dl_type_pka.keys():
#               site_data = D_dimorphite_dl_type_pka[site_name]
#               if site_data['type']=='base':
#                       base_pkas.append((site_data['pka'],smiles))
#               if site_data['type']=='acid':
#                       acid_pkas.append((site_data['pka'],smiles))
#               if site_data['type']=='diacid':
#                       diacid_pkas.append(((site_data['pka1'],site_data['pka2']),smiles))
#               
#           else:
#               print('Error:  not known dimorphite site type ' + site_type)
#               sys.exit(1)
#   return (base_pkas,acid_pkas,diacid_pkas)
