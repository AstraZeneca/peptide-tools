#!/usr/bin/python
# module to read "smarts_pKaMatcher.dat" into list_smarts_pka
import sys, os

def read_smarts_pKaMatcher():
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

            # Loop to read ionization data. Can have several ionizations for each smarts.
            pka_l=[]
            for i in range(int(len(ln[2:])/4)):
                
               pl = ln[2+i*4:2+i*4+4] # offset of 2 and 4 entries for each pKa value: ato index, pka, pka_std, type
               pd = {'pka': float(pl[1]), 'ind':int(pl[0]),'pka_std':float(pl[2]),'type':pl[3],'smarts':substr_smarts,'name':substr_name} 
               pka_l.append(pd)

            list_smarts_pka.append(pka_l)        
    return list_smarts_pka

# Read in the data
list_smarts_pka = read_smarts_pKaMatcher()

if __name__ == "__main__":
    """
    Conversion of DAT semi-structured data into JSON
    """
    import json
    data = read_smarts_pKaMatcher()
    result = {}
    for k in data:
        res_data = []
        for i in k:
            smarts = i["smarts"]
            name = i["name"]
            res_data.append({
                "pka": i["pka"],
                "pka_std": i["pka_std"],
                "idx": i["ind"],
                "type": i["type"]
            })

            result[smarts] = {
                "name": name,
                "data": res_data
            }
    with open("smarts_converted.json", "w") as f:
        json.dump({"smarts": result}, f)