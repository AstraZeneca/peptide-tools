import os
import json

from pichemist.config import SMARTS_PKA_SET_FILEPATH

def create_logic_set_from_standardised_json(filepath):
    with open(filepath) as f:
        data = json.load(f)
    
    logic_set = list()
    for smarts, v in data["smarts"].items():
        res_data = list()
        name = v["name"]
        data = v["data"]
        for d in data:
            res_data.append({
                "pka": d["pka"],
                "ind": d["idx"],
                "pka_std": d["pka_std"],
                "type": d["type"],
                "smarts": smarts,
                "name": name
            })
        logic_set.append(res_data)       
    return logic_set

# Read in the data
SMARTS_PKA_SET = create_logic_set_from_standardised_json(
    SMARTS_PKA_SET_FILEPATH)

if __name__ == "__main__":
    pass
    # """
    # Conversion of DAT semi-structured data into JSON
    # """
    # import json
    # data = read_smarts_pKaMatcher()
    # print(data)
    # result = dict()
    # for e in data:
    #     res_data = list()
    #     for i in e:
    #         smarts = i["smarts"]
    #         name = i["name"]
    #         res_data.append({
    #             "pka": i["pka"],
    #             "pka_std": i["pka_std"],
    #             "idx": i["ind"],
    #             "type": i["type"]
    #         })

    #         result[smarts] = {
    #             "name": name,
    #             "data": res_data
    #         }
    
    # with open("smarts_converted.json", "w") as f:
    #     json.dump({"smarts": result}, f, indent=4)

    # create_logic_set_from_standardised_json
    # print(result)
    # logic_set = list()

    # for smarts, v in result.items():
    #     res_data = list()
    #     name = v["name"]
    #     data = v["data"]
    #     for d in data:
    #         res_data.append({
    #             "pka": d["pka"],
    #             "ind": d["idx"],
    #             "pka_std": d["pka_std"],
    #             "type": d["type"],
    #             "smarts": smarts,
    #             "name": name
    #         })
    #     logic_set.append(res_data)
    
    # # assert logic_set == data
    # # print(logic_set)
    # with open("../expected.json") as f:
    #     exp = json.load(f)
    
    # assert exp == logic_set