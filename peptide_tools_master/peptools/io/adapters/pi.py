def _append_pichemist_props_to_mol(mol, pichemist_res):
    mol.SetProp("pI mean", "%.2f" % pichemist_res["pI"]["pI mean"])
    mol.SetProp("pI std", "%.2f" % pichemist_res["pI"]["std"])
    mol.SetProp(
        "pI interval",
        " - ".join(["%.2f" % x for x in pichemist_res["pI_interval"]]),
    )
    mol.SetProp(
        "pI interval threshold",
        "%.2f" % pichemist_res["pI_interval_threshold"],
    )


def _append_pichemist_props_to_dict(dct, pichemist_res):
    dct["pI mean"] = "%.2f" % pichemist_res["pI"]["pI mean"]
    dct["pI std"] = "%.2f" % pichemist_res["pI"]["std"]
