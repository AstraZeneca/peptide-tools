def _append_ext_coeff_props_to_mol(mol, ext_coeff_res):
    mol.SetProp("mol_name", ext_coeff_res["mol_name"])
    mol.SetProp("Sequence(FASTA)", ext_coeff_res["fasta"])
    mol.SetProp("e205(nm)", "%i" % ext_coeff_res["e205"])
    mol.SetProp("e214(nm)", "%i" % ext_coeff_res["e214"])
    mol.SetProp("e280(nm)", "%i" % ext_coeff_res["e280"])


def _append_ext_coeff_props_to_dict(dct, ext_coeff_res):
    dct["mol_name"] = ext_coeff_res["mol_name"]
    dct["Sequence(FASTA)"] = ext_coeff_res["fasta"]
    dct["e205(nm)"] = "%i" % ext_coeff_res["e205"]
    dct["e214(nm)"] = "%i" % ext_coeff_res["e214"]
    dct["e280(nm)"] = "%i" % ext_coeff_res["e280"]
