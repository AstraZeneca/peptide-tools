#!/usr/bin/python
# -*- coding: utf-8 -*
import argparse
import csv
import json
import os
import random
import re
import string
import sys
import tempfile
import time
import urllib
from operator import itemgetter
from time import gmtime
from time import strftime

from pichemist.api import pichemist_from_dict
from pichemist.io import generate_input
from rdkit import Chem

from extn_coeff_fasta import calc_extn_coeff
from pI_fasta import calc_pI_fasta


currentdir = os.getcwd()

__prog__ = "Peptide Tools Master"
__doc__ = """TODO"""


def arg_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Input filepath (or SMILES string or FASTA)",
        default=None,
    )

    ### pIChemiSt.py keys
    parser.add_argument(
        "--print_fragment_pkas",
        default=False,
        action="store_true",
        dest="print_fragment_pkas",
        help="Print the fragments with corresponding pKas used in pI calcution.",
    )
    parser.add_argument(
        "--print_pka_set",
        default=False,
        action="store_true",
        dest="print_pka_set",
        help="Print out stored pka sets.",
    )

    ### pI_fasta.py keys
    parser.add_argument(
        "--ionized_Cterm",
        dest="ionized_Cterm",
        help="is C-terminus ionized [COO-]?",
        default=True,
    )
    parser.add_argument(
        "--ionized_Nterm",
        dest="ionized_Nterm",
        help="is N-terminus ionized [N+]?",
        default=True,
    )
    parser.add_argument(
        "-p",
        action="store",
        dest="NPhosphateGroups",
        help="Number of phosphorilated residues. Phosphorilated residues must be denoted as X in the sequence. default = 0",
        default=0,
        type=int,
    )
    parser.add_argument(
        "-l",
        action="store",
        dest="NAlkylLysGroups",
        help="Number of monoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0",
        default=0,
        type=int,
    )
    parser.add_argument(
        "-d",
        action="store",
        dest="NDiAlkylLysGroups",
        help="Number of dinoalkylated Lys residues. These residues should be denoted as X in the sequence. default = 0",
        default=0,
        type=int,
    )

    args = parser.parse_args()
    return args


def read_fasta_file(inputFile):

    filename, ext = os.path.splitext(inputFile)

    # Initialize file reader
    if not ext == ".fasta":
        raise Exception(
            '!Warning: extension of file is not ".fasta". Assuming it is fasta formatted input. Continue. '
        )

    from Bio import SeqIO

    biosuppl = SeqIO.parse(open(inputFile), "fasta")

    mol_supply_json = {}
    mol_unique_ID = 0
    for biofasta in biosuppl:
        mol_unique_ID += 1
        # unique index, mol title, RDkit mol object, mol fasta
        mol_supply_json[mol_unique_ID] = {
            "mol_name": biofasta.id,
            "mol_obj": None,
            "fasta": str(biofasta.seq),
        }

    return mol_supply_json


if __name__ == "__main__":
    args = arg_parser()
    INPUT = args.input

    lPlot = True
    from peptools.io import generate_input

    (
        mol_supply_json,
        inputFile,
        outputFile,
        l_calc_extn_coeff,
        l_calc_pI_fasta,
        l_calc_pIChemiSt,
        file_extension,
        out_fext,
    ) = generate_input(INPUT)

    # print("mol_name: "+mol_name)
    # print("FASTA: "+fasta)
    # print("SMILES: "+smi)

    dict_out_extn_coeff_fasta = {}
    if l_calc_extn_coeff:
        extn_coeff_options = {
            "seq": "",
            "inputDict": mol_supply_json,
            "inputJSON": "",
            "inputFile": "",
            "outputFile": "",
            "l_json": True,
        }
        dict_out_extn_coeff = calc_extn_coeff(extn_coeff_options)

    dict_out_pI_fasta = {}
    if l_calc_pI_fasta:
        # prepare pI_fasta predictor

        if not args.ionized_Cterm:
            IonizableTerminiOfCTermRes = "''"
        else:
            IonizableTerminiOfCTermRes = "_"

        if not args.ionized_Nterm:
            IonizableTerminiOfNTermRes = "''"
        else:
            IonizableTerminiOfNTermRes = "_"

        pI_fasta_options = {
            "seq": "",
            "inputDict": mol_supply_json,
            "inputJSON": "",
            "inputFile": "",
            "outputFile": "",
            "tol": 0.001,
            "CTermRes": "_",
            "NTermRes": "_",
            "IonizableTerminiOfCTermRes": IonizableTerminiOfCTermRes,
            "IonizableTerminiOfNTermRes": IonizableTerminiOfNTermRes,
            "lCyclic": False,
            "NPhosphateGroups": args.NPhosphateGroups,
            "NAlkylLysGroups": args.NAlkylLysGroups,
            "NDiAlkylLysGroups": args.NDiAlkylLysGroups,
            "lPrintpKa": False,
            "lPlot": lPlot,
            "lIgnoreC": False,
            "plot_filename": "OUT_titration_curve.png",
            "l_json": True,
            "pka_set_list": "",
        }

        dict_out_pI_fasta = calc_pI_fasta(pI_fasta_options)

    dict_out_pIChemiSt = {}
    if l_calc_pIChemiSt:
        args = {
            "plot_ph_q_curve": lPlot,
            "print_fragments": args.print_fragment_pkas,
            "method": "pkamatcher",
        }

        dict_out_pIChemiSt = pichemist_from_dict(
            mol_supply_json,
            args["method"],
            args["plot_ph_q_curve"],
            args["print_fragments"],
        )

    dict_out_peptide_tools_master = {
        "output_extn_coeff": dict_out_extn_coeff,
        "output_pI_fasta": dict_out_pI_fasta,
        "output_pIChemiSt": dict_out_pIChemiSt,
    }

    ### ----------------------------------------------------------------------
    # Output
    if outputFile == "":  # output JSON
        print(json.dumps(dict_out_peptide_tools_master, indent=2))

    else:  # output file
        if file_extension != ".fasta":

            # for mi in mol_supply_json.keys():
            mol_list = []
            for mi in mol_supply_json.keys():
                mol = mol_supply_json[mi]["mol_obj"]

                if l_calc_pIChemiSt:
                    mol.SetProp(
                        "pI mean", "%.2f" % dict_out_pIChemiSt[mi]["pI"]["pI mean"]
                    )
                    mol.SetProp("pI std", "%.2f" % dict_out_pIChemiSt[mi]["pI"]["std"])
                    mol.SetProp(
                        "pI interval",
                        " - ".join(
                            ["%.2f" % x for x in dict_out_pIChemiSt[mi]["pI_interval"]]
                        ),
                    )
                    mol.SetProp(
                        "pI interval threshold",
                        "%.2f" % dict_out_pIChemiSt[mi]["pI_interval_threshold"],
                    )

                if l_calc_extn_coeff:
                    mol.SetProp("mol_name", dict_out_extn_coeff[mi]["mol_name"])
                    mol.SetProp("Sequence(FASTA)", dict_out_extn_coeff[mi]["fasta"])
                    mol.SetProp("e205(nm)", "%i" % dict_out_extn_coeff[mi]["e205"])
                    mol.SetProp("e214(nm)", "%i" % dict_out_extn_coeff[mi]["e214"])
                    mol.SetProp("e280(nm)", "%i" % dict_out_extn_coeff[mi]["e280"])

                mol_list.append(mol)

            if out_fext == ".sdf":
                with Chem.SDWriter(outputFile) as sdf_w:
                    for mol in mol_list:
                        sdf_w.write(mol)

            elif out_fext == ".csv":
                with open(outputFile, "w") as csv_f:
                    csv_w = csv.writer(csv_f)
                    count = 0
                    for mol in mol_list:
                        props = mol.GetPropsAsDict()

                        count += 1
                        if count == 1:
                            header = ["SMILES"] + list(props.keys())
                            csv_w.writerow(header)

                        row = [Chem.MolToSmiles(mol)]
                        for p in header[1:]:
                            row += [props[p]]
                        csv_w.writerow(row)

        else:

            dict_list = []
            for mi in mol_supply_json.keys():
                fasta = mol_supply_json[mi]["fasta"]

                D = {}

                if l_calc_pI_fasta:
                    D["pI mean"] = "%.2f" % dict_out_pI_fasta[mi]["pI"]["pI mean"]
                    D["pI std"] = "%.2f" % dict_out_pI_fasta[mi]["pI"]["std"]

                if l_calc_extn_coeff:
                    D["mol_name"] = dict_out_extn_coeff[mi]["mol_name"]
                    D["Sequence(FASTA)"] = dict_out_extn_coeff[mi]["fasta"]
                    D["e205(nm)"] = "%i" % dict_out_extn_coeff[mi]["e205"]
                    D["e214(nm)"] = "%i" % dict_out_extn_coeff[mi]["e214"]
                    D["e280(nm)"] = "%i" % dict_out_extn_coeff[mi]["e280"]

                dict_list.append(D)

            if out_fext == ".csv":
                with open(outputFile, "w") as csv_f:
                    csv_w = csv.writer(csv_f)
                    count = 0
                    for props in dict_list:

                        count += 1
                        if count == 1:
                            header = list(props.keys())
                            csv_w.writerow(header)

                        row = []
                        for p in header:
                            row += [props[p]]
                        csv_w.writerow(row)

        dict_file = {
            "outputFile": outputFile,
            "outputInfo": "Number of molecules processed:"
            + str(len(mol_supply_json.keys())),
        }
        print(json.dumps(dict_file))
