from pichemist.charges import SmartsChargeCalculator
from pichemist.config import REFERENCE_PKA_SET
from pichemist.core import calculate_isoelectric_interval_and_threshold
from pichemist.core import calculate_pI_pH_and_charge_dicts
from pichemist.core import merge_matched_and_calculated_pkas
from pichemist.fasta.matcher import FastaPKaMatcher
from pichemist.model import InputAttribute
from pichemist.model import MODELS
from pichemist.model import OutputAttribute
from pichemist.model import PKaMethod
from pichemist.molecule import MolStandardiser
from pichemist.molecule import PeptideCutter
from pichemist.pka.acd import ACDPKaCalculator
from pichemist.pka.pkamatcher import PKaMatcher
from pichemist.plot import output_ph_q_curve
from rdkit import Chem
from rdkit.Chem import Draw


class ApiException(Exception):
    pass


def fasta_pkas_from_list(smiles_list):
    """Match pKas from FASTA definitions for a SMILES list."""
    return FastaPKaMatcher().get_aa_pkas_from_list(smiles_list)


def calculated_pkas_from_list(smiles_list, method):
    """Calculates pKa values using ACD or pKa Matcher."""
    if method not in MODELS[PKaMethod]:
        raise NotImplementedError(
            "Invalid method. Only the formats " f"{MODELS[PKaMethod]} are accepted"
        )
    if method == PKaMethod.ACD.value:
        base_pkas, acid_pkas, diacid_pkas = ACDPKaCalculator().calculate_pka_from_list(
            smiles_list
        )
    if method == PKaMethod.PKA_MATCHER.value:
        base_pkas, acid_pkas, diacid_pkas = PKaMatcher().calculate_pka_from_list(
            smiles_list
        )
    return base_pkas, acid_pkas, diacid_pkas


def calc_frags_for_output_fasta(ionization_type, pkas_fasta):

    D_pka = dict()
    D_count = dict()
    pka_sets_cnt = 0
    for pka_set, list_for_pka_set in pkas_fasta.items():
        pka_sets_cnt += 1
        for v in list_for_pka_set:
            pka = v[0]
            AA = v[1]
            if AA in D_pka.keys():
                D_pka[AA].append(pka)
            else:
                D_pka[AA] = list()

            if pka_sets_cnt == 1:
                if AA in D_count.keys():
                    D_count[AA] += 1
                else:
                    D_count[AA] = 1

    frag_pkas_fasta = dict()
    idx = 0
    for k, v in D_pka.items():
        idx += 1
        frag_pkas_fasta[idx] = {
            "type": ionization_type,
            "frag": k,
            "count": D_count[k],
            "pka": sum(v) / len(v),
        }

    return frag_pkas_fasta


def smiles_to_image(smiles, image_path, image_size=(300, 300)):
    """
    Convert a SMILES string to a chemical structure image with a specified size.

    Parameters:
    - smiles (str): The SMILES string representing the molecule.
    - image_path (str): The file path where the image will be saved.
    - image_size (tuple): The size of the image (width, height).

    Returns:
    - None
    """
    try:
        # Convert SMILES to a molecule object
        molecule = Chem.MolFromSmiles(smiles)

        # Check if the molecule was successfully created
        if molecule is None:
            raise ValueError("Invalid SMILES string.")

        # Generate the image with the specified size
        image = Draw.MolToImage(molecule, size=image_size)

        # Save the image to the specified file path
        image.save(image_path)
        # print(f"Image saved successfully to {image_path}")
    except Exception as e:
        print(f"An error occurred: {e}")


def calc_frags_for_output_calc(ionization_type, pkas_calc):

    # D_pka = dict()
    # D_count = dict()
    # for k2,v2 in pkas_calc:
    #     pka = v2[0]
    #     smi = v2[1]
    #     D_pka[smi].append(pka)
    #     if smi in D_count.keys():
    #         D_count[smi]+=1
    #     else:
    #         D_count[smi]=1

    # frag_pkas_calc = dict()
    # for k,v in D_pka.items():
    #     frag_pkas_calc[k] = {'type':ionization_type,'frag':k, 'count':D_count[k], 'pka':sum(v) / len(v)}

    # TODO need to handle cases with multiple ionization in the same fragment. Apparently could only be done
    # TODO if the index of fragment or alike is stored in pka_calc dictionary. Requires, some code refurbishment.
    # Solution for now is to do not average out pka of identical fragment and display all with count 1
    frag_pkas_calc = dict()
    frg_idx = 0
    for v in pkas_calc:
        frg_idx += 1
        pka = v[0]
        smi = v[1]
        # frag_pkas_calc[frg_idx] = {'type':ionization_type,'frag':smi, 'count':1, 'pka':pka}

        image_path = "fragment_" + str(frg_idx) + ".png"
        image_size = (500, 500)  # Set the desired image size
        smiles_to_image(smi, image_path, image_size)
        frag_pkas_calc[frg_idx] = {
            "type": ionization_type,
            "frag": image_path,
            "count": 1,
            "pka": pka,
        }

    return frag_pkas_calc


def compile_frags_pkas_for_output(
    base_pkas_fasta,
    acid_pkas_fasta,
    diacid_pkas_fasta,
    base_pkas_calc,
    acid_pkas_calc,
    diacid_pkas_calc,
    net_qs_and_frags,
):
    """
    Produces dictionary with fragmets (known AA or smiles fragment), their occurences in the molecule, corresponding pKa
    (average between pKa sets in case of known AA)

    """
    frag_acid_pkas_fasta = calc_frags_for_output_fasta("acid", acid_pkas_fasta)

    frag_base_pkas_fasta = calc_frags_for_output_fasta("base", base_pkas_fasta)

    frag_acid_pkas_calc = calc_frags_for_output_calc("acid", acid_pkas_calc)

    frag_base_pkas_calc = calc_frags_for_output_calc("base", base_pkas_calc)

    frag_Qs_calc = calc_frags_for_output_calc("constant charge", net_qs_and_frags)

    # TODO: Diacid dictionaries are not used, they are
    # deprecated and should be removed from the code
    return (
        frag_acid_pkas_fasta,
        frag_base_pkas_fasta,
        frag_acid_pkas_calc,
        frag_base_pkas_calc,
        frag_Qs_calc,
    )


def pkas_and_charges_from_list(smiles_list, method):
    """
    Produces the pKa values for a list of SMILES by matching them against
    some FASTA definitions and calculating the unmatched ones.

    """
    # FASTA match
    (
        unknown_frags,
        base_pkas_fasta,
        acid_pkas_fasta,
        diacid_pkas_fasta,
    ) = fasta_pkas_from_list(smiles_list)

    # Unknown fragment calculation
    base_pkas_calc, acid_pkas_calc, diacid_pkas_calc = calculated_pkas_from_list(
        unknown_frags, method=method
    )

    # Calculate charges and their matching fragments, and return results
    net_qs_and_frags = SmartsChargeCalculator().calculate_net_qs_from_list(
        unknown_frags
    )
    return (
        base_pkas_fasta,
        acid_pkas_fasta,
        diacid_pkas_fasta,
        base_pkas_calc,
        acid_pkas_calc,
        diacid_pkas_calc,
        net_qs_and_frags,
    )


def pichemist_from_dict(
    input_dict,
    method,
    ph_q_curve_file_prefix=None,
    plot_ph_q_curve=False,
    print_fragments=False,
):
    """Runs the full logic for a given input dictionary."""
    dict_output = dict()
    for mol_idx in input_dict.keys():

        # Prepare molecule and break into fragments
        mol_name = input_dict[mol_idx][InputAttribute.MOL_NAME.value]
        mol = input_dict[mol_idx][InputAttribute.MOL_OBJECT.value]
        mol = MolStandardiser().standardise_molecule(mol)
        smiles_list = PeptideCutter().break_amide_bonds_and_cap(mol)

        # Calculate pKas and charges
        (
            base_pkas_fasta,
            acid_pkas_fasta,
            diacid_pkas_fasta,
            base_pkas_calc,
            acid_pkas_calc,
            diacid_pkas_calc,
            net_qs_and_frags,
        ) = pkas_and_charges_from_list(smiles_list, method)
        (
            base_pkas_dict,
            acid_pkas_dict,
            diacid_pkas_dict,
        ) = merge_matched_and_calculated_pkas(
            base_pkas_fasta,
            base_pkas_calc,
            acid_pkas_fasta,
            acid_pkas_calc,
            diacid_pkas_fasta,
            diacid_pkas_calc,
        )

        # Recomple fragments and pKa for table output
        (
            frag_acid_pkas_fasta,
            frag_base_pkas_fasta,
            frag_acid_pkas_calc,
            frag_base_pkas_calc,
            frag_Qs_calc,
        ) = compile_frags_pkas_for_output(
            base_pkas_fasta,
            acid_pkas_fasta,
            diacid_pkas_fasta,
            base_pkas_calc,
            acid_pkas_calc,
            diacid_pkas_calc,
            net_qs_and_frags,
        )

        # Calculate the curves
        pI_dict, q_dict, pH_q_dict = calculate_pI_pH_and_charge_dicts(
            base_pkas_dict, acid_pkas_dict, diacid_pkas_dict, net_qs_and_frags
        )

        # Calculate isoelectric interval
        interval, interval_threshold = calculate_isoelectric_interval_and_threshold(
            pH_q_dict
        )

        # Output for given molecule
        dict_output[mol_idx] = {
            OutputAttribute.MOL_NAME.value: mol_name,
            OutputAttribute.PI.value: pI_dict,
            OutputAttribute.Q_PH7.value: q_dict,
            OutputAttribute.PI_INTERVAL.value: interval,
            OutputAttribute.PI_INTERVAL_THRESHOLD.value: interval_threshold,
        }

        # Plot pH/Q curve
        ###AIF###if plot_ph_q_curve and not ph_q_curve_file_prefix:
        if plot_ph_q_curve and not isinstance(ph_q_curve_file_prefix, str):
            raise ApiException(
                "A file prefix for the pH/Q curve plots must " "be specified."
            )
        if plot_ph_q_curve:
            fig_filename = f"{ph_q_curve_file_prefix}_{mol_idx}.png"
            output_ph_q_curve(pH_q_dict, fig_filename)
            dict_output[mol_idx][OutputAttribute.PLOT_FILENAME.value] = fig_filename

        # Define set for reporting pKa of individual amino acids and fragments
        dict_output[mol_idx].update({OutputAttribute.PKA_SET.value: REFERENCE_PKA_SET})

        if print_fragments:
            # No need to include diacids pkas as they
            # are included as single apparent ionisations
            dict_output[mol_idx].update(
                {
                    OutputAttribute.BASE_PKA_FASTA.value: base_pkas_fasta,
                    OutputAttribute.ACID_PKA_FASTA.value: acid_pkas_fasta,
                    OutputAttribute.BASE_PKA_CALC.value: base_pkas_calc,
                    OutputAttribute.ACID_PKA_CALC.value: acid_pkas_calc,
                    OutputAttribute.CONSTANT_QS.value: net_qs_and_frags,
                }
            )

            dict_output[mol_idx].update(
                {
                    OutputAttribute.FRAG_BASE_PKA_FASTA.value: frag_base_pkas_fasta,
                    OutputAttribute.FRAG_ACID_PKA_FASTA.value: frag_acid_pkas_fasta,
                    OutputAttribute.FRAG_BASE_PKA_CALC.value: frag_base_pkas_calc,
                    OutputAttribute.FRAG_ACID_PKA_CALC.value: frag_acid_pkas_calc,
                    OutputAttribute.FRAG_CONSTANT_QS.value: frag_Qs_calc,
                }
            )
    return dict_output
