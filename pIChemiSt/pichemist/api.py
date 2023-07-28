from pichemist.config import REFERENCE_PKA_SET
from pichemist.charges import SmartsChargeCalculator
from pichemist.core import calculate_pI_pH_and_charge_dicts
from pichemist.core import calculate_isoelectric_interval_and_threshold
from pichemist.core import merge_matched_and_calculated_pkas
from pichemist.fasta.matcher import FastaPKaMatcher
from pichemist.molecule import MolStandardiser
from pichemist.molecule import PeptideCutter
from pichemist.pka.acd import ACDPKaCalculator
from pichemist.pka.pkamatcher import PKaMatcher
from pichemist.model import PKaMethod
from pichemist.model import OutputAttributes
from pichemist.model import MODELS
from pichemist.plot import output_titration_curve


class ApiException(Exception):
    pass


def fasta_pkas_from_list(smiles_list):
    """Match pKas from FASTA definitions for a SMILES list."""
    return FastaPKaMatcher().get_aa_pkas_from_list(smiles_list)


def calculated_pkas_from_list(smiles_list, method):
    """Calculates pKa values using ACD or pKa Matcher."""
    if method not in MODELS[PKaMethod]:
        raise NotImplementedError("Invalid method. Only the formats "
                                  f"{MODELS[PKaMethod]} are accepted")
    if method == PKaMethod.ACD.value:
        base_pkas, acid_pkas, diacid_pkas = \
            ACDPKaCalculator().calculate_pka_from_list(smiles_list)
    if method == PKaMethod.PKA_MATCHER.value:
        base_pkas, acid_pkas, diacid_pkas = \
            PKaMatcher().calculate_pka_from_list(smiles_list)
    return base_pkas, acid_pkas, diacid_pkas


def pkas_and_charges_from_list(smiles_list, method):
    """
    Produces the pKa values for a list of SMILES by matching them against
    some FASTA definitions and calculating the unmatched ones.

    """
    # FASTA match
    unknown_frags, base_pkas_fasta, acid_pkas_fasta, diacid_pkas_fasta = \
        fasta_pkas_from_list(smiles_list)

    # Unknown fragment calculation
    base_pkas_calc, acid_pkas_calc, diacid_pkas_calc = \
        calculated_pkas_from_list(unknown_frags,
                                  method=method)

    # Calculate charges and their matching fragments, and return results
    net_qs_and_frags = SmartsChargeCalculator().calculate_net_qs_from_list(
        unknown_frags)
    return base_pkas_fasta, acid_pkas_fasta, diacid_pkas_fasta, \
        base_pkas_calc, acid_pkas_calc, diacid_pkas_calc, net_qs_and_frags


def pichemist_from_list(input_dict, method,
                        titration_file_prefix=None,
                        plot_titration_curve=False,
                        print_fragments=False):
    """Runs the full logic for a given input dictionary."""
    dict_output = dict()
    for mol_idx in input_dict.keys():

        # Prepare molecule and break into fragments
        mol_name = input_dict[mol_idx]["mol_name"]
        mol = input_dict[mol_idx]["mol_obj"]
        mol = MolStandardiser().standardise_molecule(mol)
        smiles_list = PeptideCutter().break_amide_bonds_and_cap(mol)

        # Calculate pKas and charges
        base_pkas_fasta, acid_pkas_fasta, diacid_pkas_fasta, base_pkas_calc, \
            acid_pkas_calc, diacid_pkas_calc, net_qs_and_frags = \
            pkas_and_charges_from_list(smiles_list, method)
        base_pkas_dict, acid_pkas_dict, diacid_pkas_dict = \
            merge_matched_and_calculated_pkas(
                base_pkas_fasta, base_pkas_calc,
                acid_pkas_fasta, acid_pkas_calc,
                diacid_pkas_fasta, diacid_pkas_calc)

        # Calculate the curves
        pI_dict, q_dict, pH_q_dict = calculate_pI_pH_and_charge_dicts(
            base_pkas_dict, acid_pkas_dict,
            diacid_pkas_dict, net_qs_and_frags)

        # Calculate isoelectric interval
        interval, interval_threshold = \
            calculate_isoelectric_interval_and_threshold(pH_q_dict)

        # Plot titration curve
        if plot_titration_curve and not titration_file_prefix:
            raise ApiException("A file prefix for the titration plots must "
                               "be specified.")
        # TODO: ANDREY - fig_filename can be removed from dict_output
        # when the plot is not generated
        fig_filename = ""
        if plot_titration_curve:
            fig_filename = f"{titration_file_prefix}_{mol_idx}"
            output_titration_curve(pH_q_dict, fig_filename)

        # Output for given molecule
        dict_output[mol_idx] = {
            OutputAttributes.MOL_NAME.value: mol_name,
            OutputAttributes.PI.value: pI_dict,
            OutputAttributes.Q_PH7.value: q_dict,
            OutputAttributes.PI_INTERVAL.value: interval,
            OutputAttributes.PI_INTERVAL_THRESHOLD.value: interval_threshold,
            OutputAttributes.PLOT_FILENAME.value: fig_filename}

        # Define set for reporting pKa of individual amino acids and fragments
        dict_output[mol_idx].update(
            {OutputAttributes.PKA_SET.value: REFERENCE_PKA_SET})

        if print_fragments:
            # No need to include diacids pkas as they
            # are included as single apparent ionisations
            dict_output[
                mol_idx].update({
                    OutputAttributes.BASE_PKA_FASTA.value: base_pkas_fasta,
                    OutputAttributes.ACID_PKA_FASTA.value: acid_pkas_fasta,
                    OutputAttributes.BASE_PKA_CALC.value: base_pkas_calc,
                    OutputAttributes.ACID_PKA_CALC.value: acid_pkas_calc,
                    OutputAttributes.CONSTANT_QS.value: net_qs_and_frags
                    })
    return dict_output
