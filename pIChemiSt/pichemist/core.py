from pichemist.charges import PKaChargeCalculator
from pichemist.fasta.matcher import FastaPKaMatcher
from pichemist.isoelectric import CurveCalculator
from pichemist.isoelectric import IsoelectricCalculator
from pichemist.molecule import smiles_to_image
from pichemist.pka.utils import merge_pkas_lists
from pichemist.stats import mean
from pichemist.stats import stddev
from pichemist.stats import stderr


def get_net_qs_from_qs_and_frags(net_qs_and_frags):
    """Returns a list of charges from charge-fragment tuples."""
    return [v[0] for v in net_qs_and_frags]


def generate_stats_on_dict(input_dict, mean_title="mean"):
    """
    Enriches the an input dictionary of float values
    with mean, std deviation, and std error.

    """
    value_list = list()
    for k in input_dict.keys():
        value_list += [input_dict[k]]
    input_dict[mean_title] = mean(value_list)
    input_dict["std"] = stddev(value_list)
    input_dict["err"] = stderr(value_list)
    return input_dict


def get_low_and_high_from_interval_lists(interval_low_list, interval_high_list):
    """
    Processes interval lists to yield their averages.
    If molecule permanently has a charge, the interval is
    not defined and NaN are provided instead.

    """
    if len(interval_low_list) > 0:
        interval_low = mean(interval_low_list)
    else:
        interval_low = float("NaN")
    if len(interval_high_list) > 0:
        interval_high = mean(interval_high_list)
    else:
        interval_high = float("NaN")
    return interval_low, interval_high


def merge_matched_and_calculated_pkas(
    base_pkas_fasta,
    base_pkas_calc,
    acid_pkas_fasta,
    acid_pkas_calc,
    diacid_pkas_fasta,
    diacid_pkas_calc,
):
    """Merge FASTA-matched and calculated pKas."""
    base_pkas_list = [base_pkas_fasta, base_pkas_calc]
    acid_pkas_list = [acid_pkas_fasta, acid_pkas_calc]
    diacid_pkas_list = [diacid_pkas_fasta, diacid_pkas_calc]
    base_pkas_dict, acid_pkas_dict, diacid_pkas_dict = merge_pkas_lists(
        base_pkas_list, acid_pkas_list, diacid_pkas_list
    )
    return base_pkas_dict, acid_pkas_dict, diacid_pkas_dict


def calculate_pI_pH_and_charge_dicts(
    base_pkas_dict, acid_pkas_dict, diacid_pkas_dict, net_qs_and_frags
):
    """Calculates the isoelectric point, charges, and pH charges."""
    pI_dict = dict()
    q_dict = dict()
    pH_q_dict = dict()
    pka_sets_names = FastaPKaMatcher().get_pka_sets_names()
    for pka_set in pka_sets_names:

        # Merge fasta and calculated pkas
        base_pkas = base_pkas_dict[pka_set]
        acid_pkas = acid_pkas_dict[pka_set]
        diacid_pkas = diacid_pkas_dict[pka_set]

        # Calculate isoelectric point
        net_qs = get_net_qs_from_qs_and_frags(net_qs_and_frags)
        constant_q = PKaChargeCalculator().calculate_constant_charge(net_qs)
        q = PKaChargeCalculator().calculate_charge(
            base_pkas, acid_pkas, diacid_pkas, pH=7.4, constant_q=constant_q
        )
        pI = IsoelectricCalculator().calculate_pI(
            base_pkas, acid_pkas, diacid_pkas, constant_q=constant_q
        )
        pH_q = CurveCalculator().calculate_charged_curve(
            base_pkas, acid_pkas, diacid_pkas, constant_q=constant_q
        )
        pI_dict[pka_set] = pI
        q_dict[pka_set] = q
        pH_q_dict[pka_set] = pH_q

    # Generate stats
    pI_dict = generate_stats_on_dict(pI_dict, mean_title="pI mean")
    q_dict = generate_stats_on_dict(q_dict, mean_title="Q at pH7.4 mean")
    return pI_dict, q_dict, pH_q_dict


def calculate_isoelectric_interval_and_threshold(pH_q_dict):
    """
    Calculates the isoelectric interval that is
    the pH range where the charge is within the given threshold.

    """
    threshold = IsoelectricCalculator().interval_threshold
    interval_low_list = list()
    interval_high_list = list()

    # Calculate intervals across all sets
    for pka_set in FastaPKaMatcher().get_pka_sets_names():
        pH_int = IsoelectricCalculator().calculate_interval(pH_q_dict[pka_set])
        if len(pH_int) > 1:
            interval_low_list.append(pH_int[0])
            interval_high_list.append(pH_int[-1])

    # Average and return results
    interval_low, interval_high = get_low_and_high_from_interval_lists(
        interval_low_list, interval_high_list
    )
    return (interval_low, interval_high), threshold


def calculate_frags_for_output_calc(
    ionization_type, pkas_calc, generate_fragment_images=False
):
    """
    Generates a dictionary of fragments with their pKas.
    Images can additionally be included.

    """
    frag_pkas_calc = dict()
    frg_idx = 0
    for v in pkas_calc:
        frg_idx += 1
        pka = v[0]
        smi = v[1]

        frag_pkas_calc[frg_idx] = {
            "type": ionization_type,
            "frag": smi,
            "count": 1,
            "pka": pka,
        }

        base64_image = None
        if generate_fragment_images:
            frag_pkas_calc[frg_idx]["base64_image"] = smiles_to_image(smi)
    return frag_pkas_calc


def calculate_frags_for_output_fasta(ionization_type, pkas_fasta):
    """Generates a dictionary of fragments with their pKas."""
    pka_dict = dict()
    count_dict = dict()
    pka_sets_count = 0
    for pka_set, list_for_pka_set in pkas_fasta.items():
        pka_sets_count += 1
        for v in list_for_pka_set:
            pka = v[0]
            aa = v[1]
            if aa in pka_dict.keys():
                pka_dict[aa].append(pka)
            else:
                pka_dict[aa] = list()
            if pka_sets_count == 1:
                if aa in count_dict.keys():
                    count_dict[aa] += 1
                else:
                    count_dict[aa] = 1

    frag_pkas_fasta = dict()
    idx = 0
    for k, v in pka_dict.items():
        idx += 1
        frag_pkas_fasta[idx] = {
            "type": ionization_type,
            "frag": k,
            "count": count_dict[k],
            "pka": sum(v) / len(v),
        }
    return frag_pkas_fasta
