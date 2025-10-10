import numpy as np
import matplotlib.pyplot as plt

from itertools import cycle
from pichemist.config import PKA_SETS_NAMES
from pichemist.config import PLOT_LINE_WIDTHS as pld
from pichemist.config import REFERENCE_PKA_SET


def output_ph_q_curve(ph_q_dict, fig_filename):
    """Generates a pH/Q curve plot and saves it."""
    # Static parameters
    plt.matplotlib.rcParams.update({"font.size": 16})
    plt_image = plt.figure(figsize=(8, 6))
    lines = ["-", "--", "-.", ":"]
    linew = [pld["w1"], pld["w1"],
             pld["w2"], pld["w2"],
             pld["w3"], pld["w2"],
             pld["w4"], pld["w4"]]
    linecycler = cycle(lines)
    linewcycler = cycle(linew)
    plt.ylabel("peptide charge")
    plt.xlabel("pH")
    plt.minorticks_on()
    plt.grid(True)

    # Plot curve for each pKa set
    i = 0
    for pka_set in PKA_SETS_NAMES:
        i += 1
        ph_q = ph_q_dict[pka_set]
        pl = plt.plot(ph_q[:, 0], ph_q[:, 1], next(linecycler),
                      label=pka_set, linewidth=next(linewcycler))
        if pka_set == REFERENCE_PKA_SET:
            plt.setp(pl, linewidth=8, linestyle="-", color="k")
        if i == 1:
            pH = ph_q[:, 0]
            q_m = ph_q[:, 1]
        else:
            q_m = np.column_stack([q_m, ph_q[:, 1]])
    plt.plot([7, 7], [np.min(q_m), np.max(q_m)], "k-")
    plt.plot(pH, pH*0, "k-")

    # Dynamic parameters
    plt.xlim([2, 12])
    plt.ylim([np.min(q_m), np.max(q_m)])
    plt.legend(loc="center right", bbox_to_anchor=[1.1, 0.5],
               ncol=1, shadow=True,
               fontsize=10).get_frame().set_alpha(1)
    plt.savefig(fig_filename)
    plt.close(plt_image)
