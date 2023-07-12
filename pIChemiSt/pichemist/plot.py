import numpy as np
import matplotlib.pyplot as plt

from itertools import cycle
from pichemist.config import PKA_SETS_NAMES


### PLOT titration curve
def generate_and_save_titration_curve(pH_q_dict,fig_filename):
    plt.matplotlib.rcParams.update({'font.size': 16})
    lines = ["-","--","-.",":"]
    w1=4.0 ; w2=3.0 ; w3=2.0 ; w4=1.0
    linew = [w1,w1, w2,w2, w3,3, w4,w4]
    linecycler = cycle(lines)
    linewcycler = cycle(linew)

    plt.figure(figsize=(8,6))
    i=0
    for pka_set in PKA_SETS_NAMES:
        i+=1
        pH_Q = pH_q_dict[pka_set] 
        l = plt.plot(pH_Q[:,0],pH_Q[:,1],next(linecycler),label=pka_set,linewidth=next(linewcycler)) 
        if pka_set == 'IPC2_peptide': 
            plt.setp(l,linewidth=8,linestyle='-',color='k')

        # Store data for output
        if i==1: 
            pH = pH_Q[:,0]
            Q_M = pH_Q[:,1]
        else:
            Q_M = np.column_stack([Q_M,pH_Q[:,1]])

    plt.plot(pH,pH*0,'k-')
    plt.plot([7,7],[np.min(Q_M),np.max(Q_M)],'k-')
    plt.xlim([2,12])
    plt.ylim([np.min(Q_M),np.max(Q_M)])
    plt.legend(loc="center right", bbox_to_anchor=[1.1, 0.5],ncol=1, shadow=True, fontsize=10).get_frame().set_alpha(1)
    plt.ylabel('peptide charge')
    plt.xlabel('pH')
    plt.minorticks_on()
    plt.grid(True)
    plt.savefig(fig_filename)
    return