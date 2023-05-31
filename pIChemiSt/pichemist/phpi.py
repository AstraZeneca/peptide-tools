import numpy as np

from pichemist.charges import PKaChargeCalculator


# Define pH span tocalcualte itration curve and where to search for pI.
def get_pH_span():
    pH_llim = -1
    pH_hlim = 15
    return [pH_llim, pH_hlim]


def CalcChargepHCurve(base_pkas, acid_pkas, diacid_pkas, constant_q=0):
    pH_limit = get_pH_span()
    pH_step=0.1
    pH_a = np.arange(pH_limit[0], pH_limit[1]+pH_step, pH_step)
    Q_a=pH_a*0.0    
    for i in range(len(pH_a)):
        Q = PKaChargeCalculator().calculate_charge(base_pkas, acid_pkas, diacid_pkas,
                                                   pH_a[i], constant_q=constant_q)
        Q_a[i]=Q
    pH_Q = np.vstack((pH_a,Q_a))
    return pH_Q


def calculateIsoelectricPoint(base_pkas, acid_pkas, diacid_pkas, constant_q=0):   
    tolerance=0.01
    charge_tol=0.05
    na=len(acid_pkas)+len(diacid_pkas)
    nb=len(base_pkas)
    
    pH_lim = get_pH_span()
    lower_pH = pH_lim[0] 
    higher_pH = pH_lim[1] 

    while True:
        mid_pH = 0.5 * (higher_pH + lower_pH)
        charge = PKaChargeCalculator().calculate_charge(base_pkas, acid_pkas, diacid_pkas, mid_pH, constant_q=constant_q)
        
        if na == 0 and nb != 0:
            #print "---!Warning: no acidic ionizable groups, only basic groups present in the sequence. pI is not defined and thus won't be calculated. However, you can still plot the titration curve. Continue."
            refcharge = charge_tol * nb

        elif nb == 0 and na != 0:
            #print "---!Warning: no basic ionizable groups, only acidic groups present in the sequence. pI is not defined and thus won't be calculated. However, you can still plot the titration curve. Continue."
            refcharge = -charge_tol * na

        else:
            refcharge = 0.0


        if charge > refcharge + tolerance:
            lower_pH = mid_pH
        elif charge < refcharge - tolerance:
            higher_pH = mid_pH
        else:
            return mid_pH
            
        if mid_pH <= pH_lim[0]:
            return pH_lim[0]
        elif mid_pH >= pH_lim[1]:
            return pH_lim[1]