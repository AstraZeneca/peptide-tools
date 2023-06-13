import numpy as np

from pichemist.charges import PKaChargeCalculator


class CurveCalculator(object):
    """Uses pKa values to calculate the pH/Q curve."""
    def __init__(self):
        self.pH_lower_bound = -1
        self.pH_upper_bound = 15
        self.pH_step = 0.1
        self.pH_range = self._define_pH_range()
        self.q_range = self._define_charge_range()

    def _define_pH_range(self):
        """Sets the pH range for the curve (X axis)."""
        pH_limit = self.get_pH_span()
        return np.arange(pH_limit[0], pH_limit[1]+self.pH_step, self.pH_step)

    def _define_charge_range(self):
        """Sets the charge range (Y axis)."""
        return self.pH_range*0.0

    def get_pH_span(self):
        """
        Define pH span to, for example, calculate the
        titration curve or where to search for pI.

        """
        return [self.pH_lower_bound, self.pH_upper_bound]

    def calculate_charged_curve(self, base_pkas, acid_pkas,
                                diacid_pkas, constant_q=0):
        """Calculates the pH/Q curve."""
        for i in range(len(self.pH_range)):
            charge = PKaChargeCalculator().calculate_charge(
                base_pkas, acid_pkas, diacid_pkas,
                self.pH_range[i], constant_q=constant_q)
            self.q_range[i] = charge
        return np.vstack((self.pH_range, self.q_range))


def calculateIsoelectricPoint(base_pkas, acid_pkas, diacid_pkas, constant_q=0):
    tolerance=0.01
    charge_tol=0.05
    na=len(acid_pkas)+len(diacid_pkas)
    nb=len(base_pkas)
    
    pH_lim = CurveCalculator().get_pH_span()
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
