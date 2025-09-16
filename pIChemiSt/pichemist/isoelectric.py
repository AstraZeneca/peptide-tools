import numpy as np
from pichemist.charges import PKaChargeCalculator
from pichemist.config import ROUNDING_DIGITS
from pichemist.utils import get_logger

log = get_logger(__name__)


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
        return np.arange(pH_limit[0], pH_limit[1] + self.pH_step, self.pH_step)

    def _define_charge_range(self):
        """Sets the charge range (Y axis)."""
        return self.pH_range * 0.0

    def get_pH_span(self):
        """
        Define pH span to, for example, calculate the
        pH/Q curve or where to search for pI.

        """
        return [self.pH_lower_bound, self.pH_upper_bound]

    def calculate_charged_curve(self, base_pkas, acid_pkas, constant_q=0):
        """Calculates the pH/Q curve."""
        for i in range(len(self.pH_range)):
            charge = PKaChargeCalculator().calculate_charge(
                base_pkas, acid_pkas, self.pH_range[i], constant_q=constant_q
            )
            self.q_range[i] = round(charge, ROUNDING_DIGITS)
        return np.vstack((self.pH_range, self.q_range)).T


class IsoelectricCalculator(object):
    """Calculates the isoelectric point and interval."""

    def __init__(self):
        self.tolerance = 0.01
        self.charge_tolerance = 0.05
        self.interval_threshold = 0.2
        self.pH_limit = CurveCalculator().get_pH_span()
        self.lower_pH = self.pH_limit[0]
        self.higher_pH = self.pH_limit[1]

    def calculate_pI(self, base_pkas, acid_pkas, constant_q=0):
        """
        Uses the pKas and charge to iteratively calculate
        the isoelectric point of a molecule.

        """
        while True:
            self.middle_pH = 0.5 * (self.higher_pH + self.lower_pH)
            charge = PKaChargeCalculator().calculate_charge(
                base_pkas, acid_pkas, self.middle_pH, constant_q=constant_q
            )
            na = len(acid_pkas)
            nb = len(base_pkas)

            if na == 0 and nb != 0:
                log.debug(
                    "Warning: no acidic ionizable groups, "
                    "only basic groups present in the "
                    "sequence. pI is not defined and thus "
                    "won't be calculated"
                    ""
                )
                reference_charge = self.charge_tolerance * nb
            elif nb == 0 and na != 0:
                log.debug(
                    "Warning: no basic ionizable groups, "
                    "only acidic groups present in the "
                    "sequence. pI is not defined and thus "
                    "won't be calculated"
                )
                reference_charge = -self.charge_tolerance * na
            else:
                reference_charge = 0.0

            # sic - instance attributes are updates as the
            # while block keeps looping
            if charge > reference_charge + self.tolerance:
                self.lower_pH = self.middle_pH
            elif charge < reference_charge - self.tolerance:
                self.higher_pH = self.middle_pH
            else:
                return self.middle_pH

            # If middle pH ends up outside the bounds
            # then return the bounds instead
            if self.middle_pH <= self.pH_limit[0]:
                return self.pH_limit[0]
            elif self.middle_pH >= self.pH_limit[1]:
                return self.pH_limit[1]

    def calculate_interval(self, pH_q):
        """Calculates the isoelectric interval."""
        q = pH_q[:, 1]
        pH = pH_q[:, 0]
        return pH[(q > -self.interval_threshold) & (q < self.interval_threshold)]
