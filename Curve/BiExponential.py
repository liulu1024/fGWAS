import math
import numpy as np

import Curve


class BiExponential(Curve.BaseCurve):
    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(count=4, name=["a1", "r1", "a2", "r2"],
                               formula="y = a1*exp(-r1*t) + a2*exp(-r2*t)")

    def show(self):
        super().show(obj=self)

    def est_init_param(self, pheY, pheX, pheT, **options):
        r = -1 * np.nanmean((np.log(pheY.iloc[:, 1]) - np.log(pheY.iloc[:, 0])) / (pheT.iloc[:, 1] - pheT.iloc[:, 0]))
        a_double = np.nanmean(pheY.iloc[:, 0] / np.exp(-1 * r * pheT.iloc[:, 0]))

        return [a_double / 2, r, a_double / 2, r]

    def get_simu_param(self, times):
        return np.array(
            [[19.9824, 0.4699, 8.7768, 1.4699], [17.9824, 0.0699, 9.7768, 1.0699], [15.9507, 0.1836, 10.5737, 1.8836]])

    def check_param(self, par, times, *options):
        return True

    def get_gradient(self, par, times, *options):
        d_a1 = np.exp(-par[2] * times)
        d_r1 = par[1] * np.exp(-par[2] * times) * (-times)
        d_a2 = np.exp(-par[4] * times)
        d_r2 = par[3] * np.exp(-par[4] * times) * (-times)

        return [d_a1, d_r1, d_a2, d_r2]

    def get_curve_formula(self, par, times, *options):
        return par[1] * math.np.exp(-par[2] * times) + par[3] * math.np.exp(-par[4] * times)
