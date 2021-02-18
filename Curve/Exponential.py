import math
import numpy as np

import Curve
import Curve.CurveTool as tool


class Exponential(Curve.BaseCurve):
    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(count=2, name=["a", "r"], formula="y = a*exp(r*t)")

    def show(self):
        super().show(obj=self)

    def est_init_param(self, pheY, pheX, pheT, **options):
        r = np.nanmean((math.log(pheY[:, 1]) - math.log(pheY[:, 0])) / (pheT[:, 1] - pheT[:, 0]))
        w0 = np.nanmean(pheY[:, 0] / math.exp(r * pheT[:, 0]))

        return [w0, r]

    def get_simu_param(self, times):
        return np.array([[2, 0.0128], [1.8, 0.02], [1.6, 0.024]])

    def check_param(self, par, times, *options):
        return True

    def get_gradient(self, par, times, *options):
        d_a = math.exp(par[2] * times)
        d_r = par[1] * math.exp(par[2] * times) * times

        return [d_a, d_r]

    def get_curve_formula(self, par, times, **options):
        return par[[1]] * math.exp(par[[2]] * times)
