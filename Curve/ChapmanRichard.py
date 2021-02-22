import math

import numpy as np

import Curve


class ChapmanRichard(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(count=3, name=["a", "b", "r"], formula="y = a*(1-exp(-rt))^b")

    def show(self):
        super().show(obj=self)

    def get_simu_param(self, times):
        return np.array([[21.98, 0.47, 9.78], [19.98, 0.47, 8.77], [15.95, 0.48, 7.58]])

    def check_param(self, par, times, **options):
        return True

    def get_gradient(self, par, times, **options):
        d_a = (1 - np.exp(-par[3] * times)) ^ par[2]
        d_b = par[1] * (1 - np.exp(-par[3] * times)) ^ par[2] * np.log(1 - np.exp(-par[3] * times))
        d_r = par[1] * (1 - np.exp(-par[3] * times)) ^ (par[2] - 1) * par[2] * np.exp(-par[3] * times) * (
                    -1 * times)

        return [d_a, d_b, d_r]

    def get_curve_formula(self, par, times, **options):
        return par[1] * (1 - np.exp(-par[3] * times)) ^ par[2]

    def est_init_param(self, pheY, pheX, pheT, **options):
        a = np.nanmax(pheY)
        b = 1
        r = np.nanmean(np.log(1 - pheY / a) / pheT / (-1), trim=0.2)
        if math.isinf(r): r = 0

        return [a, b, r]
