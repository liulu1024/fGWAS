import math

import numpy as np

import Curve

class Power(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(count=2, name=["a", "b"], formula="y = a*t^b")

    def show(self):
        super().show(obj=self)

    def est_init_param(self, pheY, pheX, pheT, **options):
        b = np.nanmean((math.log(pheY[:, pheY.columnx.size - 1]) - math.log(pheY[:, 0])) / (
                    math.log(pheT[:, pheT.columns.size - 1]) - math.log(pheT[:, :])))
        a = math.exp(
            np.nanmean(math.log(pheY[:, pheY.columnx.size - 1]) - b * math.log(pheT[:, pheY.columnx.size - 1])))
        return [a, b]

    def get_simu_param(self, times):
        return np.array([[11.049, 1.151], [9.049, 1.251], [7.148, 1.359]])

    def check_param(self, par, times, **options):
        return True

    def get_gradient(self, par, times, **options):
        return [times ^ par[2], par[1] * (times ^ par[2]) * math.log(times)]

    def get_curve_formula(self, par, times, **options):
        return par[1] * (times ^ par[2])
