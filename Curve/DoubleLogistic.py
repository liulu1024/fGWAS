from math import exp

import numpy as np

import Curve
import Curve.Logistic as logistic

'''-----------------------------------------------------------
  curveType Double Logistic

    y = a1/(1+b1*exp(-r1*t)) +  a2/(1+b2*exp(-r2*t))

-----------------------------------------------------------'''


def get_curve_formula(par, times, **options):
    return par[1] / (1 + par[2] * exp(-par[3] * times)) + par[4] / (1 + par[5] * exp(-par[6] * times))


def get_gradient(par, times, options=list()):
    da1 = 1 / (1 + par[2] * exp(-1 * par[3] * times))
    db1 = (-1) * par[1] / ((1 + par[2] * exp(-1 * par[3] * times)) ^ 2) * exp(-1 * par[3] * times)
    dr1 = (-1) * par[1] / ((1 + par[2] * exp(-1 * par[3] * times)) ^ 2) * par[2] * exp(-1 * par[3] * times) * (
            -1 * times)

    da2 = 1 / (1 + par[5] * exp(-1 * par[6] * times))
    db2 = (-1) * par[4] / ((1 + par[5] * exp(-1 * par[6] * times)) ^ 2) * exp(-1 * par[6] * times)
    dr2 = (-1) * par[4] / ((1 + par[5] * exp(-1 * par[6] * times)) ^ 2) * par[5] * exp(-1 * par[6] * times) * (
            -1 * times)

    return list(da1, db1, dr1, da2, db2, dr2)


def get_simu_param():
    return np.array([[18.18, 9.98, 0.99], [17.08, 9.78, 0.97], [15.95, 9.88, 0.98]])


def check_param(par, times, options=list()):
    return True


def est_init_param(pheY, pheX, pheT, options=list()):
    par = logistic.est_init_param(pheX, pheT, options)
    par[1] = par[1] / 2
    return np.vstack(par, par)


def get_param_info(self,pheT,**kwargs):
    return Curve.ParamInfo( 6, ["a1", "b1", "r1", "a2", "b2", "r2"],  "y = a1/(1+b1*exp(-r1*t)) +  "
                                                                                  "a2/(1+b2*exp(-r2*t))")


class DoubleLogistic(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def show(self):
        super().show(obj=self)
