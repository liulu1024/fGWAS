from math import exp
from math import log

import numpy as np

import Curve
import Curve.Logistic as logistic

'''----------------------------------------------------------
 ABRK:

    y = a*(1+b*exp(-r*t))^(1/(1-k))

 Reference:<no>
----------------------------------------------------------'''


class Logistic(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_curve_formula(self, par, times, options):
        return par[1] * (1 + par[2] * exp(-par[3] * times)) ^ (1 / (1 - par[4]))

    def get_param_info(self):
        return {'count': 4, 'names': ["a", "b", "r", "k"], 'formula': "y = a*(1+b*exp(-r*t))^(1/(1-k))"}

    def get_gradient(self, par, times, options=list()):
        da = (1 + par[2] * exp(-par[3] * times)) ^ (1 / (1 - par[4]))
        db = par[1] * (1 / (1 - par[4])) * (1 + par[2] * exp(-par[3] * times)) ^ (1 / (1 - par[4]) - 1) * exp(
            -par[3] * times)
        dr = par[1] * (1 / (1 - par[4])) * (1 + par[2] * exp(-par[3] * times)) ^ (1 / (1 - par[4]) - 1) * par[
            2] * exp(-par[3] * times) * (-1 * times)
        dk = par[1] * (1 + par[2] * exp(-par[3] * times)) ^ (1 / (1 - par[4])) * log(
            (1 + par[2] * exp(-par[3] * times))) * (1 - par[4]) ^ (-2)

        return list(da, db, dr, dk)

    def check_param(self, par, times, options=list()):
        return True

    def get_simu_param(self):
        return np.array([[18.18, 9.98, 0.99, 2.6], [17.08, 9.78, 0.97, 2.5], [15.95, 9.88, 0.9, 2.48]])

    def est_init_param(self, pheY, pheX, pheT, options=list()):
        par = logistic.est_init_param(pheX, pheT, options)
        par_k = 2
        return np.vstack(par, par_k) * np.random.normal(0.95, 1.05, 4)
