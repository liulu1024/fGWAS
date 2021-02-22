import numpy as np

import Curve
import Curve.Logistic as logistic

'''----------------------------------------------------------
 ABRK:

    y = a*(1+b*exp(-r*t))^(1/(1-k))

 Reference:<no>
----------------------------------------------------------'''


class ABRK(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_curve_formula(self, par, times, **args):
        return par[1] * (1 + par[2] * np.exp(-par[3] * times)) ^ (1 / (1 - par[4]))

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo( 4,["a", "b", "r", "k"], "y = a*(1+b*exp(-r*t))^(1/(1-k))")

    def get_gradient(self, par, times, options=list()):
        da = (1 + par[2] * np.exp(-par[3] * times)) ^ (1 / (1 - par[4]))
        db = par[1] * (1 / (1 - par[4])) * (1 + par[2] * np.exp(-par[3] * times)) ^ (1 / (1 - par[4]) - 1) * np.exp(
            -par[3] * times)
        dr = par[1] * (1 / (1 - par[4])) * (1 + par[2] * np.exp(-par[3] * times)) ^ (1 / (1 - par[4]) - 1) * par[
            2] * np.exp(-par[3] * times) * (-1 * times)
        dk = par[1] * (1 + par[2] *np.exp(-par[3] * times)) ^ (1 / (1 - par[4])) * np.log(
            (1 + par[2] * np.exp(-par[3] * times))) * (1 - par[4]) ^ (-2)

        return list(da, db, dr, dk)

    def check_param(self, par, times, options=list()):
        return True

    def show(self):
        super().show(obj=self)

    def get_simu_param(self):
        return np.array([[18.18, 9.98, 0.99, 2.6], [17.08, 9.78, 0.97, 2.5], [15.95, 9.88, 0.9, 2.48]])

    def est_init_param(self, pheY, pheX, pheT, options=list()):
        par = logistic.Logistic().est_init_param(pheX, pheT, options)
        par_k = 2
        return np.vstack(par, par_k) * np.random.normal(0.95, 1.05, 4)
