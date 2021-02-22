import numpy as np

import Curve
import Curve.CurveTool as tool

'''-----------------------------------------------------------
   Pharmacology Curve

    y = E0 + Emax*t/(E50+t)

-----------------------------------------------------------'''


class Pharmacology(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_curve_formula(self, par, times, **options):
        return par[1] + par[3] * times / (par[2] + times)

    def show(self):
        super().show(obj=self)

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(3, ["E0", "E50", "Emax"], "y = E0 + Emax*t/(E50+t)")

    def get_gradient(self, par, times, options=list()):
        dE0 = times / times
        dE50 = (-1) * par[3] * times / (par[2] + times) ^ 2
        dEmax = times / (par[2] + times)

        return list(dE0, dE50, dEmax)

    def check_param(self, par, times, options=list()):
        return True

    def get_simu_param(self):
        return np.array([[0.9824, 15.909, 20.7768], [8.9824, 16.098, 20.7768], [6.9507, 12.090, 18.5737]])

    def est_init_param(self, pheY, pheX, pheT, options=list()):
        mc = tool.get_mean_vector(pheY, pheT)
        m = len(mc['t'])

        par = list()
        ls_max = float('Inf')
        a_rate = mc['y'][m] / mc['y'][m - 1]
        for i in range(1, 11):
            par_Emax = mc['y'][m] * a_rate ^ i
            par_E50 = mc['y'][np.where(mc['t'] == np.median(mc['t']))]
            par_E0 = mc['y'][m] - par_Emax * mc['t'][m] / (par_E50 + mc['t'][m]);

            y_ls = np.sum(abs(mc['y'] - par_E0 - par_Emax * mc['t'] / (par.E50 + mc['t'])) ^ 2, axis=1)
            if y_ls < ls_max:
                ls_max = y_ls
                par = (par_E0, par_E50, par_Emax)

        return par * np.random.uniform(0.95, 1.05, 3)
