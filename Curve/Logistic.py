import math
from math import exp
from math import log

import numpy as np

import Curve.CurveTool as tool
import  Curve

'''-----------------------------------------------------------
   Logistic curve
    y = a/(1+b*exp(-r*t))
-----------------------------------------------------------'''
class Logistic(Curve.BaseCurve):

     def __init__(self, curveType, description):
        super().__init__(curveType, description)

     def get_curve_type(self):
        return self.curveType

     def get_description(self):
        return self.description

     def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo ( 3,  ['a', 'b', 'r'], "y = a/(1+b*exp(-r*t)")

     def show(self):
         super().show(obj=self)


     def est_init_param(self,pheY, pheX, pheT, **options):
        mc =tool.get_mean_vector(pheY, pheT)
        mc['t'] = np.array(mc['t'])[np.array(mc['y']) > 0]
        mc['y'] = np.array(mc['y'])[np.array(mc['y']) > 0]
        m = len(mc['t'])
        if m == 0:
            mc = tool.get_mean_vector(pheY, pheT)
            mc['y']= mc['y']- np.min(mc['y']) * 1.01
        par = list()
        # ls_i useless
        ls_i = ls_max = float("inf")
        a_rate = np.mean(np.delete(mc['y'],0) / np.delete(mc['y'],m-1))
        if a_rate == 1:
            a_rate = 0.99
        for i in range(0, 10):
            par_a = mc['y'][m-1] * a_rate ** (i+1)
            try:
                par_r = np.log(par_a / mc['y'][0] - 1) - np.log(par_a / mc['y'][m-1] - 1) / (mc['t'][0] - mc['t'][m-1])
            except Exception:
                continue
            if np.isinf(par_r):
                continue
            par_b = (par_a / mc['y'][m-1] - 1) / np.exp(-par_r * mc['t'][m-1])
            # calculate loss to get more accurate value
            y_ls = np.nansum(abs(mc['y']- par_a / (1 + par_b * np.exp(-par_r * mc['t']))) ** 2)
            if y_ls < ls_max:
                ls_i = i
                ls_max = y_ls
                par.append(par_a)
                par.append(par_b)
                par.append(par_r)
        return par * np.random.normal(0.95, 1.05, len(par))

     def get_simu_param(self, times):
            return np.array([[18.18, 9.98, 0.99], [17.08, 9.78, 0.97], [15.95, 9.88, 0.98]])

     def check_param(self,par, times, **options):
            return True

     def get_gradient(self,par, times, **options):
            da = 1 / (1 + par[2] * exp(-1 * par[3] * times))
            db = (-1) * par[1] / ((1 + par[2] * exp(-1 * par[3] * times)) ** 2) * exp(-1 * par[3] * times)
            dr = (-1) * par[1] / ((1 + par[2] * exp(-1 * par[3] * times)) ** 2) * par[2] * exp(-1 * par[3] * times) * (
                    -1 * times)
            return [da, db, dr]

     def get_curve_formula(self, par, times, **options):
           return par[0] / (1 + par[1] * np.exp(-par[2] * times))


