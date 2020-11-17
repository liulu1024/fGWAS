import math
from math import exp
from math import log

import numpy as np

import CV.BaseCurve

'''-----------------------------------------------------------
   Logistic curve
    y = a/(1+b*exp(-r*t))
-----------------------------------------------------------'''


def get_simu_param():
    return np.array([[18.18, 9.98, 0.99], [17.08, 9.78, 0.97], [15.95, 9.88, 0.98]])


def check_param(par, times, options=list()):
    return True


def get_gradient(par, times, options=list()):
    da = 1 / (1 + par[2] * exp(-1 * par[3] * times))
    db = (-1) * par[1] / ((1 + par[2] * exp(-1 * par[3] * times)) ^ 2) * exp(-1 * par[3] * times)
    dr = (-1) * par[1] / ((1 + par[2] * exp(-1 * par[3] * times)) ^ 2) * par[2] * exp(-1 * par[3] * times) * (
            -1 * times)
    return list(da, db, dr)


def get_param_info():
    return {'count': 3, 'names': ['a', 'b', 'r'], 'formula': "y = a/(1+b*exp(-r*t))"}


def get_curve_formula(par, times, options):
    return par[0] / (1 + par[1] * math.exp(-par[3] * times))


class Logistic(CV.BaseCurve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description


def est_init_param(pheY, pheX, pheT, options=list()):
    mc =CV.BaseCurve.get_mean_vector(pheY, pheT);
    mc_t = mc['t'][mc['y'] > 0]
    mc_y = mc['y'][mc['y'] > 0]
    m = len(mc['t'])
    if m == 0:
        mc = CV.BaseCurve.get_mean_vector(pheY, pheT);
        mc['y']= mc['y']- min(mc.y) * 1.01;
    par = []
    # ls_i useless
    ls_i = ls_max = float("inf")
    a_rate = np.mean(mc.y[-1] / mc.y[-m]);
    if a_rate == 1:
        a_rate = 0.99
    for i in range(1, 11):
        par_a = mc.y[m] * a_rate ^ i
        try:
            par_r = log(par_a / mc.y[1] - 1) - log(par_a / mc.y[m] - 1) / (mc['t'][1] - mc['t'][m])
        except(Exception):
            continue
        if np.isinf(par):
            continue
        par_b = (par_a / mc.y[m] - 1) / exp(-par_r * mc['t'][m])
        # calculate loss to get more accurate value
        y_ls = np.sum(abs(mc['y']- par_a / (1 + par_b * exp(-par_r * mc['t']))) ^ 2)
        if y_ls < ls_max:
            ls_i = i
            ls_max = y_ls
            par = par.append(par_a).append(par_b).append(par_r)
    return par * np.random.normal(0.95, 1.05, len(par))
