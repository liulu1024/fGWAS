import numpy as np
import torch

import CV.Logistic


class BaseCurve(object):
    curveType = 'Base'
    description = 'BaseCurve'

    def __init__(self, curveType, description):
        self.curveType = curveType
        self.description = description

    def show(obj):
        print("     Class :", obj.__class__, "\n")
        print("Curve Type :", obj.curveType, "\n")

        info = obj.get_param_info()
        print("Parameters :", info['names'], "\n")
        print("   Formula :", info['formula'], "\n")


curves = []


def get_curve_formula(obj, *args):
    return obj.get_curve_formula(*args)


def get_param_info(obj):
    return obj.get_param_info()


def get_gradient(obj, *args):
    return obj.get_gradient(*args)


def check_param(obj, *args):
    return obj.check_param(*args)


def get_simu_param(obj, *args):
    return obj.get_simu_param(*args)


def est_init_param(obj, *args):
    return obj.est_init_param(*args)

    '''
       pheT: 1:time.points np.array 便于标量向量运算
       pheY: numpy
       '''


def get_mean_vector(self, pheY, pheT):
    # mean t
    t_count = len(np.unique(pheT))
    if t_count >= 20:
        t_min = min(pheT)
        t_max = max(pheT)
    pheT = np.round((pheT - t_min) / (t_max - t_min) * 20) / 20 * (t_max - t_min) + t_min
    t_all = np.sort(np.unique(pheT))
    # mean y
    y_all = []
    pheT = torch.tensor(pheT)
    for t in t_all:
        y_all = y_all.append(np.mean(pheY[pheT == t]))
    return {'t': t_all, 'y': y_all}


def get_curve(type):
    curve = None
    if isinstance(type, str):
        for obj in get_all_curve():
            if obj.curveType.upper() == type.upper():
                curve = obj

    if isinstance(type, int):
        curve = get_all_curve()[type]

    curve.show(curve)
    return curve


'''
all curve 
'''


def get_all_curve():
    if len(curves) == 0:
        curves.append(CV.Logistic.Logistic(curveType="Logistic", description="logistic curve"))
        # self.curves.append(BaseCurve(curveType="Bi-Logistic", description="Double logistic curve"))
        # self.curves.append(BaseCurve(curveType="ABRK", description="ABRK model"))
        #
        # self.curves.append(BaseCurve(curveType="Pharmacology", description="Pharmacology curve"))
        # self.curves.append(BaseCurve(curveType="Exponential", description="Exponential curve"))
        # self.curves.append(BaseCurve(curveType="Bi-Exponential", description="Bi-exponential curve"))
        #
        # self.curves.append(BaseCurve(curveType="Power", description="power curve"))
        #
        # self.curves.append(BaseCurve(curveType="Legendre2", description="Legendre Polynomial(2nd-order)"))
        # self.curves.append(BaseCurve(curveType="Legendre3", description="Legendre Polynomial(3rd-order)"))
        # self.curves.append(BaseCurve(curveType="Legendre4", description="Legendre Polynomial(4th-order)"))

    return curves


def get_curve_count():
    return len(get_all_curve())