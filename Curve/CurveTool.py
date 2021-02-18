import numpy as np
import torch
import Curve.Logistic
import Curve.BiExponential
import Curve.ABRK
import Curve.ChapmanRichard
import Curve.DoubleLogistic
import Curve.Exponential
import Curve.Legendre
import Curve.Pharmacology
import Curve.Power

import Curve

curves = []


def get_param_info(obj):
    return obj.get_param_info()


def get_gradient(obj, *args):
    return obj.get_gradient(*args)


def check_param(obj, *args):
    return obj.check_param(*args)


def get_simu_param(obj, *args):
    return obj.get_simu_param(times=args)


def est_init_param(obj, *args):
    return obj.est_init_param(*args)

    '''
       pheT: 1:time.points np.array 便于标量向量运算
       pheY: numpy
       '''


def get_mean_vector(pheY, pheT):
    # mean t
    t_count = len(np.unique(pheT))
    if t_count >= 20:
        t_min = min(pheT)
        t_max = max(pheT)
        pheT = np.round((pheT - t_min) / (t_max - t_min) * 20) / 20 * (t_max - t_min) + t_min
    t_all = np.sort(np.unique(pheT))
    # mean y
    y_all = list()
    for t in t_all:
        y_all.append(np.nanmean(np.array(pheY)[np.array(pheT) == t]))
    return {'t': t_all, 'y': y_all}


def get_curve(type):
    curve = None
    if isinstance(type, str):
        for obj in get_all_curve():
            if obj.curveType.upper() == type.upper():
                curve = obj
                curve.show()
                return curve

    if isinstance(type, int):
        curve = get_all_curve()[type]
        curve.show()
        return curve


'''
all curve 
'''


def get_all_curve():
    if len(curves) == 0:
        curves.append(Curve.Logistic.Logistic(curveType="Logistic", description="logistic curve"))

        curves.append(Curve.BiExponential.BiExponential(curveType="Bi-Logistic", description="Double logistic curve"))
        curves.append(Curve.ABRK.ABRK(curveType="ABRK", description="ABRK model"))

        curves.append(Curve.Pharmacology.Pharmacology(curveType="Pharmacology", description="Pharmacology curve"))
        curves.append(Curve.Exponential.Exponential(curveType="Exponential", description="Exponential curve"))
        curves.append(Curve.BiExponential.BiExponential(curveType="Bi-Exponential", description="Bi-exponential curve"))

        curves.append(Curve.Power.Power(curveType="Power", description="power curve"))

        curves.append(Curve.Legendre.Legendre2(curveType="Legendre2", description="Legendre Polynomial(2nd-order)"))
        curves.append(Curve.Legendre.Legendre3(curveType="Legendre3", description="Legendre Polynomial(3rd-order)"))
        curves.append(Curve.Legendre.Legendre4(curveType="Legendre4", description="Legendre Polynomial(4th-order)"))

    return curves


def get_curve_count():
    return len(get_all_curve())
