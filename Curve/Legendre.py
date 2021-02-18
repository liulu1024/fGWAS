import math

import numpy as np

import Curve

class Legendre2(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(count=3, name=["u0","u1","u2"], formula="y = u0 + u1*t + u2*1/2*(3*t^2-1)")

    def show(self):
        super().show(obj=self)

    def get_simu_param(self, times):
        # QQ2 = [simu_u0=11.049, simu_u1=1.551, simu_u2=-8.019, simu_u3=3.151, simu_u4=0.652, simu_u5=-0.597,
        #         simu_u6=0.821]
        # Qq1 = [simu_u0=9.049, simu_u1=1.151, simu_u2=-6.019, simu_u3=2.651, simu_u4=0.652, simu_u5=-0.797,
        #         simu_u6=0.621]
        # qq0 = [simu_u0=7.148, simu_u1=1.379, simu_u2=-4.489, simu_u3=2.004, simu_u4=0.662, simu_u5=-0.836,
        #         simu_u6=0.432)
        return np.array([[11.049,1.551,-8.019], [9.049,1.151,-6.019], [7.148,1.379,-4.489]])

    def check_param(self, par, times, **options):
        return True

    def get_gradient(self, par, times, **options):
        ti = -1 + 2 * (times - options["min_time"]) / (options["max_time"] - options["min_time"])
        return [np.ones(ti.shape[0],ti.shape[1]), ti, 0.5 * (3 * ti * ti - 1)]

    def get_curve_formula(self, par, times, **options):
        ti = -1 + 2 * (times - options["min_time"]) / (options["max_time"] - options["min_time"])
        return par[1] + ti * par[2] + 0.5 * (3 * ti * ti - 1) * par[3]


def est_init_param(self, pheY, pheX, pheT, **options):
        ti = -1 + 2 * (pheT - options["min_time"]) / (options["max_time"] - options["min_time"])
        y1 = pheY[:, 1]
        y2 = pheY[:, 2]
        y3 = pheY[:, 3]

        u2 =  np.nanmean(((y1 - y2) / (ti[:, 0]-ti[:, 1]) - (y1 - y3) / (ti[:, 0]-ti[:, 2])) / (
                    (ti[:, 0]+ti[:, 1]) - (ti[:, 0]+ti[:, 2])) / 1.5)
        u1 =  np.nanmean((y1 - y2) / (ti[:, 0]-ti[:, 1]) - u2 * 1.5 * (ti[:, 0]+ti[:, 1]))
        u0 =  np.nanmean(y3 + u2 * 0.5 - u1 * ti[:, 2] - u2 * 1.5 * ti[:, 2] ^ 2)

        return [u0, u1, u2]

'''
===========================================================================
'''
class Legendre3(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(count=4, name=["u0", "u1", "u2", "u3"],
                               formula="y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) ")

    def show(self):
        super().show(obj=self)

    def est_init_param(self, pheY, pheX, pheT, **options):
        u3 =est_init_param(object, pheY, pheX, pheT, options)
        return [u3, 0.0001]

    def get_simu_param(self, times):
        # QQ2 = [simu_u0=11.049, simu_u1=1.551, simu_u2=-8.019, simu_u3=3.151, simu_u4=0.652, simu_u5=-0.597,
        #         simu_u6=0.821]
        # Qq1 = [simu_u0=9.049, simu_u1=1.151, simu_u2=-6.019, simu_u3=2.651, simu_u4=0.652, simu_u5=-0.797,
        #         simu_u6=0.621]
        # qq0 = [simu_u0=7.148, simu_u1=1.379, simu_u2=-4.489, simu_u3=2.004, simu_u4=0.662, simu_u5=-0.836,
        #         simu_u6=0.432)
        return np.array([[11.049, 1.551, -8.019, 3.151], [9.049, 1.151, -6.019, 2.651], [7.148, 1.379, -4.489, 2.004]])

    def check_param(self, par, times, **options):
        return True

    def get_gradient(self, par, times, **options):
        ti = -1 + 2 * (times - options["min_time"]) / (options["max_time"] - options["min_time"])
        return [np.ones(ti.shape[0], ti.shape[1]), ti, 0.5 * (3 * ti * ti - 1), 0.5 * (5 * ti ^ 3 - 3 * ti)]

    def get_curve_formula(self, par, times, **options):
        ti = -1 + 2 * (times - options["min_time"]) / (options["max_time"] - options["min_time"])
        return (par[1] + ti * par[2] + 0.5 * (3 * ti * ti - 1) * par[3] + 0.5 * (5 * ti ^ 3 - 3 * ti) * par[4])

'''========================================================================='''
class Legendre4(Curve.BaseCurve):

    def __init__(self, curveType, description):
        super().__init__(curveType, description)

    def get_curve_type(self):
        return self.curveType

    def get_description(self):
        return self.description

    def get_param_info(self,pheT,**kwargs):
        return Curve.ParamInfo(count=5, name=["u0","u1","u2","u3", "u4"], formula="y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) + u4*1/8*(35*ti^4-30*ti^2+3)")

    def show(self):
        super().show(obj=self)

    def est_init_param(self, pheY, pheX, pheT, **options):
        u3 =est_init_param(object, pheY, pheX, pheT, options)
        return [u3, 0.0001, 0.0001]

    def get_simu_param(self, times):
        # QQ2 = [simu_u0=11.049, simu_u1=1.551, simu_u2=-8.019, simu_u3=3.151, simu_u4=0.652, simu_u5=-0.597,
        #         simu_u6=0.821]
        # Qq1 = [simu_u0=9.049, simu_u1=1.151, simu_u2=-6.019, simu_u3=2.651, simu_u4=0.652, simu_u5=-0.797,
        #         simu_u6=0.621]
        # qq0 = [simu_u0=7.148, simu_u1=1.379, simu_u2=-4.489, simu_u3=2.004, simu_u4=0.662, simu_u5=-0.836,
        #         simu_u6=0.432)
        return np.array([[11.049, 1.551, -8.019, 3.151,0.652], [9.049, 1.151, -6.019, 2.651,0.652], [7.148, 1.379, -4.489, 2.004,0.662]])

    def check_param(self, par, times, **options):
        return True

    def get_gradient(self, par, times, **options):
        ti = -1 + 2 * (times - options["min_time"]) / (options["max_time"] - options["min_time"])
        return [np.ones(ti.shape[0], ti.shape[1]), ti, 0.5 * (3 * ti * ti - 1), 0.5 * (5 * ti ^ 3 - 3 * ti), 0.125*(35*ti^4-30*ti^2+3)]

    def get_curve_formula(self, par, times, **options):
        ti = -1 + 2 * (times - options["min_time"]) / (options["max_time"] - options["min_time"])
        return (par[1] + ti * par[2] + 0.5 * (3 * ti * ti - 1) * par[3] + 0.5 * (5 * ti ^ 3 - 3 * ti) * par[4]+  0.125*(35*ti^4-30*ti^2+3)* par[5])