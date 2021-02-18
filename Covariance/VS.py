import math

import Covariance.CovarianceTool
import numpy as np
import Covariance


class VS(Covariance.BaseCovariance):
    def __init__(self, covarType, description):
        super().__init__(covarType, description)

    def get_matrix(self, par, times, **kwargs):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        s2 = abs(par)

        sigma = np.diag(s2)
        return sigma

    def get_gradient(self, par, times, **kwargs):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        s2 = abs(par)

        d_ret = list()
        for i in range(0, len(s2)):
            A = np.zeros((n, n))
            A[i, i] = 1
            d_ret[i] = A

        return d_ret

    def get_param_info(self, par, times, **options):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return Covariance.ParamInfo(count=1 + n, name=["rho"] + ['sigma2_'.format(i) for i in range(0, n)])

    def check_param(self, par, times, **kwargs):
        return  True

    def get_simu_param(self, times, **options):

        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return np.arange(1.2,1.0+n*0.2,0.2)

    def est_init_param(self, pheY, pheX, pheT, **options):
        return  np.std(pheY,ddof=1)^2
