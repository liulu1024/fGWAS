import math

import Covariance.CovarianceTool
import numpy as np
import Covariance


class AR1(Covariance.BaseCovariance):
    def __init__(self, covarType, description):
        super().__init__(covarType, description)

    def get_matrix(self, par, times, *options):
        rho = par[0]
        s2 = abs(par[1])
        n = times.shape[1] if isinstance(times, np.ndarray) else times.columns.size
        matrix = np.reshape(np.tile(range(0, n), n), (n, n))
        matrix1 = np.reshape(np.tile(range(0, n), n), (n, n), order='F')
        sigma = abs(s2) * rho ** abs(matrix - matrix1)
        return sigma

    def get_gradient(self, par, times, options=list()):
        rho =par[0]
        s2 =abs(par[1])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        matrix = np.reshape(np.tile(range(0, n), n), (n, n))
        matrix1 = np.reshape(np.tile(range(0, n), n), (n, n), order='F')
        rho_times =abs(matrix -matrix1)

        d_rho = abs(s2) * rho_times * rho ** (rho_times - 1)
        d_s2 = rho ** rho_times
        return [d_rho, d_s2]

    def get_param_info(self, times, **options):
        return Covariance.ParamInfo(count=2, name=["rho", "sigma2"])

    def check_param(self, par, times, **kwargs):
        if par[0] > 1 or  par[0] < 0:
            return  False
        else:
            return True

    def get_simu_param(self, times, **options):
       return [0.75, 1.2]

    def est_init_param(self, pheY, pheX, pheT, **options):
        def func(i):
            sel = np.isnan(pheY[:, i]) | np.isnan(pheY[:, i + 1])
            return np.corrcoef(pheY[~sel, i], pheY[~sel, i + 1])
        rho = np.nanmean([func(i) for i in range(0,pheY.shape[1]-1)])

        s2 = np.nanstd(pheY, ddof=1) ** 2

        return np.append(rho, s2 * np.random.uniform(0.9, 1.1,1))
