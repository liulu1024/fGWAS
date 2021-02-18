import math

import Covariance.CovarianceTool
import numpy as np
import Covariance


class CS(Covariance.BaseCovariance):
    def __init__(self, covarType, description):
        super().__init__(covarType, description)

    def get_matrix(self, par, times, **kwargs):
        rho = par[0]
        s2 = abs(par[1])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size

        sigma = abs(s2) * rho ^ abs(np.ones((n, n)) - np.identity(n))
        return sigma

    def get_gradient(self, par, times, **kwargs):
        rho = par[0]
        s2 = abs(par[1])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size

        d_rho = abs(s2) * rho ^ np.zeros((n, n))
        d_s2 = rho ^ abs(np.ones((n, n)) - np.identity(n))

        return [d_rho, d_s2]

    def get_param_info(self, par, times, **options):
        n = len(times)
        return Covariance.ParamInfo(count=2, name=["rho","sigma2"])

    def check_param(self, par, times, **kwargs):
        if par[0] > 1 or par[0] < 0:
            return False
        else:
            return True

    def get_simu_param(self, times, **options):

        return [0.75, 1.2]

    def est_init_param(self, pheY, pheX, pheT, **options):
        sel =np.isnan(pheY[:, 0]) | np.isnan(pheY[:, 1])
        rho =np.corrcoef(pheY[not sel, 0], pheY[not sel, 1])
        s2 =np.nanstd(pheY, ddof=1)^ 2

        return [rho, s2 * np.random.uniform(0.9, 1.1,1)]
