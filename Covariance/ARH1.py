import math

import Covariance.CovarianceTool
import numpy as np
import Covariance


class ARH1(Covariance.BaseCovariance):
    def __init__(self, covarType, description):
        super().__init__(covarType, description)

    def get_matrix(self, par, times, *options):
        rho = par[1]
        s2 = abs(par[-1])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        matrix = np.reshape(np.repeat(range(0, n)), (n, n), order='F')
        matrix2 = np.reshape(np.repeat(range(0, n)), (n, n), order='C')
        sigma = rho ^ abs(matrix - matrix2) * math.sqrt(
            np.reshape(s2, (n, n), order='F') * np.reshape(s2, (n, n), order='C'))

        return sigma

    def get_gradient(self, par, times, options=list()):
        rho = par[0]
        s2 = abs(par[0:len(par)-1])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size

        matrix = np.reshape(np.repeat(range(0, n)), (n, n), order='F')
        matrix2 = np.reshape(np.repeat(range(0, n)), (n, n), order='C')
        rho_times = abs(matrix - matrix2)
        sn_mat = math.sqrt(np.reshape(s2, (n, n), order='F') * np.reshape(s2, (n, n), order='C'))

        d_ret = list()
        d_ret[1] = d_rho = rho_times * rho ^ (rho_times - 1) * sn_mat

        A = np.c_[np.tile(0.5 * math.sqrt(s2))[0:n], np.zeros((n, n - 1))]
        for i in range(0, len(s2)):
            d_ret[i + 1] = rho ^ rho_times * (A + A.T) / math.sqrt(s2[i])
            A = np.c_[np.zeros((A.shape[0], 1)),A[:,0:-1]]

        return d_ret

    def get_param_info(self, par, times, **options):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return Covariance.ParamInfo(count=1 + n, name=["rho"] + ['sigma2_'.format(i) for i in range(0, n)])

    def check_param(self, par, times, **kwargs):
        if par[0] > 1 or par[0] < 0:
            return False
        else:
            return True

    def get_simu_param(self, times, **options):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return [0.75, np.arange(1.2, 1 + n * 0.2, 0.2)]

    def est_init_param(self, pheY, pheX, pheT, **options):
        n = pheT.shape[0] if isinstance(pheT, np.ndarray) else pheT.columns.size

        def func(i):
            sel = np.isnan(pheY[:, i-1]) | np.isnan(pheY[:, i])
            return np.corrcoef(pheY[sel==False, i-1], pheY[sel==False, i])

        rho = np.nanmean([func(i) for i in range(0, pheY.shape[1])])

        s2 = np.std(pheY, axis=0, ddof=1) ** 2

        return [rho, s2 * np.random.uniform(0.9, 1.1, n)]
