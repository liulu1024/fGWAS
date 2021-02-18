import math

import Covariance.CovarianceTool
import numpy as np
import Covariance


class FA1(Covariance.BaseCovariance):
    def __init__(self, covarType, description):
        super().__init__(covarType, description)

    def get_matrix(self, par, times, **kwargs):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        d =abs(par[0])
        s2 =abs(par[0,len(par)-1])

        sigma = d * np.identity(n) + math.sqrt(np.reshape(s2, (n,n), order='F') *
                                      np.reshape(s2, (n,n), order='C'))
        return sigma

    def get_gradient(self, par, times, **kwargs):

        rho = par[0]
        s2 = abs(par[0:len(par)-1])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        d_ret =list()
        d_ret[0] = d_d = np.identity(n)

        d_ret = list()
        d_ret[0] = d_d = np.identity(n)
        A =np.c_(np.tile(0.5 * math.sqrt(s2))[0:n],np.zeros((n, n - 1)))
        for  i in range(0,len(s2)):
            d_ret[i+1] = (A + A.T) / math.sqrt(s2[i])
            A = np.c_[np.zeros((A.shape[0], 1)),A[:,0:-1]]

        return d_ret

    def get_param_info(self, par, times, **options):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return Covariance.ParamInfo(count=1 + n, name=["d"]+['sigma2_'.format(i) for i in range(0, n)])

    def check_param(self, par, times, **kwargs):
            return True

    def get_simu_param(self, times, **options):

        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return [0.5,np.arange(1.2, 1.0+ n*0.2, 0.2)]

    def est_init_param(self, pheY, pheX, pheT, **options):
        return [np.nanstd(pheY,ddof=1) ^ 2 / 100, np.std(pheY,axis=0,ddof=1) ^ 2 * np.random.uniform(0.8, 1.2,1)]
