import math

import Covariance.CovarianceTool
import numpy as np
import Covariance


class CSH(Covariance.BaseCovariance):
    def __init__(self, covarType, description):
        super().__init__(covarType, description)

    def get_matrix(self, par, times, **kwargs):
        rho =par[0]
        s2 =abs(par[0:len(par)])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size

        sigma = rho** abs(np.ones((n,n))-np.identity(n) )*math.sqrt(np.reshape(s2,(n,n), order='F') *
            np.reshape(s2,(n,n), order='C'))

        return sigma

    def get_gradient(self, par, times, **kwargs):
        rho = par[0]
        s2 = abs(par[0:len(par)-1])
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size

        rho_times = abs(np.ones((n,n)) - np.identity(n))
        s2_mat = math.sqrt(np.reshappe(s2, (n,n), order='C') * np.reshape(s2, (n,n), order='F'))

        d_ret = list()
        d_ret[0] = d_rho = rho_times * rho** (rho_times - 1) * s2_mat
        np.c_(np.tile(0.5 * math.sqrt(s2))[0:n], np.zeros((n, n - 1)))
        A =np.c_(np.tile(0.5 * math.sqrt(s2))[0:n],np.zeros((n, n - 1)))
        for  i in range(0,len(s2)):
            d_ret[i+1] =  rho** rho_times * (A + A.T) / math.sqrt(s2[i])
            A = np.c_[np.zeros((A.shape[0], 1)),A[:,0:-1]]

        return d_ret

    def get_param_info(self, par, times, **options):
        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return Covariance.ParamInfo(count=1 + n, name=["rho"]+['sigma2_'.format(i) for i in range(0, n)])

    def check_param(self, par, times, **kwargs):
        if par[0] > 1 or par[0] < 0:
            return False
        else:
            return True

    def get_simu_param(self, times, **options):

        n = times.shape[0] if isinstance(times, np.ndarray) else times.columns.size
        return [0.75, np.arange(1.2,1+n*0.2,0.2)]

    def est_init_param(self, pheY, pheX, pheT, **options):
        #(list, pd.core.series.Series, np.ndarray)
        n = pheT.shape[0] if isinstance(pheT, np.ndarray) else pheT.columns.size

        sel = np.isnan(pheY[:, 0]) | np.isnan(pheY[:, 1])
        rho =np.corrcoef(pheY[sel ==False, 0], pheY[sel==False, 1])
        s2 =np.nanstd(pheY, ddof=1)** 2

        return[rho, s2 * np.random.uniform(0.9, 1.1,n)]
