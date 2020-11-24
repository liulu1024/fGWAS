import Covariance.CovarianceTool
import numpy as np
import Covariance


class ARMA1(Covariance.BaseCovariance):
    def __init__(self, covarType, description):
        super().__init__(covarType, description)

    def get_matrix(self, par, times, *options):
        # rho,sigma2,phi
        rho = par[0]
        phi = abs(par[1])
        s2 = abs(par[2])
        n = len(times) if isinstance(times, np.ndarray) else np.size(times, 0)

        sigma = abs(s2) * ((np.ones((n, n)) - np.eye(n)) * phi + np.eye(n)) * np.power(rho,abs(
            np.full((n, n), range(0, n)) - np.full((n, n), range(0, n)).T))
        return sigma

    def get_gradient(self, par, times, options=list()):
        rho = par[1]
        phi = abs(par[2])
        s2 = abs(par[3])
        n = len(times) if isinstance(times, list) else np.size(times, 1)

        rho_times = abs(np.arange(1, n + 1).reshape(n, n) - np.arange(1, n + 1).reshape(n, n).T)

        drho = abs(s2) * ((np.ones(n, n) - np.eye(n)) * phi + np.eye(n)) * rho_times * rho ^ (
                rho_times - 1)
        dphi = abs(s2) * ((np.ones(n, n) - np.eye(n)) + 0) * rho ^ rho_times
        ds2 = ((np.ones(n, n) - np.eye(n)) * phi + np.eye(n)) * rho ^ rho_times

        return list(drho, dphi, ds2)

    def get_param_info(self, par, times, *options):
        return {'count': 3, 'names': ["rho", "phi", "sigma2"]}

    def check_param(self, par, times, options=list()):
        if par[1] > 1 or par[1] < 0 or par[2] > 1 or par[2] < 0:
            return False
        else:
            return True

    def get_simu_param(self, times, *options):
        return [0.75, 0.9, 1.2]

    def est_init_param(self, pheY, pheX, pheT, *options):
        # TODO nan 值处理
        s2 = np.var(pheY[:, 1])
        # ! & !  and !(&)
        rho = np.corrcoef(pheY[not (np.isnan(pheY[:, 3]) & np.isnan(pheY[:, 2])), 3],
                          pheY[not (np.isnan(pheY[:, 3]) & np.isnan(pheY[:, 2])), 2])
        phi = np.corrcoef(pheY[not (np.isnan(pheY[:, 2]) & np.isnan(pheY[:, 1])), 2],
                          pheY[not (np.isnan(pheY[:, 2]) & np.isnan(pheY[:, 1])), 1]) / rho

        return [rho, phi, s2 * np.random.normal(0.9, 1.1, 1)]
