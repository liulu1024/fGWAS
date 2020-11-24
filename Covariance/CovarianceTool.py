import Covariance.ARMA1

covariances = []


def get_matrix(obj,par,times, *args):
    obj.get_matrix(par,times,args)


def get_gradient(obj, *args):
    obj.get_gradients(*args)


def get_param_info(obj, *args):
    obj.get_param_info(*args)


def check_param(obj, *args):
    obj.check_param(*args)


def get_simu_param(obj, *args):
    return obj.get_simu_param(times=args)


def est_init_param(obj, *args):
    obj.est_init_param(*args)


def get_all_covariances():
    return covariances


def getCovariance(type):
    covar = None
    if isinstance(type, str):
        for obj in get_all_covariances():
            if obj.covarType.upper() == type.upper():
                covar = obj

    if isinstance(type, int):
        covar = get_all_covariances()[type]

    return covar


def get_all_covariances():
    if len(covariances) == 0:
        covariances.append(Covariance.ARMA1.ARMA1(covarType="ARMA1", description="none"))
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

    return covariances


def get_covariances_count():
    return len(get_all_covariances())
