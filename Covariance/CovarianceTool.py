import Covariance.AR1
import Covariance.ARH1
import  Covariance.ARMA1
import  Covariance.CS
import Covariance.CSH
import Covariance.FA1
import  Covariance.FAH1
import Covariance.VS

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
        covariances.append(Covariance.AR1.AR1(covarType="AR1", description="none"))
        covariances.append(Covariance.ARH1.ARH1(covarType="ARH1", description="none"))
        covariances.append(Covariance.CS.CS(covarType="CS", description="none"))
        covariances.append(Covariance.CSH.CSH(covarType="CSH", description="none"))
        covariances.append(Covariance.FA1.FA1(covarType="FA1", description="none"))
        covariances.append(Covariance.FAH1.FAH1(covarType="FAH1", description="none"))
        covariances.append(Covariance.VS.VS(covarType="VS", description="none"))
    return covariances


def get_covariances_count():
    return len(get_all_covariances())
