'''
err
err_info
'''


class NormalResultVo(object):
    def __init__(self, error, err_info):
        self.error = error
        self.err_info = err_info


class SnpCurveVo(object):
    def __init__(self, error=None, result=None):
        self.error = error
        self.result = result


class MinDist(object):
    def __init__(self, min_grp=None, dist=None):
        self.min_grp = min_grp
        self.dist = dist


class EstCurve(object):
    def __init__(self, type=None, intercept=None,param=None, param_lower=None, param_upper=None, parX=None, parX_lower=None,
                 parX_upper=None, R2=None):
        self.type = type
        self.intercept=intercept
        self.param = param
        self.param_lower = param_lower
        self.param_upper = param_upper
        self.parX = parX
        self.parX_lower = parX_lower
        self.parX_upper = parX_upper
        self.R2 = R2
class EstCovar(object):
    def __init__(self,type=None,param=None):
        self.type=type
        self.param=param

# est
class FitCovarVo(object):
    def __init__(self, error=None, type=None, par=None, summary=None):
        self.error = error
        self.type = type
        self.par = par
        self.summary = summary


class FitCurveVo(object):
    def __init__(self, error=None, type=None, y_resd=None, par=None, summary=None, est_curve=None):
        self.error = error
        self.typr = type
        self.y_resd = y_resd
        self.par = par
        self.summary = summary
        self.est_curve = est_curve


class CurveRange(object):
    def __init__(self, lower=None, upper=None):
        self.lower = lower
        self.upper = upper


class R2List(object):
    def __init__(self, error=None, par_count=None, AIC=None,
                 AICc=None, BIC=None, SSE=None, MSE=None, RMSE=None, R2=None,
                 par=None, y_resd=None, val=None):
        self.error = error
        self.par_count = par_count
        self.AIC = AIC
        self.AICc = AICc
        self.BIC = BIC
        self.SSE = SSE
        self.MSE = MSE
        self.RMSE = RMSE
        self.R2 = R2
        self.par = par
        self.y_resd = y_resd
        self.val = val


class H0(object):
    def __init__(self, x=None, fun=None, success=None):
        self.x = x
        self.fun = fun
        self.success = success


class MleLoop(object):
    def __init__(self, error=None, par=None, value=None):
        self.error = error
        self.par = par
        self.value = value


class ProcCurve(object):
    def __init__(self, error=None, AIC=None, BIC=None, par=None, logL=None, type=None):
        self.error = error
        self.AIC = AIC
        self.BIC = BIC
        self.par = par
        self.logL = logL
        self.type = type


class ProcCovar(object):
    def __init__(self, error=None, par=None, err_info=None):
        self.error = error
        self.par = par
        self.err_info = err_info


class control(object):
    def __init__(self, maxit=None, reltol=None):
        self.maxit = maxit
        self.reltol = reltol


class AllInOneVo(object):
    def __init__(self, x=None, fun=None, R2=None,pv=None,ssr=None):
        self.x = x
        self.fun = fun
        self.R2 = R2
        self.pv=pv
        self.ssr=ssr
