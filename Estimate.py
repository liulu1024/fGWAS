import random

import numpy as np
import pandas as pd
import warnings

from reportlab.pdfgen.canvas import Canvas

import Covariance.CovarianceTool as covarTool
import Curve.CurveTool as curveTool
import Curve as curve
import scipy.stats
import math
import ResultVo as result
import FgCommon as comm

from reportlab.graphics.shapes import *
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics import renderPDF
from sklearn.linear_model import LinearRegression
from scipy.optimize import OptimizeResult
from scipy.optimize import minimize as minimize


def fg_dat_est(obj_phe, curve_type, covariance_type, file_plot_pdf, **kwargs):
    default_options = {"max_optim_failure": 100, "min_optim_success": 20, "R2_loop": 5, "verbose": False}
    default_options.update(kwargs)
    kwargs = default_options

    intercept = False if obj_phe.intercept is None else obj_phe.intercept
    pheX = obj_phe.pheX
    if intercept:

        if pheX is None:
            pheX = np.full((obj_phe.pheY.index.size, 1), 1)
        else:
            pheX = np.insert(pheX, 0, values=1, axis=1)

    r_test = None
    if obj_phe.obj_curve is None or obj_phe.obj_curve.curveType.upper() != curve_type.upper():

        print("  No curve or new curve type is specified, curve fitting is being performed.\n")
        r_test = fg_fit_curve(obj_phe.pheY, pheX, obj_phe.pheT, curve_type=curve_type, file_plot_pdf=file_plot_pdf,
                              **kwargs)

        if r_test.error:
            sys.exit("? No curve is fitted to this data.")

        obj_phe.obj_curve = curveTool.getCurve(r_test.type)
        obj_phe.summary_curve = r_test

    print("  Curve Type==> %s <==\n" % obj_phe.obj_curve.curveType)
    par_init = r_test.par if r_test is not None else None
    r = proc_est_curve(obj_phe.pheY, pheX, obj_phe.pheT, obj_phe.obj_curve, par_init=par_init, **kwargs)
    if r.error:
        sys.exit("? Can't estimate the parameters of mean vector according to the curve function")

    parX_len = pheX.columns.size
    if pheX is None: parX_len = 0

    print(" Parameter range estimation ...... \n")
    kwargs["max_optim_failure"] = 2
    kwargs["min_optim_success"] = 2
    range = proc_est_curve_range(obj_phe.pheY, pheX, obj_phe.pheT, obj_phe.obj_curve, r.par, **kwargs)
    obj_phe.est_curve =result.EstCurve(obj_phe.obj_curve.curveType,
                         intercept,(r.par if parX_len == 0 else r.par[parX_len:len(r.par)]),
                         (r.par if parX_len == 0 else range.lower[parX_len:len(range.lower)]),
                         (r.par if parX_len == 0 else range.upper[parX_len:len(range.upper)]),
                                       (list() if parX_len == 0 else r.par[0:parX_len]),
                                        (list() if parX_len == 0 else range.lower[0:parX_len]),
                                         (list() if parX_len == 0 else range.upper[0:parX_len]),
                         r.R2)

    if r.R2 < - 0.1:
        print(
            "! The R2 value for curve fitting indiprintes the curve type is inappropriate for the phenotype data.(R2=",
            r.R2, ")\n")

    if obj_phe.obj_covar is None or obj_phe.obj_covar.covarType.upper() != covariance_type.upper():

        r_test = fg_fit_covar(obj_phe.pheY, pheX, obj_phe.pheT, r.y_resd, obj_phe.obj_curve, covariance_type,
                              **kwargs)

        if r_test.error:
            sys.exit("? No covariance is fitted to this data.")
        else:
            obj_phe.summary_covariance = r_test
            obj_phe.obj_covar = covarTool.getCovariance(r_test.type)

    print("  covariance_type==>%s <==\n" % obj_phe.obj_covar.covarType)

    kwargs["max_optim_failure"] = 2
    kwargs["min_optim_success"] = 2
    r_est = proc_est_covar(r.y_resd, None, obj_phe.pheT, obj_phe.obj_curve, obj_phe.obj_covar, **kwargs)
    if r_est.error:
        sys.exit("? Can't estimate the parameters of covariance structure according to the curve function")

    obj_phe.est_covar = result.EstCovar(obj_phe.obj_covar.covarType, r_est.par)
    obj_phe.error = False

    return obj_phe


def fn_get_resd(pheY, pheX, pheT, obj_curve, parin):
    par_X = list()
    par_c = parin
    if pheX is not None:
        par_X = parin[0:pheX.columns.size]
        if len(par_c) > pheX.columns.size:
            par_c = par_c[pheX.columns.size:len(par_c)]

    mu_gen = obj_curve.get_curve_formula(par_c, pheT, max_time=np.nanmax(pheT.values), min_time=np.nanmin(pheT.values))
    if mu_gen is None:
        return None

    if pheX is not None:
        X = np.repeat(np.dot(pheX.values, par_X), pheY.columns.size).reshape((-1, pheY.columns.size), order='F')
    else:
        X = 0

    # todo 此处有不同操作后续可进行更改
    y_resd = pheY.values - mu_gen.values - X

    return y_resd


def get_R2(pheY, pheX, pheT, obj_curve, h_par, snp_vec=None):
    pheY_hat = np.full((pheY.index.size, pheY.columns.size), np.nan)
    pheY_hat_vec = pheY_hat

    pheT_vec = pheT.values
    pheY_vec = pheY.values
    for ti in np.unique(pheT_vec):
        if ti is not np.nan:
            pheY_hat_vec[pheT_vec == ti] = np.nanmean(pheY.values[np.where(pheT_vec == ti)])
    pheY_hat = pheY_hat_vec

    if snp_vec is None:
        y_resd = fn_get_resd(pheY, pheX, pheT, obj_curve, h_par)
    else:
        y_resd = np.array([])
        par_X = None if pheX is None else h_par[0:pheX.columns.size]
        par_curve = h_par if pheX is None else h_par[pheX.columns.size:len(h_par)]
        len_curve = obj_curve.get_param_info(pheT).count

        for i in range(0, 3):

            idx = np.where(snp_vec == i)[0].tolist()
            if len(idx) > 0:
                hx_par = par_X
                hx_par = np.insert(hx_par,hx_par.size,par_curve[0:len_curve])
                par_curve = par_curve[len_curve:len(par_curve)]
                if y_resd.size==0:
                    y_resd=fn_get_resd(pheY.iloc[idx,:], pheX.iloc[idx,:], pheT.iloc[idx,:], obj_curve, hx_par)
                else:
                    y_resd = np.r_[y_resd, fn_get_resd(pheY.iloc[idx,:], pheX.iloc[idx,:], pheT.iloc[idx,:], obj_curve, hx_par)]
    R2 = 1 - (np.nansum(y_resd ** 2) / (
        1 if np.nansum(pheY.values - pheY_hat) == 0 else np.nansum(pheY.values - pheY_hat) ** 2))
    return R2


def proc_est_curve(pheY, pheX, pheT, obj_curve, par_init, **kwargs):
    kwargs["max_time"] = np.nanmax(pheT)
    kwargs["min_time"] = np.nanmin(pheT)
    if "R2_loop" not in kwargs.keys(): kwargs["R2_loop"] = None

    def get_init_covariate_par(pheY, pheX):
        pheX.columns = ["x{0}".format(i) for i in range(0, pheX.shape[1])]
        # list
        fit = LinearRegression(fit_intercept=False).fit(pheX, np.mean(pheY, axis=1))
        if fit.coef_ is not None:
            return fit.coef_
        else:
            return np.nanmean(pheY, axis=1) / np.nanmean(pheX, axis=0)

    def get_init_curve_par(pheY, pheX, pheT, f_obj):

        par_X = []
        YX = 0
        if pheX is not None:
            par_X = get_init_covariate_par(pheY, pheX)
            YX = np.repeat([np.dot(pheX, par_X)], pheY.shape[1])

        par_curve = f_obj.est_init_param(pheY, pheX, pheT, **kwargs)
        return np.append(par_X, par_curve)

    def get_rand_curve_par(pheY, pheX, pheT, f_obj, parin):

        uc = check_fittness(pheY, pheX, pheT, f_obj, parin)

        if uc["fit"]:

            return parin * np.random.uniform(0.9, 1.1, (len(parin), 1))

        else:

            if uc["over_0.05"] < pheY.columns.size / 4:
                return parin * np.random.uniform(0.5, 1.5, (len(parin), 1))
            else:

                par_curve = get_init_curve_par(pheY, pheX, pheT, f_obj)
                return par_curve * np.random.uniform(0.5, 1.5, (len(par_curve, 1)))

    def fn_mle_est(parin, extra_par):
        y_resd = fn_get_resd(extra_par["pheY"], extra_par["pheX"], extra_par["pheT"], extra_par["obj_curve"], parin)
        if y_resd is None:
            return None
        A = np.nansum(y_resd ** 2)
        return A

    if par_init is None:
        par_init = get_init_curve_par(pheY, pheX, pheT, obj_curve)

    if comm.DEBUG: print("Initial parameters for curve: ", par_init, "\n")

    R2_max = float("-inf")
    R2_success = 0
    R2_failed = 0
    R2_list = None
    # todo R2_fail<500 测试先改为5

    while (R2_success < (2 if kwargs["R2_loop"] is None else kwargs["R2_loop"])) and R2_failed <= 1:

        comm.reset_seed()
        h0 = proc_mle_loop(pheY, pheX, pheT, obj_curve, fn_mle_est,
                           mle_extra_par={"pheY": pheY, "pheX": pheX, "pheT": pheT, "obj_curve": obj_curve},
                           parin=par_init, fn_init_par=get_init_curve_par, fn_rand_par=get_rand_curve_par,
                           options=kwargs)

        if h0 is None or h0.value is None or h0.value == float("inf"):
            R2_failed = R2_failed + 1
            continue

        else:

            R2 = get_R2(pheY, pheX, pheT, obj_curve, h0.par)
            if comm.DEBUG: print("MLE results for curve: R2=", R2, h0.value, h0.par, "\n")

            if R2 > 1 or R2 < -0.1:
                R2_failed = R2_failed + 1
            else:
                R2_success = R2_success + 1

            if R2 <= R2_max: continue

            R2_max = R2
            par_init = h0.par * np.random.uniform(0.99, 1, (1, len(h0.par)))

            y_resd = fn_get_resd(pheY, pheX, pheT, obj_curve, h0.par)
            ## The sum of the squared differences between each observation and its predicted value.
            y_SSE = np.nansum(y_resd ** 2)
            ##Gives the average of the squares of the errors of each value.
            y_MSE = y_SSE / len(np.isnan(pheY.values) == False)
            ##The square root of the MSE that estimates the standard deviation of the random error.
            y_RMSE = math.sqrt(y_MSE)
            ## pramater count
            K = obj_curve.get_param_info(pheT).count
            n_sample = pheY.index.size

            AIC = 2 * K + n_sample * math.log(y_SSE / n_sample)
            AICc = math.log(y_SSE / n_sample) + (n_sample + K) / (n_sample - K - 2)
            BIC = n_sample * math.log(y_SSE / n_sample) + K * math.log(n_sample)

            R2_list = result.R2List(error=False, par_count=K, AIC=AIC, AICc=AICc, BIC=BIC, SSE=y_SSE, MSE=y_MSE,
                                    RMSE=y_RMSE, R2=R2,
                                    par=h0.par, y_resd=y_resd)

    if R2_list is None:
        return result.R2List(error=True, par=None, val=None, R2=None)
    else:
        return R2_list


def proc_est_curve_range(pheY, pheX, pheT, f_curve, par_init, **kwargs):
    n_obs = pheY.index.size
    mu_pars_temp = list()

    loop = 0
    while loop < comm.MU_RANGE_LOOP:
        num = np.around(n_obs * np.random.uniform(0.5, 0.9, 1), 0).astype(int)
        y_sub = random.sample(range(0, n_obs), n_obs)[0:num[0]]

        pheX0 = None
        # drop =F
        pheY0 = pheY.iloc[y_sub]
        if pheX is not None: pheX0 = pheX.iloc[y_sub]
        pheT0 = pheT.iloc[y_sub] if pheT is not None else pheT

        r = proc_est_curve(pheY0, pheX0, pheT0, f_curve, par_init, **kwargs)
        if r.error:
            continue

        mu_pars_temp.append(r.par)
        loop = loop + 1

    mu_pars = np.array(mu_pars_temp)
    mu_lower = []
    mu_upper = []

    for i in range(0, mu_pars.shape[1]):
        mu_lower.append(np.min(mu_pars[:, i]))
        mu_upper.append(np.max(mu_pars[:, i]))

    return result.CurveRange(mu_lower, mu_upper)


def proc_est_covar(y_resd, pheX, pheT, obj_curve, obj_covar, par_init=None, **kwargs):
    def get_init_covar_par(y_resd, pheX, pheT, f_covar):

        return f_covar.est_init_param(y_resd, pheX, pheT, **kwargs)

    def get_rand_covar_par(y_resd, pheX, pheT, f_covar, parin):

        return parin * np.random.uniform(0.9, 1.1, len(parin))

    # parin:
    # phi1, s1, phi2, s2
    def fn_mle_est(parin, extra_par):

        y_resd = extra_par["y_resd"]
        pheX = extra_par["pheX"]
        pheT = extra_par["pheT"]
        f_covar = extra_par["f_covar"]

        cov_mat = f_covar.get_matrix(parin, pheT)
        # 判断是否为半正定矩阵
        if(not np.all(np.linalg.eigvals(cov_mat)>0)):
            return None
        if np.isnan(cov_mat).any():
            return None

        pv = comm.dmvnorm_fast(y_resd, np.repeat(0, y_resd.shape[1]), cov_mat, None,log=True)
        if np.isinf(pv).any():
            return None

        A = np.nansum(pv)
        return -A

    if par_init is None:
        par_init = get_init_covar_par(y_resd, pheX, pheT, obj_covar)

    if comm.DEBUG:  print("Initial parameters for covariance: ", par_init, "\n")

    h0 = proc_mle_loop(y_resd, pheX, pheT, obj_covar, fn_mle_est,
                       mle_extra_par={"y_resd": y_resd, "pheX": pheX, "pheT": pheT, "f_covar": obj_covar},
                       parin=par_init, fn_init_par=get_init_covar_par, fn_rand_par=get_rand_covar_par, options=kwargs)

    if (comm.DEBUG):  print("MLE results for curve covariance:", h0.value, h0.par, "\n")

    if h0.value is not None and (not math.isinf(h0.value)):
        K = obj_covar.get_param_info(pheT).count
        AIC = 2 * K + 2 * h0.value
        BIC = 2 * h0.value + K * np.log(y_resd.shape[0])
        return result.ProcCurve(False, AIC=AIC, BIC=BIC, par=h0.par, logL=-h0.value)

    else:
        return result.ProcCovar(error=True, par=None, err_info="Failed to estimate covariance structure.")


# vo
# to do  options={"min_optim_success": 5, "max_optim_failure": 100}
def proc_mle_loop(pheY, pheX, pheT, f_obj, fn_mle, mle_extra_par, parin, fn_init_par, fn_rand_par,
                  options={"min_optim_success": 5, "max_optim_failure": 100}):
    init_optim = 1
    fun = float("inf")
    #达到1e-8精度耗时
    # control = result.control(50000, 1e-8)
    control = result.control(500, 1e-1)
    h0 = scipy.optimize.OptimizeResult()
    while math.isinf(fun):
        try:
            # fun 最小值目标函数值
            # x
            # success
            #method = ("nelder-mead" if random.uniform(0, 1) > 0.75 else "bfgs"), tol = control.reltol
            h0 = minimize(fun=fn_mle, x0=parin, args=((mle_extra_par)),
                          method="nelder-mead", tol=control.reltol,
                          options={"maxiter": control.maxit})
            if len(h0.keys()) == 0: continue
            comm.TRY_SILENT = False
            if h0.fun is None or h0.success is None or (not h0.success):
                raise Exception
            fun = h0.fun
        except:
            if h0 is not None and len(h0.keys()) != 0:
                if h0.status == 1:
                    control.maxit = control.maxit * 2
                    if control.maxit > 500 * 4096:
                        control.maxit = 500 * 4096
                        #control.reltol = control.reltol * 2

                else:
                    if comm.DEBUG: print("h0.success=", h0.success, "\n")

            comm.reset_seed()
            parin = fn_init_par(pheY, pheX, pheT, f_obj)

            init_optim = init_optim + 1
            if init_optim > 100: return result.MleLoop(True)
            continue
        comm.reset_seed()
    if comm.DEBUG: print("MLE[0]", h0.fun, h0.x, "\n")

    parin0 = h0.x
    n_optim_failure = 0
    n_optim_success = 0

    while n_optim_failure < options["max_optim_failure"] and n_optim_success < options["min_optim_success"]:

        parinx = fn_rand_par(pheY, pheX, pheT, f_obj, parin0)
        h2 = scipy.optimize.OptimizeResult()
        try:
            #method = "nelder-mead" if (n_optim_failure + n_optim_success) % 2 == 1 else "bfgs"
            h2 = minimize(fun=fn_mle, x0=parinx, args=((mle_extra_par)),
                          method="nelder-mead",
                          tol=control.reltol, options={"maxiter": control.maxit})
            comm.TRY_SILENT = False
            if len(h2.keys()) == 0: continue
            if h2 is None or (not h2.success):
                raise Exception
            if h2.fun < h0.fun:
                h0 = h2
            n_optim_success = n_optim_success + 1
            comm.reset_seed()
        except:
            if len(h2.keys()) != 0:
                if not h2.success:
                    control.maxit = control.maxit * 1.5
                    #500*4096
                    if control.maxit > 500*4096:
                        control.maxit = 500*4096
                        control.reltol = control.reltol * 2
            n_optim_failure = n_optim_failure + 1
            comm.reset_seed()
            continue

    if comm.DEBUG: print("MLE[F]", h0.fun, h0.x, "\n")
    return result.MleLoop(error=False, par=h0.x, value=h0.fun)


def check_fittness(pheY, pheX, pheT, f_obj, parin):
    if not issubclass(f_obj.__class__, curve.BaseCurve):
        return {"fit": True}

    y_resd = fn_get_resd(pheY, pheX, pheT, f_obj, parin)
    if y_resd is None:
        return {"fit": True}
    # 按列计算标准差
    y_sd = np.std(pheY, 0)
    py_prob = 1 - scipy.stats.norm(0, y_sd).cdf(np.mean(y_resd, axis=0))
    py_05 = len(np.where(py_prob < 0.05))

    return {"fit": False if py_05 > 1 else True, "over_0.05": py_05}


def fg_fit_curve(pheY, pheX, pheT, curve_type, file_plot_pdf, **kwargs):
    print(" Searching curve type ....s..\n")

    obj_curves = list()
    if curve_type == "auto" or curve_type == "" or curve_type is None:
        for i in range(0, curveTool.get_curve_count()):
            obj_curves.append(curveTool.get_curve(i))

    else:
        for i in range(0, len(curve_type)):
            obj_curves.append(curveTool.get_curve(curve_type[i]))
    est_curve = []
    index = 0
    c = Canvas(file_plot_pdf)
    for obj_curve in obj_curves:
        print("* [%d/%d] try curve: %s\n" % (index, len(obj_curves), obj_curve.curveType))
        r = proc_est_curve(pheY, pheX, pheT, obj_curve, None, **kwargs)
        if r.error:
            warnings.warn(
                "Can't estimate the parameters of mean vector according to the curve function[ curve_type=%s]" % obj_curve.curveType)
            index = index + 1
        else:
            r.type = obj_curve.curveType
            est_curve.append(r)
            index = index + 1

            if file_plot_pdf is not None:
                drawing = Drawing(400, 200)
                lp = LinePlot()
                lp.x = 50
                lp.y = 50
                lp.height = 400
                lp.width = 200
                lp.xValueAxis._valueMin = np.nanmin(pheT.values)
                lp.xValueAxis._valueMax = np.nanmax(pheT.values)
                lp.yValueAxis._valueMin = np.nanmin(pheY.values)
                lp.yValueAxis._valueMax = np.nanmax(pheY.values)

                lp.xValueAxis.valueSteps = [i for i in range(np.natnmin(pheT.values), np.nanmax(pheT.values) + 1)]
                lp.yValueAxis.valueSteps = [i for i in range(np.nanmin(pheY.values), np.nanmax(pheY.values) + 1)]
                # 标题
                title = String(200, 220, r.type.astype(str))
                title.fontName = 'msyh'
                title.fontSize = 14
                title.fillColor = colors.gray
                title.textAnchor = 'middle'

                # drawing.add(lp)

                pheT_vec = pheT.values
                pheY_vec = pheY.values
                ti = np.unique(pheT.values)
                ti = ti[np.isnan(ti) == False]
                y_mu = [np.nan for i in range(0, len(ti))]
                for i in range(0, len(ti)):
                    y_mu[i] = np.nanmean(pheY_vec[pheT_vec == ti[i]])

                # draw.add(PolyLine(list(zip(ti,y_mu)),strokerColor=colors.black,strokerWidth=1))

                mu_X = 0
                par_c = r.par
                if pheX is not None:
                    par_X = r.par[0:pheX.columns.size]
                    par_c = r.par[pheX.columns.size:len(r.par)]
                    mu_X = np.nanmean(np.dot(pheX.values, par_X))

                ti = np.arange(np.nanmin(pheT.values), np.nanmax(pheT.values) + 1, 1)
                mu_gen = obj_curve.get_curve_formula(par_c, ti, max_time=np.nanmax(pheT), min_time=np.nanmin(pheT))

                # data
                data = []
                for i in range(0, pheY.index.size):
                    data.append(list(zip(pheT[i, :], pheY[i, :])))
                data.append(list(zip(ti, y_mu)))
                data.append(list(zip(ti, mu_gen + mu_X)))
                lp.data = data
                for i in range(0, pheY.index.size):
                    lp.lines[i].strokerColor = colors.gray
                    lp.lines[i].strokerWidth = 0.2
                lp.lines[-2].strokerColor = colors.black
                lp.lines[-2].strokerWidth = 1
                lp.lines[-1].strokerColor = colors.red
                lp.lines[-1].strokerWidth = 1.5
                drawing.add(lp)
                renderPDF.draw(drawing, c, 100, 100, showBoundary=False)
                c.showPage()
    c.save()
    if len(est_curve) > 0:
        data = {"type": [est_curve[i].type for i in range(0, len(est_curve))],
                "parm": [est_curve[i].par_count for i in range(0, len(est_curve))],
                "AIC": [est_curve[i].AIC for i in range(0, len(est_curve))],
                "AICc": [est_curve[i].AICc for i in range(0, len(est_curve))],
                "BIC": [est_curve[i].BIC for i in range(0, len(est_curve))],
                "SSE": [est_curve[i].SSE for i in range(0, len(est_curve))],
                "MSE": [est_curve[i].MSE for i in range(0, len(est_curve))],
                "RMSE": [est_curve[i].RMSE for i in range(0, len(est_curve))],
                "R2": [est_curve[i].R2 for i in range(0, len(est_curve))]}
        est_summary = pd.DataFrame(data)
    else:
        est_summary = None
        sys.exit("Failed to do curve fitting.\n")
    ## use AIC to determine the best curve type.
    print("  Curve Fitting Summary:\n")
    print(est_summary.head(132))
    fit_idx = est_summary.iloc[:, 2].idxmin()
    return result.FitCurveVo(False, est_curve[fit_idx].curveType, est_curve[fit_idx].y_resd, est_curve[fit_idx].par,
                             est_summary, est_curve)


def fg_fit_covar(pheY, pheX, pheT, y_resd, obj_curve, covariance_type="auto", **kwargs):
    print(" Searching covariance matrix .......\n")

    obj_covars = list()
    if covariance_type == "auto" or covariance_type == "" or covariance_type is None:
        for i in range(0, covarTool.get_covariances_count()):
            obj_covars[i] = covarTool.getCovariance(i)

    else:
        obj_covars.append(covarTool.getCovariance(covariance_type))
        # for i in range(0, len(covariance_type)):
        #     obj_covars[i] = covarTool.getCovariance(covariance_type[i])

    est_covar = list()
    index = 0
    for obj_covar in obj_covars:

        print(
            " *[%d/%d]: try covariance matrix: %s\n" % (index, covarTool.get_covariances_count(), obj_covar.covarType))

        r_est = proc_est_covar(y_resd, None, pheT, obj_curve, obj_covar, **kwargs)
        if r_est.error:
            warnings.warn(
                "Can't estimate the parameters of covariance structure according to the covariance function[ covariance_type=%s]" % obj_covar.covarType)
        else:
            r_est.type = obj_covar.covarType
            est_covar.append(r_est)
            index = index + 1

    if len(est_covar) >= 1:
        data = {"type": [est_covar[i].type for i in range(0, len(est_covar))],
                "L": [est_covar[i].logL for i in range(0, len(est_covar))],
                "AIC": [est_covar[i].AIC for i in range(0, len(est_covar))],
                "BIC": [est_covar[i].BIC for i in range(0, len(est_covar))]}
        est_summary = pd.DataFrame(data)
        print("Covariance Estimation Summary:\n")
        print(est_summary.head(132))

        fit_idx = est_summary.iloc[:, 2].idxmin()
        return result.FitCovarVo(False, est_covar[fit_idx].type, est_covar[fit_idx].par, est_summary)

    else:
        return result.FitCovarVo(error=True)
