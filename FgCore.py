import math
import sys

import numpy as np
import pandas as pd
import scipy

import FgCommon as comm
import ResultVo as result
import Estimate as estimate
from scipy.stats import chi2_contingency
from scipy.stats import multivariate_normal
import scipy.optimize
import ResultVo as resultVo
import numdifftools as nd
from scipy.optimize import OptimizeResult


def proc_est_h0(obj_phe, **kwargs):
    pheY = obj_phe.pheY.values
    pheT = obj_phe.pheT.values
    pheX = obj_phe.pheX.values
 # todo 为NOne
    if obj_phe.intercept:

        if obj_phe.pheX is not None:
            pheX = np.c_[pheX, np.ones((pheX.shape[0], 1))]
        else:
            pheX = np.ones((obj_phe.pheY.index.size, 1))

    kwargs["max_time"] = np.nanmax(pheT)
    kwargs["min_time"] = np.nanmin(pheT)
    if kwargs["min_optim_success"] is None: kwargs["min_optim_success"] = 5
    kwargs["b_permu"] = False

    h0 = proc_full_h0(obj_phe, pheY, pheX, pheT,h0_ref=None, **kwargs)

    n_par_curve = obj_phe.obj_curve.get_param_info(pheT, **kwargs).count

    parin_X = None
    parin_curve = h0.x
    if pheX is not None:
        parin_X = h0.x[0:pheX.shape[1]]
        parin_curve = h0.x[pheX.shape[1]:len(h0.x)]

    parin_covar = parin_curve[n_par_curve:len(parin_curve)]
    parin_curve = parin_curve[0:n_par_curve]

    X = 0
    if pheX is not None:
        X = np.tile(np.dot(pheX, parin_X).flatten(order='F'), pheY.shape[1]).reshape((-1, pheY.shape[1]), order='F')

    mu = obj_phe.obj_curve.get_curve_formula(parin_curve, pheT, **kwargs)
    Y_delt = pheY - mu - X
    mat_cov = obj_phe.obj_covar.get_matrix(parin_covar, pheT)

    h0.pv = comm.dmvnorm_fast(Y_delt, np.repeat(0, mat_cov.shape[0]), mat_cov, nna_vec=None, log=True),

    h0.ssr = np.nansum(Y_delt ** 2)

    return h0


def proc_mle(index, snp_info, snp_vec, obj_phe, intercepted, **kwargs):
    default_kwargs = {"b_permu": False}
    default_kwargs.update(kwargs)
    kwargs = default_kwargs

    gen_par = len(np.unique(snp_vec[~np.isnan(snp_vec)]))

    pheY = np.array(obj_phe.pheY)
    pheT = np.array(obj_phe.pheT)
    pheX = np.array(obj_phe.pheX)

    if intercepted:
        if obj_phe.pheX is not None:
            pheX = np.c_(np.ones((pheX.shape[0], 1)), pheX)
        else:
            pheX = np.ones((pheY.shape[0], 1))
    # find h0
    h0_ref = np.nansum(-1 * obj_phe.h0.pv)
    snp_miss = np.where((np.isnan(snp_vec) == True) + (snp_vec < 0) + (snp_vec > 2))[0]
    snp_miss0 = snp_miss
    snp_vec0 = np.where(snp_vec == 0)[0]
    snp_vec1 = np.where(snp_vec == 1)[0]
    snp_vec2 = np.where(snp_vec == 2)[0]
    ## if the genotye is less than 5, curve fitting is difficult.
    if len(snp_vec0) < 5:
        snp_miss = np.hstack((snp_miss, snp_vec0))
    if len(snp_vec1) < 5:
        snp_miss = np.hstack((snp_miss, snp_vec1))
    if len(snp_vec2) < 5:
        snp_miss = np.hstack((snp_miss, snp_vec2))

    if len(snp_miss) > 0:
        if pheX is not None:
            pheX = np.delete(pheX, snp_miss, axis=0)
            pheY = np.delete(pheY, snp_miss, axis=0)
            pheT = np.delete(pheT, snp_miss, axis=0)
            snp_vec = np.delete(snp_vec, snp_miss)
            h0_ref = np.nansum(-1 * np.delete(obj_phe.h0.pv, snp_miss))

    kwargs["max_time"] = np.nanmax(pheT)
    kwargs["min_time"] = np.nanmin(pheT)
    kwargs["max_optim_failure"] = 20

    if kwargs["opt_method"].upper() == "FAST" or kwargs["opt_method"].upper() == "FAST-NORM":
        if kwargs["min_optim_success"] is None:
            kwargs["min_optim_success"] = 1

        h0 = proc_fast_h0(obj_phe, pheY, pheX, pheT, **kwargs)

        if kwargs["min_optim_success"] is None:
            kwargs["min_optim_success"] = 2
        h1 = proc_fast_h1(obj_phe, pheY, pheX, pheT, snp_vec, h0, **kwargs)

    bOptim = True
    if (h0 is None) or (h0.fun is None) or (math.isinf(h0.fun)):
        bOptim = False
    if isinstance(h1,resultVo.AllInOneVo):
         h1_value=h1.fun
    else:
         h1_value=h1.value
    if (h1 is None) or (h1_value is None) or (math.isinf(h1_value)):
        bOptim = False

    if bOptim:
        #todo 未测试
        r_val = (h0.fun - h1.value) * 2
        # 卡方检验
        r_pv = 1 - scipy.stats.chi2.cdf(r_val, df=(gen_par - 1) * obj_phe.obj_curve.get_param_info(pheT).count)

    else:
        r_val = r_pv = None

    maf = 1 - np.mean(snp_vec) / 2 if (np.mean(snp_vec) / 2 > 0.5) else np.mean(snp_vec) / 2

    col=len(snp_info[0:3])
    p_1=pd.DataFrame(np.array([index]))
    p_1=pd.concat([p_1,pd.DataFrame(snp_info[0:3].reshape((1,col)))],axis=1)

    p_2=np.array([maf,len(snp_miss0), len(snp_vec0), len(snp_vec1), len(snp_vec2),gen_par,r_val, r_pv])
    p_2=np.insert(p_2,p_2.size,h0.fun)
    p_2 = np.insert(p_2, p_2.size, h0.x)
    p_2 = np.insert(p_2, p_2.size, h1_value)
    p_2 = np.insert(p_2, p_2.size, h1.x)
    p_2=pd.DataFrame(p_2.reshape((1,-1)))

    re=pd.concat([p_1,p_2],axis=1)

    if comm.DEBUG:
        print("snp[", re[0, 1].astype(str), "] loci=", re[0, 2], re[0, 3], "LR2=", re[0, 10], re[0, 11], "Param=",
              re.iloc[0, 12:re.columns.size], "\n")

    return re


def proc_fast_h1_allinone(obj_phe, pheY, pheX, pheT, snp_vec, h0, parX, **kwargs):
    nna_vec = comm.get_non_na_number(pheY)
    n_par_curve = obj_phe.obj_curve.get_param_info(pheT, **kwargs).count
    if len(np.unique(snp_vec)) == 1:
        h1 = h0

        parin_X = h0.x[0:pheX.shape[1]]
        parin_covar = h0.x[(1 + pheX.shape[1] + n_par_curve):len(h0.x)]
        parin_curve = np.full((3, n_par_curve), np.nan)
        parin_curve[np.unique(snp_vec).astype(int), :] = h0.x[pheX.shape[1]:pheX.shape[1] + n_par_curve]
        h1.x = parin_X
        h1.x=np.insert(h1.x,(h1.x).size,parin_curve.flatten(order='C'))
        h1.x=np.insert(h1.x,(h1.x).size, parin_covar)
        return h1

    ## H1 is not better than H0, try to get all parameters of H1 in one optim function
    if pheX is None or parX is not None:
        parin_X = parX
        parin_curve = h0.x[0:len(obj_phe.est_curve.param)]
        parin_covar = h0.x[(len(obj_phe.est_curve.param)): len(h0.x)]

    else:
        parin_X = h0.x[0:pheX.shape[1]]
        parin_curve = h0.x[(len(parin_X)):(len(parin_X) + len(obj_phe.est_curve.param))]
        parin_covar = h0.x[(len(parin_X) + len(obj_phe.est_curve.param)): len(h0.x)]
    geno_type = len(np.unique(snp_vec))

    def optim_func(parin, sum_type=0):
        if pheX is not None:
            parin_X = parin[0:pheX.shape[1]]
            parin = np.delete(parin, range(0, pheX.shape[1]))

        parin_curve = parin[0:(3 * n_par_curve)]
        parin_covar = parin[(3 * n_par_curve):len(parin)]

        if pheX is not None:
            Y_delt = pheY - np.tile(np.dot(pheX, parin_X).flatten(order='F'), pheY.shape[1]).reshape((-1, pheY.shape[1]), order='F')
        else:
            Y_delt = pheY

        col=math.floor(parin_curve.size/geno_type)
        parin_curve = parin_curve[0:geno_type*col].reshape((geno_type, col))
        snp_idx = 0
        for k in range(0, 3):
            idx_k = np.where(snp_vec == k)[0].tolist()
            if len(idx_k) > 0:
                snp_idx = snp_idx + 1
                Y_delt[idx_k, :] = Y_delt[idx_k, :] - obj_phe.obj_curve.get_curve_formula(parin, pheT[snp_idx, :], **kwargs)

        if sum_type == 0:
            mat_cov = obj_phe.obj_covar.get_matrix(parin_covar, pheT)
            try:
                if not np.all(np.linalg.eigvals(mat_cov) > 0):
                    return None
                pv = comm.dmvnorm_fast(Y_delt, np.repeat(0, mat_cov.shape[1]), mat_cov, nna_vec, log=True)
                if pv is None or np.isnan(pv).any():
                    raise Exception
            except:
                return None

            return -np.sum(pv)


        else:
            return np.nansum(Y_delt ** 2)

    try:
        x0=parin_X
        x0=np.insert(x0,x0.size,np.tile(parin_curve, geno_type))
        x0=np.insert(x0,x0.size,parin_covar)
        # todo minimize 仍有问题
        hx = scipy.optimize.minimize(fun=optim_func, args=((0)),
                                     x0=x0,
                                     method="Nelder-Mead")
        if (hx is None) or (hx.fun is None) or (math.isinf(hx.fun)):
            raise Exception
    except:
        x0 = parin_X
        x0 = np.insert(x0, x0.size, np.tile(parin_curve.flatten(order='C'), 3))
        x0 = np.insert(x0, x0.size, parin_covar)
        return result.AllInOneVo(x=x0, fun=None, R2=None)

    if geno_type < 3:
        snp_miss = np.where(np.isnan(np.in1d([0, 1, 2], np.unique(snp_vec))) == True)[0] - 1
        nstart = len(parin_X) + snp_miss * n_par_curve
        hx.x = [hx.x[0:nstart], np.full(n_par_curve, np.nan), hx.x[(nstart):len(hx.x)]]

    hx.ssr = optim_func(hx.x, sum_type=1)
    hx.R2 = estimate.get_R2(pd.DataFrame(pheY), pd.DataFrame(pheX), pd.DataFrame(pheT), obj_phe.obj_curve, hx.x, snp_vec)

    return hx


def proc_fast_h0(obj_phe, pheY, pheX, pheT, **kwargs):
    def proc_h0_curve(parin, pheY, pheT, pheX, snp_vec, parX, obj_curve, sum_type):
        if pheX is not None:

            parin_X = parin[0:pheX.shape[1]]
            parin = parin[pheX.shape[1]:len(parin)]
            X = np.tile(np.dot(pheX, parin_X).flatten(order='F'), pheY.shape[1]).reshape((pheX.shape[0], pheY.shape[1]),
                                                                          order='F')

        else:
            X = np.zeros((pheY.shape[0], 1))

        mu = obj_curve.get_curve_formula(parin, pheT, **kwargs)
        Y_delt = pheY - mu - X

        if sum_type == 0:
            Y_var = np.var(Y_delt, axis=1)
            Y_var[np.isnan(Y_var) == True] = 0

            pv = np.log(multivariate_normal.pdf(Y_delt, cov=np.diag(Y_var)))
            A = -np.sum(pv)

        else:
            A = Y_delt[np.isnan(Y_delt) == False]
        return A


    def proc_h0_covmatrix(parin, Y_delt, pheT, obj_covariance, nna_vec, **kwargs):
        mat_cov = obj_covariance.get_matrix(parin, pheT)

       # 判断是否为半正定矩阵
        if not np.all(np.linalg.eigvals(mat_cov)>0) :
             return None
        try:

           pv = comm.dmvnorm_fast(Y_delt, np.repeat(0, mat_cov.shape[0]), mat_cov, nna_vec, log=True)
           if pv is None or np.isnan(pv).any():
             raise Exception
        except:
            return None

        A = np.sum(pv)
        return -A

    parinx = np.array(obj_phe.est_curve.parX)
    parinx=np.insert(parinx,len(parinx),obj_phe.est_curve.param)
    h0_curve = optim_least_square(parinx, proc_h0_curve, pheY=pheY, pheT=pheT, pheX=pheX, snp_vec=None, parX=None,
                                  obj_curve=obj_phe.obj_curve, method=(kwargs["opt_method"] == "FAST"))

    if h0_curve is None or h0_curve.value is None or math.isinf(h0_curve.value):
        return result.AllInOneVo(x=[parinx, obj_phe.est_covar.param], fun=None, R2=None)

    X = 0
    parin_curve = h0_curve.x
    if pheX is not None:
        parin_X = h0_curve.x[0:pheX.shape[1]]
        parin_curve = h0_curve.x[pheX.shape[1]:len(h0_curve.x)]
        X = np.tile(np.dot(pheX, parin_X).flatten(order='F'), pheY.shape[1]).reshape((-1, pheY.shape[1]), order='F')

    nna_vec = comm.get_non_na_number(pheY)
    mu = obj_phe.obj_curve.get_curve_formula(parin_curve, pheT, **kwargs)
    Y_delt = pheY - mu - X
    h0_cm= comm.optim_BFGS(np.array([]), np.array([]), obj_phe.est_covar.param,
                            proc_h0_covmatrix,None,Y_delt=Y_delt, pheT=pheT,
                            obj_covariance=obj_phe.obj_covar,
                            nna_vec=nna_vec,funType='matrix', **kwargs)

    if h0_cm is None or h0_cm.fun is None or math.isinf(h0_cm.fun):
        x=h0_curve.x
        x=np.insert(x,x.size,obj_phe.est_covar.param)
        h0_cm = result.AllInOneVo(x=x, fun=None, R2=None,ssr=None)
    else:
        mat_cov = obj_phe.obj_covar.get_matrix(h0_cm.x, pheT)
        h0_cm.pv = comm.dmvnorm_fast(Y_delt, np.repeat(0, mat_cov.shape[0]), mat_cov, nna_vec, log=True)
        h0_cm.x =h0_curve.x
        h0_cm.x =np.insert(h0_cm.x,(h0_cm.x).size,h0_cm.x)
        h0_cm.ssr = np.nansum(Y_delt ** 2)
        h0_cm.R2 = estimate.get_R2(pd.DataFrame(pheY), pd.DataFrame(pheX), pd.DataFrame(pheT), obj_phe.obj_curve, h0_cm.x)

    if not kwargs["b_permu"]:
        print("H0[F]", h0_cm.fun, h0_cm.ssr, "Param=", h0_cm.x, "\n")

    return h0_cm


def proc_fast_h1_comb(obj_phe, pheY, pheX, pheT, snp_vec, h0, snp_comb=np.array([0, 1, 2]), parX=None, **kwargs):
    def proc_h1_curve(parin, pheY, pheT, pheX, snp_vec, parX, obj_curve, sum_type):
        if pheX is not None and parX is  None:

            parin_X = parin[0:pheX.shape[1]]
            parin_curve = parin[pheX.shape[1]:len(parin)]
        else:
            parin_curve = parin
            parin_X = parX

        n_par_curve = obj_curve.get_param_info(pheT, **kwargs).count
        Y_delt = np.array([])

        for k in range(0, 3):
            idx_k = np.where(snp_vec == k)[0].tolist()
            if len(idx_k) > 0:
                if pheX is not None:
                    pheX_k = pheX[idx_k, :]
                    X = np.tile(np.dot(pheX_k, parin_X[0] if len(parin_X)==1 else parin_X).flatten(order='F'), pheY.shape[1]).reshape((-1, pheY.shape[1]), order='F')
                else:
                    X = np.repeat(0, len(idx_k))

                parin_curve0 = parin_curve[0:n_par_curve]
                if Y_delt.size==0:
                    Y_delt=(pheY[idx_k, :] - obj_curve.get_curve_formula(parin_curve0, pheT[idx_k, :], **kwargs) - X)
                else:
                    Y_delt=np.r_[Y_delt,(pheY[idx_k, :] - obj_curve.get_curve_formula(parin_curve0, pheT[idx_k, :], **kwargs) - X)]
                parin_curve = parin_curve[n_par_curve:len(parin_curve)]

        Y_delt=Y_delt.reshape((-1,X.shape[1]),order='C')
        if sum_type == 0:
            Y_var = np.var(Y_delt, axis=1)
            Y_delt[np.isnan(Y_delt) == True] = 0
            # todo
            pv = np.log(multivariate_normal.pdf(Y_delt, cov=np.diag(Y_delt)))
            A = -np.sum(pv)
        else:
            A = Y_delt[np.isnan(Y_delt) == False]

        return A

    def proc_h1_cov(parin, Y_delt, pheT, obj_covariance, nna_vec, **kwargs):
        mat_cov = obj_covariance.get_matrix(parin, pheT)
        try:
            pv = comm.dmvnorm_fast(Y_delt, np.repeat(0, mat_cov.shape[1]), mat_cov, nna_vec, log=True)
            if pv is None: raise Exception
        except:
            return None

        return - np.sum(pv)

    def add_miss_snp_parameter(parin):
        if geno_type == 3:
            return parin
        else:
            n_par_curve = obj_phe.obj_curve.get_param_info(pheT, **kwargs).count
            par0 = par1 = par2 = np.full(n_par_curve, np.nan)
            if 0 in snp_vec:
                par0 = parin[0:n_par_curve]
                parin = np.delete(parin, range(0, n_par_curve))
            if 1 in snp_vec:
                par1 = parin[0:n_par_curve]
                parin = np.delete(parin, range(0, n_par_curve))
            if 2 in snp_vec:
                par2 = parin[0:n_par_curve]
                parin = np.delete(parin, range(0, n_par_curve))

            result=par0
            result=np.insert(result,result.size,par1)
            result=np.insert(result,result.size,par2)

            return result

    org_snp = snp_vec
    snp_vec = np.choose(snp_vec.astype(int),snp_comb)
    geno_type = len(np.unique(snp_vec))

    if parX is None:
        parin=np.array(obj_phe.est_curve.parX)
        parin=np.insert(parin,len(parin), np.tile(obj_phe.est_curve.param, geno_type))
    else:
        parin = np.tile(obj_phe.est_curve.param, geno_type)

    h1_curve= optim_least_square(parin, proc_h1_curve, pheY=pheY, pheT=pheT, pheX=pheX, snp_vec=snp_vec, parX=parX,
                                  obj_curve=obj_phe.obj_curve, method=(kwargs["opt_method"] == "FAST"), ssr_ref=h0.ssr)

    if (h1_curve is None) or (h1_curve.value is None) or (math.isinf(h1_curve.value)):
        return result.AllInOneVo(x=[add_miss_snp_parameter(parin), obj_phe.est_covar.param], fun=None, R2=None)

    Y_delt = np.full((pheY.shape[0], pheY.shape[1]), np.nan)
    if (pheX is None) or parX is not None:

        parin_curve = add_miss_snp_parameter(h1_curve.x).reshape((3, -1), order='C')
        parin_X = list
        if parX is not None:
            parin_X = parX

    else:

        parin_curve = add_miss_snp_parameter(np.array(h1_curve.x)[pheX.shape[1]:len(h1_curve.x)]).reshape((3, -1),order='C')
        parin_X = np.array(h1_curve.x)[0:pheX.shape[1]]

    for k in range(0, 3):
        idx_k = np.where(snp_vec == k)
        if len(idx_k) > 0:
            if pheX is not None:
                pheX_k = pheX[idx_k, :]
                X = np.tile(np.dot(pheX_k, parin_X).flatten(order='F'), pheY.shape[1]).reshape((-1, pheY.shape[1]), order='F')

            else:
                X = 0

            Y_delt[idx_k,] = pheY[idx_k, :] - obj_phe.obj_curve.get_curve_formula(parin_curve[k, :], pheT[idx_k, :],
                                                                          **kwargs) - X

    parin_covar = np.array(h0.x)[(len(obj_phe.est_curve.parX) + len(obj_phe.est_curve.param)):len(h0.x)]
    mat_cov = obj_phe.obj_covar.get_matrix(parin_covar, pheT)
    # todo 非半正定情况
    if  np.all(np.linalg.eigvals(mat_cov) > 0):
        pv = comm.dmvnorm_fast(Y_delt, np.repeat(0, mat_cov.shape[1]), mat_cov,None, log=True)
    else :
        pv=np.array([0])
    h1_cm = h1_curve
    h1_cm.pv = pv
    h1_cm.value = -np.sum(pv)
    h1_cm.ssr = np.nansum(Y_delt ** 2)

    parin_curve = parin_curve[snp_comb, :]
    if not (0 in org_snp):
        parin_curve[0, :] = np.full(parin_curve.shape[1], np.nan)
    if not (1 in org_snp):
        parin_curve[1, :] = np.full(parin_curve.shape[1], np.nan)
    if not (2 in org_snp):
        parin_curve[2, :] = np.full(parin_curve.shape[1], np.nan)

    h1_cm.x = parin_X
    h1_cm.x=np.insert(h1_cm.x,(h1_cm.x).size,parin_curve.flatten(order='C'))
    h1_cm.x=np.insert(h1_cm.x,(h1_cm.x).size,parin_covar)
    h1_cm.R2 = estimate.get_R2(pd.DataFrame(pheY), pd.DataFrame(pheX), pd.DataFrame(pheT), obj_phe.obj_curve, h1_cm.x, snp_vec)

    if not kwargs["b_permu"]:
        print("H1[F]", h1_cm.value, h1_cm.ssr, "Param=", h1_cm.x, "\n")

    return h1_cm


def proc_fast_h1(obj_phe, pheY, pheX, pheT, snp_vec, h0, **kwargs):
    h1 = proc_fast_h1_comb(obj_phe, pheY, pheX, pheT, snp_vec, h0, np.array([0, 1, 2]), **kwargs)

    pv_diff = h1.pv - (0 if h0.pv is None else h0.pv)

    pv_neg = np.where(
        np.array([np.sum(pv_diff[snp_vec == 0]), np.sum(pv_diff[snp_vec == 1]), np.sum(pv_diff[snp_vec == 2])]) < 0)
    if len(pv_neg) == 0:
        return h1

    h1x = [h1]
    h1x.append(proc_fast_h1_comb(obj_phe, pheY, pheX, pheT, snp_vec, h0, np.array([0, 1, 1]), **kwargs))
    h1x.append(proc_fast_h1_comb(obj_phe, pheY, pheX, pheT, snp_vec, h0, np.array([0, 0, 2]), **kwargs))
    h1x.append(proc_fast_h1_comb(obj_phe, pheY, pheX, pheT, snp_vec, h0, np.array([0, 1, 0]), **kwargs))
    h1x.append(proc_fast_h1_allinone(obj_phe, pheY, pheX, pheT, snp_vec, h0,None, **kwargs))
   # proc_fast_h1_allinone 有两种返回值 后续可优化
    value=[h1.value for h1 in h1x[0:-1]]
    value.append(h1x[-1].fun if h1x[-1].fun is not None else float("inf"))
    idx_min = np.argmin(value)
    print(idx_min)

    return h1x[idx_min]


def optim_least_square(parin, fn_mle, pheY, pheT, pheX, snp_vec=None, parX=None, obj_curve=None, method=None,
                       ssr_ref=None):
    # def f_optim():
    #     parinx = parin * np.random.uniform(0.9, 1.1, len(parin))
    #     control = result.control(maxit=500 + len(parin) * 200, reltol=1e-8)
    #     loop_optim = 1
    #     # to do
    #     while True:
    #         try:
    #             h0 = scipy.optimize.minimize(x0=parinx, func=fn_mle,
    #                                          args=(pheY, pheT, pheX, snp_vec, parX, obj_curve, 0),
    #                                          method="Nelder-Mead" if loop_optim % 2 == 0 else "BFGS",
    #                                          tol=control.reltol, maxiter=control.maxit)
    #             if h0 is None or h0.status == 1:
    #                 raise Exception
    #             break
    #         except:
    #             if h0.status == 1:
    #                 control.maxit = control.maxit * 2
    #                 if control.maxit > 500 * 496:
    #                     control.maxit = 500 * 496
    #                     control.reltol = control.reltol * 2
    #             else:
    #                 print("optim() convergence=", h0.status, " ", h0.message, "\n")
    #
    #             comm.reset_seed()
    #             loop_optim = loop_optim + 1
    #             continue
    #     return h0

    def f_minpack_lm():

        parinx = parin * np.random.uniform(0.9, 1.1, parin.size)
        control = result.control(maxit=500 + parin.size * 100)
        if control.maxit > 4000:
            control.maxit = 4000

        temp_ret =list()
        n_try = 0
        while True:
            h0_nls=OptimizeResult({"value":None})
            try:
                # The purpose of nls.lm is to minimize the sum square of the vector returned by the function fn
                # by a modifiprintion of the Levenberg-Marquardt algorithm.
                # 非线性规划最小二乘
                # todo 损失函数返回值
                h0_nls = scipy.optimize.least_squares(fun=fn_mle, x0=parinx, method="lm",
                                                      args=([pheY, pheT, pheX,snp_vec, parX,
                                                             obj_curve, 1]), max_nfev=control.maxit)

                if len(h0_nls.keys())==1 or (h0_nls is None) or (h0_nls.status < 0):
                    raise Exception
                comm.reset_seed()
                n_try=n_try+1
                # to do nls.lm value 的值 fvec的值 todo
                h0_nls["value"]=np.sum(h0_nls.fun ** 2)
               # h0_nls_value = np.sum(h0_nls.fun ** 2)
                temp_ret.append(h0_nls)
                if ssr_ref is None:
                    break
                else:
                    if h0_nls.value < ssr_ref * ((1.1 + (n_try - 10) * 0.02) if n_try > 10 else 1.1):
                        index = np.argmin([x.value for x in temp_ret])
                        h0_nls = temp_ret[index]
                        break

            except:
                # exceed maximum iteration
                if len(h0_nls.keys())!=1:
                    if h0_nls.status == 0:
                        control.maxit = control.maxit * 2
                        if control.maxit > 4000:
                            control.maxit = 4000
                    else:
                        print("nls.lm(), status=", h0_nls.status, "", h0_nls.message, "\n")
                        parinx = parin * np.random.uniform(0.9, 1.1, len(parin))
                        control.maxit = 500 + len(parin) * 100 + n_try * 50
                        if control.maxit > 4000: control.maxit = 4000
                comm.reset_seed()
                n_try = n_try + 1
                continue
        return h0_nls

    A1 = f_minpack_lm()
    return A1


def proc_full_h0(obj_phe, pheY, pheX, pheT, h0_ref, **kwargs):
    def gradient_h0_cov(parin, pheY, pheT, pheX, obj_curve, obj_covariance, nna_vec, n_par_curve, **kwargs):

        if pheX is not None:

            parin_X = parin[0:pheX.columns.size]
            parin0 = parin[pheX.columns.size, len(parin)]
            parin_curve = parin0[0:n_par_curve]
            parin_covar = parin0[n_par_curve, len(parin0)]

        else:

            parin_X = None
            parin_curve = parin[0:n_par_curve]
            parin_covar = parin[n_par_curve, len(parin)]

        d = np.repeat(0, len(parin))
        mat_cov = obj_covariance.get_matrix(parin_covar, pheT)

        for nna in np.unique(nna_vec):

            nna_idx = np.where(nna_vec == nna)
            col_idx = ~np.isnan(pheY.loc[nna_idx[0], :])

            mat_cov_nna = mat_cov[col_idx, col_idx]
            # todo
            mat_inverse = np.linalg.solve(mat_cov_nna)

            X = np.repeat(0, len(nna_idx))
            if pheX is not None:
                X = np.repeat(np.dot(pheX[nna_idx], parin_X), np.sum(col_idx)).reshape((-1, np.sum(col_idx)), order='F')

            mu = obj_curve.get_curve(parin_curve, pheT.values[nna_idx, col_idx], **kwargs)
            Y_delt = pheY.values[nna_idx, col_idx] - mu - X
            d_u = np.dot(mat_inverse, Y_delt.T)

            d_X = list()
            if pheX is not None:

                for i in range(0, pheX.columns.size):
                    d_X.append(np.nansum(-1 * d_u.T * pheX.values[nna_idx, i]))

            d_curve = list()
            d.list = obj_curve.get_gradient(parin_curve, pheT.values[nna_idx, col_idx], **kwargs)
            for i in range(0, len(parin_curve)):
                d_curve.append(np.nansum(-1 * d.u * d_list[i].T))

            Y_delt[np.isnan(Y_delt) == True] = 0
            d_list = obj_covariance.get_gradient(parin_covar, pheT.values[nna_idx, col_idx])
            d_sigma = -0.5 * (mat_inverse * Y_delt.shape[0] - np.dot(np.dot(mat_inverse, (np.dot(Y_delt, Y_delt))),
                                                                     mat_inverse))
            d_cov = list()
            for i in range(0, len(parin_covar)):
                d_cov.append(np.nansum(d_sigma * d_list[i]))
            d_X.append(d_curve)
            d_X.append(-1 * d_cov)
            d = d + d_X

        return d

    def proc_h0_cov(parin, pheY, pheT, pheX, obj_curve, obj_covariance, nna_vec, n_par_curve, **kwargs):
        if pheX is not None:

            parin_X = parin[0:pheX.shape[1]]
            parin = parin[pheX.shape[1]:len(parin)]
            parin_curve = parin[0:n_par_curve]
            parin_covar = parin[n_par_curve:len(parin)]

        else:

            parin_X = None
            parin_curve = parin[0:n_par_curve]
            parin_covar = parin[n_par_curve:len(parin)]

        X = np.repeat(0, pheY.shape[0])
        if pheX is not None:
            X = np.tile(np.dot(pheX, parin_X).flatten(order='F'), pheY.shape[1]).reshape((-1, pheY.shape[1]), order='F')

        mu = obj_curve.get_curve_formula(parin_curve, pheT, **kwargs)
        Y_delt = pheY - mu - X

        mat_cov = obj_covariance.get_matrix(parin_covar, pheT)

        if (not np.all(np.linalg.eigvals(mat_cov) > 0)):
            return None

        try:

            pv = comm.dmvnorm_fast(Y_delt, np.repeat(0, mat_cov.shape[0]), mat_cov, nna_vec, log=True)
            if np.isnan(pv).any():
                return None
        except:
            return None

        A = np.nansum(pv)
        return -A

    nna_vec = comm.get_non_na_number(pheY)
    n_par_curve = obj_phe.obj_curve.get_param_info(pheT, **kwargs).count

    ## TEST gradient function
    if kwargs["use_gradient"] is not None and kwargs["use_gradient"]:

        parin = np.array(obj_phe.est_curve.parX)
        parin = np.c_[parin, obj_phe.est_curve.param]
        parin = np.c_[parin, obj_phe.est_covar.param]

        r0 = gradient_h0_cov(parin, pheY, pheT, pheX, obj_phe.obj_curve, obj_phe.obj_covar, nna_vec, n_par_curve,
                             **kwargs)
        r1 = [nd.Gradient(proc_h0_cov)(u, pheY, pheT, pheX, obj_phe.obj_curve, obj_phe.obj_covar, nna_vec, n_par_curve,
                                       **kwargs) for u in parin]
        if np.max(abs(r0 - r1)) > 0.1:
            print("Gradient=", r0, "\n")
            print("numDeriv=", r1, "\n")
            sys.exit("The gradient function doesn't work well.\n")

        comm.RW("graident_test", True)

    h0_failure = 0
    h0_best = None
    while h0_failure < 5:
        h0_best = comm.optim_BFGS(obj_phe.est_curve.parX, obj_phe.est_curve.param, obj_phe.est_covar.param,
                                  proc_h0_cov, (gradient_h0_cov if kwargs["use_gradient"] else None),
                                  pheY=pheY, pheT=pheT, pheX=pheX,
                                  obj_curve=obj_phe.obj_curve, obj_covariance=obj_phe.obj_covar,
                                  nna_vec=nna_vec, n_par_curve=n_par_curve, **kwargs)

        if h0_ref is not None and h0_best.value > h0_ref * 1.01:

            h0_failure = h0_failure + 1
            kwargs["min_optim_success"] = kwargs["min_optim_success"] * 2
            h0_best = None
            comm.reset_seed()

        else:
            break

    par = np.array(obj_phe.est_curve.parX)
    par = np.insert(par, len(par),obj_phe.est_curve.param)
    par = np.insert(par, len(par), obj_phe.est_covar.param)
    if h0_best is None or math.isinf(h0_best.fun) or np.isnan(h0_best.fun):
        h0_best = resultVo.AllInOneVo(x=par,fun=None, R2=None)
    else:
        h0_best.R2 = estimate.get_R2(pd.DataFrame(pheY), pd.DataFrame(pheX), pd.DataFrame(pheT), obj_phe.obj_curve, h0_best.x)

    print("H0[F]", h0_best.fun, "Param=", h0_best.x, "\n")
    return h0_best
