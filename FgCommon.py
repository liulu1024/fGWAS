import random
import sys

import numpy as np
import pandas as pd

from scipy.stats import multivariate_normal
import FgDmObject as fgObj
import ResultVo as result
from scipy.optimize import OptimizeResult
from scipy.optimize import minimize as minimize

# common
N_SEED = 1
DEBUG = False
TRY_SILENT = True
MU_RANGE_LOOP = 1
RW = {"n_seed": 0}


def reset_seed(**kwargs):
    RW["n_seed"] = RW["n_seed"] + random.uniform(1, 50)
    random.seed(N_SEED)
    if DEBUG: print(".")


def paste(start,list, seq):
    if not isinstance(list,np.ndarray):
        list=np.array(list)
    return  [start+seq+i for i in list.astype(str)]


def get_non_na_number(y_resd):
    nna_mat = np.ones((y_resd.shape[0], y_resd.shape[1]))
    nna_mat[np.where(np.isnan(y_resd) == True)] = 0

    nna_vec = np.repeat(0, nna_mat.shape[0])
    for i in range(0, nna_mat.shape[1]):
        nna_vec = nna_vec * 2 + nna_mat[:, i]

    return nna_vec


def dmvnorm_fast(y_resd, mu, cov_mat, nna_vec, log=True):
    if nna_vec is None:
        nna_vec = get_non_na_number(y_resd)
    pv = np.repeat(np.nan, nna_vec.shape[0])

    for nna in np.unique(nna_vec):
        if nna > 0:
            nna_idx = np.where(nna_vec == nna)
            col_idx = ~ np.isnan(y_resd[nna_idx[0][0:1], :])
            col_idx = col_idx[0, :]
            var = multivariate_normal(mean=mu[col_idx], cov=cov_mat[col_idx, :][:, col_idx])
            pv[nna_idx] = np.log(var.pdf(y_resd[nna_idx[0], :][:, col_idx])) if log else np.log(var.pdf(y_resd[nna_idx[0], :][:, col_idx]))

    return pv


def impute_simple_snp(snpmat):
    def f_imputed(snp):
        s_miss = np.where(snp == np.nan)
        if len(s_miss) > 0:
            n_AA = len(np.where(snp == 0))
            n_Aa = len(np.where(snp == 1))
            n_aa = len(np.where(snp == 2))
            n_s = n_AA + n_Aa + n_aa

            r_miss = np.random.uniform(size=len(s_miss))
            r_snp = np.repeat(2, len(s_miss))
            r_snp[r_miss <= n_AA / n_s] = 0
            r_snp[r_miss <= (n_AA + n_Aa) / n_s] < -1
            snp[s_miss] = r_snp

        if np.mean(snp) / 2 > 0.5:
            snp = 2 - snp

        return fgObj.Snp(snp=snp, NMISS=len(s_miss), MAF=np.mean(snp) / 2)

    r_imp = []
    total_miss = len(np.where(snpmat == np.nan))

    for i in snpmat[:, 0:snpmat.shape[1]]:
        r_imp.append(f_imputed(snpmat[:, i]))
    data = []
    for i in range(0, snpmat.columns.size):
        data.appent(pd.DataFrame(r_imp[[i]].snp))

    snpmat = pd.DataFrame(data, axis=1)
    NMISS = [r_imp[[i]].NMISS for i in range(0, snpmat.columns.size)]
    MAF = [r_imp[[i]].MAF for i in range(0, snpmat.columns.size)]
    print("* Missing SNPs are imputed(", total_miss, "SNPs).\n")

    return fgObj.SnpMat(snpmat=snpmat, NMISS=NMISS, MAF=MAF)


def optim_BFGS(par_x, par_curve, par_covar, proc_fn, proc_gr, **kwargs):
    mle_control = {"optim_success": 0, "optim_loop": 0}
    parinx = np.array(par_x)
    parinx = np.insert(parinx, len(parinx), par_curve)
    parinx = np.insert(parinx, len(parinx), par_covar)
    n_par = len(parinx)
    control = result.control(maxit=500 + n_par * 200, reltol=1e-2)
    h0_best = OptimizeResult({"fun": float("Inf"), "x": np.repeat(np.nan, len(parinx))})
    h0 = OptimizeResult()

    # prepare param
    pheY = kwargs["pheY"] if "pheY" in kwargs.keys() else None
    pheT = kwargs["pheT"] if "pheT" in kwargs.keys() else None
    pheX = kwargs["pheX"] if "pheX" in kwargs.keys() else None
    Y_delt = kwargs["Y_delt"] if "Y_delt" in kwargs.keys() else None
    funType = kwargs["funType"] if "funType" in kwargs.keys() else None
    obj_curve = kwargs["obj_curve"] if "obj_curve" in kwargs.keys() else None
    obj_covariance = kwargs["obj_covariance"] if "obj_covariance" in kwargs.keys() else None
    nna_vec = kwargs["nna_vec"] if "nna_vec" in kwargs.keys() else None
    n_par_curve = kwargs["n_par_curve"] if "n_par_curve" in kwargs.keys() else None

    while mle_control["optim_success"] < kwargs["min_optim_success"] and mle_control["optim_loop"] < kwargs[
        "max_optim_failure"]:
        alter_method = "Nelder-Mead" if proc_gr is None else "CG"
        try:
            # method="BFGS" if mle_control["optim_loop"] % 2 == 0 else alter_method,
            # 对应关系 x:par value:fun
            #papare args
            # todo 对于不同的参数需要输入不同的变量值 目前是硬编码
            if funType=="matrix":args1=(Y_delt, pheT, obj_covariance, nna_vec)
            else: args1=((pheY, pheT, pheX, obj_curve, obj_covariance, nna_vec, n_par_curve))
            h0 = minimize(fun=proc_fn, x0=parinx,
                          args=(args1), jac=proc_gr,
                          method="Nelder-Mead",options={"maxiter": control.maxit}, tol=control.reltol)
            TRY_SILENT = False
            mle_control["optim_loop"] = mle_control["optim_loop"] + 1
            if len(h0.keys()) == 0 or (not h0.success):
                raise Exception
            if h0.fun < h0_best.fun: h0_best = h0

            mle_control["optim_success"] = mle_control["optim_success"] + 1

            reset_seed()

            ## avoid overflow, fine tuning
            parinx = np.array(par_x)
            parinx = np.insert(parinx,len(parinx),par_curve * np.random.uniform(0.99, 1.0, len(par_curve)))
            parinx = np.insert(parinx, len(parinx), par_covar)
        except:
            mle_control["optim_loop"] = mle_control["optim_loop"] + 1
            if len(h0.keys()) != 0:
                if h0.status == 2:
                    control.maxit = control.maxit * 2
                    if control.maxit > 500 * 512:
                        control.maxit = 500 * 512
                else:
                    print("optim, convergence=", h0.status, "", h0.message, "\n")

            reset_seed()
            parinx = np.array(par_x * np.random.uniform(0.95, 1.05, len(par_x)))
            parinx = np.insert(parinx, len(parinx), par_curve)
            parinx = np.insert(parinx, len(parinx), par_covar)

            continue

    return h0_best


'''====pheno=='''


def select_individuals(object, ids_used):
    idx_match = np.array(ids_used) == np.array(object.ids)
    if any(idx_match == False):
        sys.exit("Some IDs are not matached in the SNP data file")

    object.n_ind = len(ids_used)
    object.ids = ids_used
    object.pheY = object.pheY.iloc[idx_match, :]
    object.pheX = object.pheX.iloc[idx_match, :]
    # object.pheZ = object.pheZ.iloc[idx_match, :]

    return object
