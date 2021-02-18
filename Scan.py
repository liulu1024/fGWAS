import sys

import numpy as np
import  pandas as pd

import Curve as curve
import Covariance as covar
import math
from FgCommon import paste
from  FgCommon import select_individuals
import FgDmObject as fgObj
import  FgCore as core
from Estimate import fg_dat_est
import time

from multiprocessing import Pool
import ResultVo as result


class ScanObj(object):
    def __init__(self, ret_fast=None, obj_gen=None, obj_phe=None, systime_time=None, options=None):
        self.ret_fast = ret_fast
        self.obj_gen = obj_gen
        self.obj_phe = obj_phe
        self.system_time = systime_time
        self.options = options


def check_ids(obj_gen, obj_phe):
    if any(obj_phe.pheT.index.values != obj_phe.pheX.index.values) and any(obj_phe.pheX.index.values != obj_phe.pheY.index.values):
        sys.exit()
    ids_phe = obj_phe.pheY.index.values
    ids_gen = obj_gen.reader.get_individuals()

    m_phe = ids_phe == ids_gen
    if not np.size(np.where(m_phe == False)) == 0:
        print("!", np.size(np.where(m_phe == False)), "IDs, can not be found in the genotype file.\n")
        print("! First 5 IDs:\n")
        print(ids_gen[np.where(m_phe == False)][0:5])

    m_phe = ids_phe == ids_gen
    if np.size(np.where(m_phe==False)) > 0:
        print("!", len(np.where(m_phe == False)), "IDs, can not be found in the phenotype file.\n")
        print("! First 5 IDs:\n")
        print(ids_phe[np.where(m_phe == False)][0:5])
    # 交集
    ids_com = [i for i in ids_gen if i in ids_phe]
    return ids_com


# 基因型和表现型
def fg_snpscan(obj_gen, obj_phe, method, curve_type=None, covariance_type=None, snp_sub=None, permutation=None,
               **kwargs):
    if not isinstance(obj_gen, fgObj.FgGenObj):
        sys.exit("The first paramater should be genotype object.")
    if not isinstance(obj_phe, fgObj.FgPheObj):
        sys.exit("The second paramater should be phenotype object.")

    default_kwargs = {"max_optim_failure": 2, "min_optim_success": 2, "intercept": False, "degree": 3, "ncores": 1,
                      "verbose": False, "use_gradient": False, "piecewise": 1000, "use_snowfall": True}
    default_kwargs.update(kwargs)
    kwargs = default_kwargs
    kwargs["opt_method"] = method.upper()

    if kwargs["verbose"]:
        print("[ SNP Scanning ] \n")

    ids_com = check_ids(obj_gen, obj_phe)
    # 选取样本
    obj_phe = select_individuals(obj_phe,ids_com)
    obj_gen.reader.select_individuals(ids_com)

    if curve_type is None and issubclass(obj_phe.obj_curve, curve.BaseCurve):
        curve_type = obj_phe.obj_curve.curveType
    if curve_type is None and issubclass(obj_phe.obj_covar, covar.BaseCovariance):
        covariance_type = obj_phe.obj_covar.covarType
    if curve_type.upper() != obj_phe.obj_curve.curveType.upper() or covariance_type.upper() != obj_phe.obj_covar.covarType:
        # test obj_phe = fg_dat_est(obj_phe, curve_type, covariance_type, "testDraw.pdf", **kwargs)
        obj_phe = fg_dat_est(obj_phe, curve_type, covariance_type, None, **kwargs)

    if kwargs["verbose"]:
        print("Genetic Effect Analysis by '", method, "' method......\n", sep="")

    ret = {"error": False}
    if snp_sub is None: snp_sub = range(0, obj_gen.n_snp)

    if isinstance(snp_sub, str):
        snp_sub = obj_gen.reader.get_snpindex(snp_sub)
        snp_sub = snp_sub[np.where(np.isnan(snp_sub)==False)]

    if len(snp_sub) == 0:
        sys.exit("No SNPs are found in the genotype file.")

    r_time = None

    methods = ["FAST", "FAST-NORM", "FGWAS", "OPTIM-FGWAS"]
    if method.upper() in methods:
        if obj_phe.est_curve is None:
            curveType = obj_phe.obj_curve.curveType if obj_phe.obj_curve.curveType is not None else "auto"
            covarType = obj_phe.obj_covar.covarType if obj_phe.obj_covar.covarType is not None else "auto"
            obj_phe = fg_dat_est(obj_phe, curveType, covarType, None, **kwargs)
        obj_phe.h0=core.proc_est_h0(obj_phe,**kwargs)

    ret = ScanObj()

    if method.upper() in ["FAST", "FAST-NORM"]:
        start_time = time.time()
        r = snpsnp_curve(obj_gen, obj_phe, snp_sub, **kwargs)
        end_time = time.time()
        r_time = end_time - start_time
        if r.error:
            sys.exit()(r.err_info)

        ret.ret_fast = r

    ret.obj_gen = obj_gen
    ret.obj_phe = obj_phe
    ret.system_time = r_time
    ret.options = kwargs

    return ret


def snpsnp_curve(obj_gen, obj_phe, snp_idx=None, **kwargs):
    #test
    kwargs["ncores"]=1
    if snp_idx is None: snp_idx = range(0, obj_gen.n_snp)

    snp_len = len(snp_idx)
    snp_sect0 = np.arange(0, snp_len, kwargs["piecewise"])
    snp_sect1 = [snp_sect0[-1] - 1, snp_len]

    intercept = False if obj_phe.intercept is None else obj_phe.intercept
    r_list = list()

    for i in range(0, len(snp_sect0)):

        print("  Calculated SNP Range =", snp_idx[snp_sect0[i]], snp_idx[snp_sect1[i]], "\n")

        # sub,snpmat snpinfo
        snp_idx_sub = snp_idx[snp_sect0[i]:snp_sect1[i]]
        snp_mat = obj_gen.reader.get_snpmat(snp_idx_sub)
        snp_info = obj_gen.reader.get_snpinfo(snp_idx_sub)

        snps_percore = math.ceil(len(snp_idx_sub) / kwargs["ncores"])
        used_ncores = math.ceil(len(snp_idx_sub) / snps_percore)

        def cpu_fun(k):
            sub_rows = len(snp_idx_sub)
            i_start = sub_rows if (k * snps_percore + 1 > sub_rows) else k * snps_percore + 1
            i_stop = sub_rows if (k+1) * snps_percore > sub_rows else (k+1) * snps_percore
            r_mle = list()
            for i in range(i_start, i_stop):
                r = None
                #try:
                r = core.proc_mle(snp_idx_sub[i], np.array(snp_info.loc["SNP_{0}".format(i)]), np.array(snp_mat.snpmat.loc[:,"SNP_{0}".format(i)]), obj_phe, intercept, **kwargs)
                # except:
                #     continue

                r_mle.append(r)

            return r_mle

        # and  require(snowfall) ) 并行处理


        if kwargs["ncores"] > 1:
            # todo 并行处理未完善
            if kwargs["verbose"]: print("Starting parallel computing......\n")

            #pool = Pool(kwargs["ncores"])
            pool = Pool(2)

            r_cluster = pool.apply_async(cpu_fun, range(0, used_ncores))
            if kwargs["verbose"]: print("Stopping parallel computing(snowfall/snow)\n")
        else:
            if kwargs["verbose"]:  print("Starting piecewise analysis(parallel package)......\n")
            # apply 的并行版本 ??
            r_cluster = [cpu_fun(i) for i in range(0, used_ncores)]

            if kwargs["verbose"]: print("Stopping piecewise.\n")
        r_list.append(r_cluster)
    # 拼接为一个矩阵
    r_fgwas=None
    for i in range(0,len(r_list)):
        # 含n个DataFrame的list
        cluster=r_list[i]
        for j in range(0,len(cluster)):
            df_1 = cluster[i]
            for z in range(0,len(df_1)):
              df=df_1[i]
              if r_fgwas is None:
                 r_fgwas=df
              else:
                 #按列拼接
                 r_fgwas=pd.concat([r_fgwas,df])


    curve_par_names = obj_phe.obj_curve.get_param_info(obj_phe.pheT).name
    covar_pars_names = obj_phe.obj_covar.get_param_info(obj_phe.pheT).name
    par_curve_len = obj_phe.obj_curve.get_param_info(obj_phe.pheT).count
    par_covar_len = obj_phe.obj_covar.get_param_info(obj_phe.pheT).count

    re_names = r_fgwas.columns.values.astype(str)
    re_names[0:4] = ["INDEX", "NAME", "CHR", "POS"]
    re_names[4:13] = ["MAF", "NMISS", "SNP0", "SNP1", "SNP2", "GENO", "LR", "pv", "L0"]

    if intercept:
        par_X_len = 1 if obj_phe.pheX is None else 1 + obj_phe.pheX.columns.size
    else:
        par_X_len = 0 if obj_phe.pheX is None else obj_phe.pheX.columns.size

    idx_start = 13

    if par_X_len > 0:
        num=np.arange(1,par_X_len)- (1 if intercept else 0)
        re_names[idx_start:(idx_start + par_X_len)] = paste("h0.X", num,"")
        idx_start = idx_start + par_X_len


    re_names[idx_start:(idx_start + par_curve_len)] = paste("h0", curve_par_names,".")
    idx_start = idx_start + par_curve_len
    re_names[idx_start:(idx_start + par_covar_len)] = paste("h0X", covar_pars_names,".")
    idx_start = idx_start + par_covar_len
    re_names[idx_start] = "L1"
    idx_start = idx_start + 1
    if par_X_len > 0:
        re_names[idx_start:(idx_start + par_X_len)] = paste("h1.X", np.arange(1, par_X_len) - (1 if intercept else 0),
                                                                "")
        idx_start = idx_start + par_X_len

    re_names[idx_start:(idx_start + par_curve_len)] = paste("h1.G0", curve_par_names, ".")
    idx_start = idx_start + par_curve_len
    re_names[idx_start:(idx_start + par_curve_len)] = paste("h1.G1", curve_par_names, ".")
    idx_start = idx_start + par_curve_len
    re_names[idx_start:(idx_start + par_curve_len)] = paste("h1.G2", curve_par_names, ".")
    idx_start = idx_start + par_curve_len
    re_names[idx_start:(idx_start + par_covar_len)] = paste("h1X", covar_pars_names, ".")

    re_names[idx_start + par_covar_len] = "h0.R2"
    re_names[idx_start + par_covar_len] = "h1.R2"

    r_fgwas.columns = re_names

    return result.SnpCurveVo(False, r_fgwas)


def snp_cluster(obj_gen, snp_idx=None, dist=3):
    if snp_idx is None: snp_idx = range(0, obj_gen.n_snp)
    snp_len = len(snp_idx)
    snp_sect0 = np.arange(1, snp_len, 1000)
    snp_sect1 = [snp_sect0[-1] - 1, snp_len]

    snp_idx.sub = snp_idx[snp_sect0[0]:snp_sect1[0]]
    snp_mat = obj_gen.reader.get_snpmat(snp_idx.sub)

    r_cluster = list()
    r_cluster[0] = {"snp0": snp_mat.snpmat[:, 0], "grp": [snp_idx[snp_sect0[0]]]}

    for i in range(0, len(snp_sect0)):

        # if(kwargs.verbose)
        print("  Calculated SNP Range =", snp_idx[snp_sect0[i]], snp_idx[snp_sect1[i]], "\n")

        snp_idx.sub = snp_idx[snp_sect0[i]:snp_sect1[i]]
        snp_mat = obj_gen.reader.get_snpmat(snp_idx.sub)

        def get_min_dist(snp0):
            def min_dist_a0a1a2(a0, a1, a2):
                snpx = snp0
                snpx[np.where(snp0 == 0)] = a0
                snpx[np.where(snp0 == 1)] = a1
                snpx[np.where(snp0 == 2)] = a2

                dists = [np.where(
                    snpx - r_cluster[i].snp0[(snpx - r_cluster[i].snp0 != 0) & (snpx - r_cluster[i].snp0 != np.nan)])
                         for i in range(0, len(r_cluster))]

                return result.MinDist(np.where(dists == np.min(dists)), dist=np.nanmin(dists))

            d1 = min_dist_a0a1a2(0, 1, 2)
            d2 = min_dist_a0a1a2(0, 2, 1)
            d3 = min_dist_a0a1a2(1, 0, 2)
            d4 = min_dist_a0a1a2(1, 2, 0)
            d5 = min_dist_a0a1a2(2, 0, 1)
            d6 = min_dist_a0a1a2(2, 1, 0)

            min_grp = [d1.min_grp, d2.min_grp, d3.min_grp, d4.min_grp, d5.min_grp, d6.min_grp]
            min_dist = [d1.dist, d2.dist, d3.dist, d4.dist, d5.dist, d6.dist]

            idx = np.where(min_dist == np.min(min_dist))

            return result.MinDist(min_grp=min_grp[idx], dist=min.dist[idx])

        for k in range(0, snp_idx.sub.index.size):
            min_dist = get_min_dist(snp_mat.snpmat[:, k])
            print(k, min_dist.dist, min_dist.min_grp, "\n")
            if min_dist.dist <= dist:
                r_cluster[[min_dist.min_grp]].grp = [r_cluster[min_dist.min_grp].grp, snp_idx.sub[k]]
            else:
                r_cluster[[len(r_cluster) + 1]] = {"snp0": snp_mat.snpmat[:, k], "grp": [snp_idx.sub[k]]}

    return r_cluster
