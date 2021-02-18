
import random

import numpy as np
import pandas as pd

import FgDmObject as fgObj
import Curve.CurveTool as tool
import Covariance.CovarianceTool as covarTool
import FgDmObject


# '''---------------------
#  fg_simulate
#
#    return fg_obj
#
#       str(fg_obj)
#       _ obj_gen : 'fgwas_gen_obj':
#       _ obj_phe : 'fgwas_phe_obj':
#       _ error   : logi FALSE
# ---------------------'''
def proc_dat_simu(n_obs, n_snp, par_X, par0, par1, par2, par_covar, curve, covar, times, sig_idx, snp_missing,
                  phe_missing):
    obj_gen = fgObj.FgGenObj()
    obj_phe = fgObj.FgPheObj()

    # generate SNPs
    tmp_matrix = proc_simu_geno(n_obs, n_snp, sig_idx, prob_miss=snp_missing)
    tmp_snpinfo = pd.DataFrame(
        {'SNPNAME': np.array(['SNP_' + str(i) for i in range(0, n_snp)]), 'CHR': np.ones(n_snp),
         'POS': np.arange(0, n_snp, 1), 'RefBase': 'A', 'AltBase': 'T'})
    # to do change to get set method
    obj_gen.n_snp = n_snp
    obj_gen.n_ind_total = n_obs

    obj_gen.n_ind_used = n_obs
    obj_gen.reader = FgDmObject.FgDmSimple(type="SIMPLE", description="Simple geno table", n_snp=n_snp,
                                           n_ind_total=n_obs, n_ind_used=n_obs, ids_used=tmp_matrix.columns.values,
                                           file_simple_snp="",
                                           rawData=pd.concat([tmp_snpinfo, tmp_matrix], axis=1),
                                           snpData=fgObj.SnpData(snp_info=tmp_snpinfo, snp_mat=tmp_matrix),
                                           )

    pheX = None

    if len(par_X) > 0:
        data = np.zeros((len(par_X), n_obs))
        for i in range(0, len(par_X)):
            if i == 0:
                data[i] = np.round(np.random.normal(1, 2, n_obs))
            else:
                data[i] = np.random.normal(-1, 1, n_obs)
        pheX = pd.DataFrame(data.T, columns=['X_{0}'.format(i) for i in range(0, len(par_X))],
                            index=['N_{0}'.format(i) for i in range(0, n_obs)])

    # generate traits
    # max_time,min_time
    # options = {"max_times": np.nanmax(times), "min_times": np.nanmin(times)}
    sim_mu = np.zeros((3, len(times)))
    sim_mu[0, :] = curve.get_curve_formula(par0, times, max_time=np.nanmax(times),min_times=np.nanmin(times))
    sim_mu[1, :] = curve.get_curve_formula(par1, times, max_time=np.nanmax(times),min_times=np.nanmin(times))
    sim_mu[2, :] = curve.get_curve_formula(par2, times, max_time=np.nanmax(times),min_times=np.nanmin(times))

    d_g = get_gencode(tmp_matrix.iloc[[sig_idx]], tmp_snpinfo.iloc[[sig_idx]], n_obs)
    # d_g1 = get_gencode(tmp_matrix.iloc[[0]], tmp_snpinfo.iloc[[0]], n_obs)
    # d_g2 = get_gencode(tmp_matrix.iloc[[1]], tmp_snpinfo.iloc[[1]], n_obs)
    sim_covar = covar.get_matrix(par_covar, times)
    pheY_data = np.zeros((n_obs, len(times)))
    for i in range(0, n_obs):
        if d_g[:, i] == 9: d_g[:, i] = np.round(np.random.uniform(0, 2, 1))
        pheY_data[i] = np.random.multivariate_normal(sim_mu[d_g[:, i], :][0], sim_covar, 1)
        if pheX is not None:
            pheY_data[i] = pheY_data[i] + ((pheX.iloc[[i]] * par_X).sum(axis=1))[0]
    pheY = pd.DataFrame(pheY_data, columns=['Y_{0}'.format(i) for i in times],
                        index=['N_{0}'.format(i) for i in range(0, n_obs)])
    obj_phe.obj_curve = curve
    obj_phe.obj_covar = covar
    obj_phe.ids = pheY._stat_axis.values.tolist()
    obj_phe.pheY = pheY
    obj_phe.pheX = pheX

    columns = list('T_{0}'.format(i) for i in range(0, len(times)))
    # 不加ID 列名
    data = np.repeat(times, n_obs).reshape(len(times), n_obs)
    # print(pheY._stat_axis.values.tolist())
    obj_phe.pheT = pd.DataFrame(data.T,
                                columns=columns,
                                index=pheY._stat_axis.values.tolist())

    obj = fgObj.SimuObj(obj_gen, obj_phe, False)

    print(" Data simulation is done![Sig={0}]\n".format(sig_idx))

    return obj


def write_csv(fg_obj):
    fg_obj.obj_phe.pheY.to_csv(fg_obj.obj_phe.file_pheY, sep=' ', header=True, index=True)
    fg_obj.obj_phe.pheT.to_csv(fg_obj.obj_phe.file_pheT, sep=' ', header=True, index=True)
    if fg_obj.obj_phe.file_pheX is not None:
        fg_obj.obj_phe.pheX.to_csv(fg_obj.obj_phe.file_pheX, sep=' ', header=True, index=True)
    else:
        fg_obj.obj_phe.file_pheX = None


# time_points type np.ndarray
def fg_simulate(curveType, covarianceType, n_obs, n_snp, time_points, par0=None, par1=None, par2=None, par_covar=None,
                par_X=None, phe_missing=0.03, snp_missing=0.03, sig_pos=None, plink_format=False, file_prefix=None):
    if not isinstance(time_points, np.ndarray):
        time_points = np.array(time_points)
    curve = tool.get_curve(curveType)
    covar = covarTool.getCovariance(covarianceType)

    curve_simu_param = curve.get_simu_param(time_points)
    covar_simu_param = covar.get_simu_param(time_points)

    par0 = curve_simu_param[0]
    par1 = curve_simu_param[1]
    par2 = curve_simu_param[2]

    par_covar = covar_simu_param

    fg_obj = proc_dat_simu(n_obs, n_snp, par_X, par0, par1, par2, par_covar, curve, covar, time_points,
                           sig_pos,
                           snp_missing, phe_missing)

    if file_prefix is not None:
        fg_obj.obj_phe.file_pheX = file_prefix + "_pheX.csv"
        fg_obj.obj_phe.file_pheY = file_prefix + "_pheY.csv"
        fg_obj.obj_phe.file_pheT = file_prefix + "_pheT.csv"

    write_csv(fg_obj)

    if not plink_format:
        file_gen_dat = file_prefix + ".geno.tab"
        tb_gen = pd.concat([fg_obj.obj_gen.reader.get_snpinfo(None),
                            fg_obj.obj_gen.reader.get_snpmat(None, impute=False, allel=True).snpmat], axis=1)
        # todo fg_obj.obj_phe.file_pheX is str, rownames(str)=NULL
        # tb_gen.columns =tb_gen.columns.values,fg_obj.obj_phe.file_pheX.index.values)
        tb_gen.to_csv(file_gen_dat, sep="\t")

        fg_obj.obj_gen.files = file_gen_dat
    else:
        snp_mat = fg_obj.obj_gen.reader.get_snpmat(impute=False, allel=False).snpmat
        snp_info = fg_obj.obj_gen.reader.get_snpinfo()
        # todo
        # r =convert_simpe_to_plink(data_frame(snp_info[, c(2, 1)], 0, snp_info[, c(3:5)]), snp_mat, paste(file_prefix,
        #                                                                                                     "_geno",
        #                                                                                                     sep="") );

        # fg_obj_obj_gen_files =(file_plink_bed = r_file_plink_bed,
        #                        file_plink_bim = r_file_plink_bim,
        #                                           file_plink_fam = r_file_plink_fam);
    return fg_obj


# class FgDmSimple(object):
#     def __init__(self, type, description, n_snp, n_ind_total, n_ind_used, ids_used, file_simple_snp, rawdata,
#                  snpdata):
#         self.type = type
#         self.description = description
#         self.n_snp = n_snp
#         self.n_ind_total = n_ind_total
#         self.n_ind_used = n_ind_used
#         self.ids_used = ids_used
#         self.file_simple_snp = file_simple_snp
#         self.rawdata = rawdata
#         self.snpdata = snpdata


def get_gencode(gen, snp_info, n_obs):
    d_g = np.full((1, n_obs), 9)
    gen2 = gen.values.astype(str)
    tmp = snp_info.values
    snpB = tmp[0, 4]
    snpA = tmp[0, 3]

    QQ2 = snpB + snpB
    qq0 = snpA + snpA
    Qq1 = (snpA + snpB, snpB + snpA)

    d_g[gen2 == QQ2] = 2
    d_g[gen2 == qq0] = 0
    d_g[gen2 == Qq1[0]] = 1
    d_g[gen2 == Qq1[1]] = 1

    return d_g


# aim 检查snp_mat 与snp_info 行数不同
def generate_bc_marker(n_obs, dist):
    if dist[0] != 0:
        cm = np.hstack((0, dist)) / 100
    else:
        cm = dist / 100

    n = len(cm)
    rs = 1 / 2 * (np.exp(2 * cm) - np.exp(-2 * cm)) / (np.exp(2 * cm) + np.exp(-2 * cm))

    mk = np.zeros((n_obs, n))

    for j in range(0, n_obs):
        mk[j, 0] = np.random.uniform(0, 1, 1) > 0.5

    for i in range(1, n):
        for j in range(1, n_obs):
            if mk[j, i - 1] == 1:
                mk[j, i] = np.random.uniform(0, 1, 1) > rs[i]
            else:
                mk[j, i] = np.random.uniform(0, 1, 1) < rs[i]

    return mk


def proc_simu_geno(n_obs, n_snp, sig_idx, prob_miss=0.03):
    # runif取500个值， cumsum 累加
    dist = np.random.uniform(0.05, 0.12, n_snp).cumsum()
    snp1 = generate_bc_marker(n_obs, dist).T
    snp2 = generate_bc_marker(n_obs, dist).T
    for i in range(0, n_snp):
        n_miss = random.uniform(0, n_obs * prob_miss)
        if n_miss >= 1:
            snp1[i, random.sample(range(0, n_obs), n_obs)[0:round(n_miss)]] = 9

    for i in range(0, n_snp):
        n_miss = random.uniform(0, n_obs * prob_miss)
        if n_miss >= 1:
            snp2[i, random.sample(range(0, n_obs), n_obs)[0:round(n_miss)]] = 9

    cors = np.zeros(n_snp)
    snpx = snp1 + snp2
    for i in range(0, n_snp):
        cors[i] = np.corrcoef(snpx[i].tolist(), snpx[sig_idx].tolist())[0, 1]
    cor_high = np.where(abs(cors) > 0.75)
    if len(cor_high) >= 2:
        for i in range(0, len(cor_high)):
            if cor_high[i] != sig_idx:
                snp1[cor_high[i]] = snp1[cor_high[i] - 1]
                snp2[cor_high[i]] = snp2[cor_high[i] - 1]
    snp1_s = snp1.astype(str).flatten('F')
    snp2_s = snp2.astype(str).flatten('F')

    if len(np.where(snp1_s == "0.0")) > 0:
        snp1_s[np.where(snp1_s == "0.0")] = 'A'
    if len(np.where(snp1_s == "1.0")) > 0:
        snp1_s[np.where(snp1_s == "1.0")] = "T"
    if len(np.where(snp1_s == "9.0")) > 0:
        snp1_s[np.where(snp1_s == "9.0")] = "."

    if len(np.where(snp2_s == "0.0")) > 0:
        snp2_s[np.where(snp2_s == "0.0")] = "A"
    if len(np.where(snp2_s == "1.0")) > 0:
        snp2_s[np.where(snp2_s == "1.0")] = "T"
    if len(np.where(snp2_s == "9.0")) > 0:
        snp2_s[np.where(snp2_s == "9.0")] = "."
    #
    snp_s = np.char.add(snp1_s, snp2_s)[0:n_snp * n_obs]
    gen = pd.DataFrame(np.reshape(snp_s[0:n_snp * n_obs], (n_snp, n_obs), 'F'),
                       columns=['N_{0}'.format(i) for i in range(0, n_obs)],
                       index=['SNP_{0}'.format(i) for i in range(0, n_snp)])
    return gen
