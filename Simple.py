import os
import sys
import pandas as pd
import FgDmObject as fgObj
import numpy as np
from numpy import NaN
import FgCommon as comm


def create_simple_obj(file_simple_snp, **args):
    if file_simple_snp is None:
        os.exit("! file_simple_snp  must be assigned with the valid values.")

    objref = fgObj.FgDmSimple(file_simple_snp, pd.DataFrame(), list())

    tb_gen = pd.read_table(objref.file_simple_snp, header=0)
    objref.n_snp = tb_gen.index.size
    objref.n_ind_total = tb_gen.columns.size - 5
    objref.n_ind_used = tb_gen.columns.size - 5
    objref.ids_used = tb_gen.columns.values[5:tb_gen.columns.size]

    load_simple_snp(objref, args["verbose"])

    return objref


def simple_shrink(objref):
    return objref


def simple_get_snpmat(objref, snp_idx=None, impute=False, allel=False):
    if snp_idx is None:
        allel_mat = objref.snpData.snp_mat[objref.ids_used]
        allel_info = objref.snpData.snp_info

    else:
        allel_mat = objref.snpData.snp_mat.loc[['SNP_{0}'.format(i) for i in snp_idx], objref.ids_used]
        allel_info = objref.snpData.snp_info.loc[['SNP_{0}'.format(i) for i in snp_idx]]

    ##BB/AB/AA/A./B.==>2/1/0/NA
    def covert_geno(allel_info, allel_geno):
        snpB = allel_info[[4]].values
        snpA = allel_info[[3]].values
        QQ2 =snpB+snpB
        qq0 =snpA+snpA
        Qq1 =[snpA+snpB,snpB+snpA]

        d_g = np.full((len(allel_geno), 1),np.nan)
        d_g[allel_geno == QQ2[0]] = 2
        d_g[allel_geno == qq0[0]] = 0
        d_g[allel_geno == Qq1[0][0]] = 1
        d_g[allel_geno == Qq1[1][0]] = 1

        if np.nanmean(d_g) / 2 > 0.5:
            d_g = 2 - d_g

        return d_g

        ## if data element is AA/AB/BB

    # AT AA Category
    if isinstance(allel_mat.iloc[0, 0], str):
        mat = []
        for i in range(0, allel_mat.index.size):
            mat.append(pd.DataFrame(covert_geno(allel_info.iloc[i, :], (allel_mat.iloc[i, :]))))

        snpmat = pd.concat(mat, axis=1)

    ## if data element is 0/1/2/NA
    else:
        # 转置
        snpmat = allel_mat.swapaxes(0, 1)

    snpmat.index = allel_mat.columns.values
    snpmat.columns = allel_mat.index.values

    if not allel:
        if impute:
            snpmat_imp = comm.impute_simple_snp(snpmat)
            return fgObj.SnpMat(snpmat=snpmat_imp.snpmat, NMISS=snpmat_imp.NMISS, MAF=snpmat_imp.MAF)
        else:
            NMISS = []
            for i in range(0,snpmat.index.size):
                snp=snpmat.iloc[i, :]
                NMISS.append(np.isnan(snp).sum())
            return fgObj.SnpMat(snpmat=snpmat, NMISS=NMISS, MAF=np.nanmean(snpmat.values / 2, axis=0))
    else:
        return fgObj.SnpMat(snpmat=allel_mat, NMISS=None, MAF=None)


def simple_get_snpinfo(objref, snp_idx=None):
    if snp_idx is None:
        return objref.snpData.snp_info
    else:
        return objref.snpData.snp_info.loc[['SNP_{0}'.format(i) for i in snp_idx]]


def simple_get_individuals(objref):
    return objref.snpData.snp_mat.columns.values


def simple_select_individuals(objref, ids_used):
    objref.ids_used = ids_used

    load_simple_snp(objref,True)

    ids_all = objref.snpData.snp_mat.columns.values
    if not (ids_used == ids_all).all:
        sys.exit("Some IDs are not matached in the SNP data file")

    objref.n_ind_used = len(ids_used)
    return objref


def simple_get_used_individuals(objref):
    return objref.ids_used


def simple_get_snpindex(objref, snp_names):
    return np.where((snp_names == objref.snpData.snp_info.index.values) == True)



def fg_load_simple(file_simple_snp=None, **args):
    if file_simple_snp is None:
        sys.exit("! file_simple_snp must be assigned with the valid file name..")

    if args["verbose"]:
        print("[ Loading Simple ] \n")
        print(" Checking the parameters ......\n")
        print("* SIMPLE SNP File = %s\n" % file_simple_snp)

    if not file_simple_snp.exists():
        sys.exit("Failed to open Simple data files.")

    params = {"file_simple_snp": file_simple_snp}
    obj_gen = create_simple_obj(file_simple_snp, args)
    ret = fgObj.FgGenObj(obj_gen.n_snp, obj_gen.n_ind_total, obj_gen.n_ind_used, obj_gen, None, args, params)

    return ret


# TODO plink command
def plink_command(plink_path, plink_parms):
    # t1    = try(system(paste(c(plink.path, plink.parms), collapse=" "), intern = TRUE))
    # ## show(t1)
    #
    # return(t1)
    pass


# get Identity-By-State matrix using PLINK command
def fg_simple_getPCA(objref, plink_path):
    snp_mat = objref.reader.get_snpmat(None, impute=False, allel=False).snpmat
    snp_info = objref.reader.get_snpinfo(None)
    snp_info.columns.values = ["SNP", "CHR", "POS", "A1", "A2"]


# TODO 因子
# snp_mat    = snp_mat[, with(snp.info, order(CHR, POS))]
# snp_info    = snp_info[ with(snp.info, order(CHR, POS)),]
#
# #change chromosome string to number
# snp.info.CHR=as.numeric(factor(snp.info.CHR))
#
# snp.file.base    = tempfile(fileext = ".plink")
# r    = convert_simpe_to_plink( data.frame(snp.info[,c(2,1)], 0, snp.info[,c(3:5)]),  snp.mat, snp.file.base)
#
# plink.out.pca    = paste(snp.file.base, "pca", sep=".")
# t0    = plink_command( plink.path, c ( "--bfile ", snp.file.base, "--pca --out ", plink.out.pca)  )
#
# tb    = try(read.table(paste(plink.out.pca, "eigenvec", sep=".")))
# if (class(tb)=="try-error")
# R
#    show(t0)
#    sys.exit("Failed to call PLINK.")
#
#
# unlink( paste(snp.file.base, c("bim", "bed", "fam", "pca.eigenvec"), sep=".") )
#
# rownames(tb)    = tb[,1]
# tb    = tb[, -c(1,2)]
#     colnames(tb)    = paste0("PCA", 1:NCOL(tb) )
#
# return(tb)

# scaffoldId, Loci, RefBase, AltBase, P1, P2, P3,....
def load_simple_snp(objref, verbose=True):
    if not objref.file_simple_snp == "":
        tb_gen = pd.read_table(objref.file_simple_snp)
    else:
        tb_gen = objref.rawData

    if verbose:
        print("gen.table[", tb_gen.shape, "]\n")

    snp_info = tb_gen.iloc[:, 0:5]
    snp_mat = tb_gen.iloc[:, 5:tb_gen.shape[1]]
    snp_info.index = snp_info.iloc[:, 0]
    objref.snpData.snp_info = snp_info

    objref.snpData.snp_mat = snp_mat.iloc[:, objref.ids_used==snp_mat.columns.values]
    objref.snpData.snp_mat = snp_mat.iloc[:, objref.ids_used==snp_mat.columns.values]
