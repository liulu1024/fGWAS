# S4 class for storing TFBS.
# sequence preferences of a set of TFs.
import numpy as np
import os
import Simple as simple

class FgGenObj(object):
    def __init__(self, n_snp=None, n_ind_total=None, n_ind_used=None, reader=None, files=None, options=None, params=None):
        self.n_snp = n_snp
        self.n_ind_total = n_ind_total
        self.n_ind_used = n_ind_used
        self.reader = reader
        self.files = files
        self.options = options
        # (file_plink_bed,file_plink_bim, file_plink_fam,command )
        self.params = params


class FgPheObj(object):
    def __init__(self, obj_curve=None, obj_covar=None, ids=None, pheY=None,
                 pheX=None, pheT=None, n_ids=None, intercept=None,params=None,
                 est_curve=None,est_covar=None,summary_curve=None,summary_covar=None,h0=None):
        self.obj_curve = obj_curve
        self.obj_covar = obj_covar
        self.ids = ids
        self.pheY = pheY
        self.pheX = pheX
        self.pheT = pheT
        self.n_ids = n_ids
        self.intercept = intercept
        self.params = params
        self.est_curve = est_curve
        self.est_covar = est_covar
        self.summary_curve = summary_curve
        self.summary_covariance = summary_covar
        self.h0=h0
       # self.pheZ=pheZ if pheZ is None else pd.DataFrame()



# for function proc_dat_simu
class SimuObj(object):
    def __init__(self, obj_gen=None, obj_phe=None, error=None):
        self.obj_gen = obj_gen
        self.obj_phe = obj_phe
        self.error = error


class SnpData(object):
    def __init__(self, snp_mat=None, snp_info=None):
        self.snp_mat = snp_mat
        self.snp_info = snp_info


class FgGenoData(object):
    """
    type: plink,cvf,simple
    desc:plink,cvf,simple
    snp: number of snp in data file
    ind_total:Number of individuals in data file
    ind_used:Number of individuals used
    ids_used:individual used.
    """

    # type, description, snp_num, ind_total_num, ind_used_num, ids_used
    def __init__(self, type, description, snp_num, ind_total_num, ind_used_num, ids_used):
        self.type = type
        self.description = description
        self.snp_num = snp_num
        self.ind_total_num = ind_total_num
        self.ind_used_num = ind_used_num
        self.ids_used = ids_used

    def show(self):
        print("Reference matrix object of class:")
        print(type(self))
        print("  Data type: " + self.type)
        print("  Description: " + self.description)
        print("  SNP Count:0}".format(self.snp_num))
        print("  Individual Count:0}".format(self.ind_total_num))
        print("  Individual Used:0}".format(self.ind_used_num))

    def shrink(self):
        pass

    def get_snpinfo(self, snp_idx):
        pass

    def get_snpmat(snp_idx, impute=False, allel=False):
        pass

    def select_individuals(self, ids_used):
        pass

    def get_individuals(self):
        pass

    def get_used_individuals(self):
        pass

    def get_snpindex(self, snp_names):
        pass


'''
Simple obj
'''


class Snp(object):
    def __init__(self, snp=None, NMISS=None, MAF=None):
        self.snp = snp
        self.NMISS = NMISS
        self.MAF = MAF


class SnpMat(object):
    def __init__(self, snpmat=None, NMISS=None, MAF=None):
        self.snpmat = snpmat
        self.NMISS = NMISS
        self.MAF = MAF


class FgDmSimple(FgGenoData):

    def __init__(self, type, description, n_snp, n_ind_total, n_ind_used, ids_used, file_simple_snp, rawData,
                     snpData):
            self.type = type
            self.description = description
            self.n_snp = n_snp
            self.n_ind_total = n_ind_total
            self.n_ind_used = n_ind_used
            self.ids_used = ids_used
            self.file_simple_snp = file_simple_snp
            self.rawData = rawData
            self.snpData = snpData

    def shrink(self):
        super().shrink(self)

    def get_snpinfo(self, snp_idx):
        return simple.simple_get_snpinfo(self, snp_idx)

    def get_snpmat(self, snp_idx, impute=False, allel=False):
        return simple.simple_get_snpmat(self, snp_idx, impute=impute, allel=allel)

    def select_individuals(self, ids_used):
        return simple.simple_select_individuals(self, ids_used)

    def get_used_individuals(self):
        return simple.simple_get_used_individuals(self)

    def get_individuals(self):
        return simple.simple_get_individuals(self)

    def get_snpindex(self, snp_names):
        return simple.simple_get_snpindex(self, snp_names)

# def select_individuals(object, ids_used):
#     idx_match = np.where(ids_used.astype(str), object.ids)
#     if not len(idx_match) == len(object.ids):
#         os.exit("Some IDs are not matached in the SNP data file"
# 
#     object.n_ind = len(ids_used)
#     object.ids = ids_used.astype(str)
#     object.pheY = object.pheY.iloc[idx_match]
#     object.pheX = object.pheX.iloc[idx_match]
#     object.pheZ = object.pheZ.iloc[idx_match]
# 
#     return object
