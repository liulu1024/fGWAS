class FgGenObj(object):
    n_snp = None
    n_ind_total = None
    n_ind_used = None
    reader = None
    files = None
    options = None
    #(file_plink_bed,file_plink_bim, file_plink_fam,command )
    params = None


class FgPheObj(object):
    obj_curve = None
    obj_covar = None
    ids = None
    pheY = None
    pheX = None
    pheT = None
    n_ids = None


# todo add get set method

# for function proc_dat_simu
class SimuObj(object):
    def __init__(self, obj_gen, obj_phe, error):
        self.obj_gen = obj_gen
        self.obj_phe = obj_phe
        self.error = error


class SnpData(object):
    def __init__(self, snp_mat, snp_info):
        self.snp_mat = snp_mat
        self.snp_infp = snp_info
