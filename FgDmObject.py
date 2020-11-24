# S4 class for storing TFBS.
# sequence preferences of a set of TFs.
import numpy as np
import os


class FgGenoData(object):
    """
    type: plink,cvf,simple
    desc:plink,cvf,simple
    snp: number of snp in data file
    ind_total:Number of individuals in data file
    ind_used:Number of individuals used
    ids_used:individual used.
    """

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
        print("  SNP Count: {0}".format(self.snp_num))
        print("  Individual Count: {0}".format(self.ind_total_num))
        print("  Individual Used: {0}".format(self.ind_used_num))

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


class FgDmPlink(FgGenoData):
    """
    file_bed   character
    file_bim   character
    file._fam   character
    command      character
    chromosome    numeric
    snp_block_size   numeric
    snp_data         list
    """

    def __init__(self, file_bed, file_bim, file_fam, command, chromosome, snp_block_size, snp_data):
        self.file_bed = file_bed
        self.file_bim = file_bim
        self.file_fam = file_fam
        self.command = command
        self.chromosome = chromosome
        self.snp_block_size = snp_block_size
        self.snp_data = snp_data

    def show(self):
        if not self.chromosome == -1:
            print("  Chromosome: {0}".format(self.chromosome))
        else:
            print("  Chromosome: all")
        super(FgDmPlink, self).show()
        print("  Plink Command: " + self.command)

    def shrink(self):
        pass


# shrink = function()
# {
#
#
# return (plink.shrink(.self));
# },
#
# get_snpinfo = function(snp.idx)
# {
# return (plink.get.snpinfo(.self, snp.idx));
# },
#
# get_snpmat = function(snp.idx, impute=F, allel=F)
# {
# return (plink.get.snpmat(.self, snp.idx, impute=impute, allel=allel));
# },
#
# select_individuals = function(ids.used)
# {
# return (plink.select.individuals(.self, ids.used));
# },
#
# get_used_individuals = function()
# {
# return (plink.get.used.individuals(.self));
# },
#
# get_individuals = function()
# {
# return (plink.get.individuals(.self));
# },
#
# get_snpindex = function(snp.names)
# {
# return (plink.get.snpindex(.self, snp.names));
# } )


class FgDmSimple(FgGenoData):
    def __init__(self, file_simple_snp, rawData, snpData, *args):
        super(FgDmSimple, self).__init__(args)
        self.file_simple_snp = file_simple_snp
        self.rawData = rawData
        self.snpData = snpData


# shrink = function()
# {
# return (simple.shrink(.self));
# },
#
# get_snpinfo = function(snp.idx)
# {
# return (simple.get.snpinfo(.self, snp.idx));
# },
#
# get_snpmat = function(snp.idx, impute=F, allel=F)
# {
# return (simple.get.snpmat(.self, snp.idx, impute=impute, allel=allel));
# },
#
# select_individuals = function(ids.used)
# {
# return (simple.select.individuals(.self, ids.used));
# },
#
# get_used_individuals = function()
# {
# return (simple.get.used.individuals(.self));
# },
#
# get_individuals = function()
# {
# return (simple.get.individuals(.self));
# },
#
# get_snpindex = function(snp.names)
# {
# return (simple.get.snpindex(.self, snp.names));
# } )
#
# )


def select_individuals(object, ids_used):
    idx_match = np.where(ids_used.astype(str), object.ids)
    if not len(idx_match) == len(object.ids):
        os.exit("Some IDs are not matached in the SNP data file");

    object.n_ind = len(ids_used)
    object.ids = ids_used.astype(str)
    object.pheY = object.pheY.iloc[idx_match]
    object.pheX = object.pheX.iloc[idx_match]
    object.pheZ = object.pheZ.iloc[idx_match]

    return object
