# import Simulation.simulate
#
# # base.get_all_curve()
# # base.get_all_curve()
# # base.get_all_curve()
# Simulation.simulate.fg_simulate("Logistic", "ARMA1", 200, 100, range(0, 8), phe_missing=0.05, snp_missing=0.05,
#                                 sig_pos=51, plink_format=False, file_prefix='test', par_X=(2, 3.2))
from pyplink import PyPlink

with PyPlink("plink_file_prefix") as bed:
    # Getting the BIM and FAM
    bim = bed.get_bim()
    fam = bed.get_fam()

    # Iterating over all loci
    for loci_name, genotypes in bed:
        pass

    # Getting the genotypes of a single marker (numpy.ndarray)
    genotypes = bed.get_geno_marker("rs12345")
