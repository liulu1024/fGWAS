# Phenotype File( CSV Format, Header=T  Comma seprated )
#
# Format 1:
# shareid, PHE1, PHE2, ...
#
# Covariate File( CSV Format, Header=T  Comma seprated  )
# Format 2:
#
# shareid, COV1,....
#
# Format 3:
# shareid, TIME1, TIME2, ...
import sys
from os.path import exists
import Simulation.FgObj as fgObj
import ResultVo as result
from Estimate import fg_dat_est
from FgCommon import paste
import pandas as pd
import numpy as np


def check_pheno_file(file_phe_long, file_phe_cov, file_phe_time, verbose=False):
    err_info = "Error!"
    if verbose:
        print("Checking phenotype file......\n")
        print("* phe_long =", file_phe_long, "\n")

    if not exists(file_phe_long):
        return result.NormalResultVo(True, paste(["Longitudinal file is not found, file=", file_phe_long], sep=""))

    try:
        phe_long = pd.read_csv(file_phe_long)
    except:
        return result.NormalResultVo(True, err_info)
    if verbose:
        print("* Individuals =", phe_long.index.size, "\n")
        print("* Times =", phe_long.columns.size, "\n")
        print("* Mean =", phe_long.mean(), "\n")
        print("* SD =", phe_long.std(), "\n")

        print("  First 5 Items:\n")
        print(phe_long.head())

    if verbose: print("* phe_time =", file_phe_time, "\n")

    phe_time = None
    if file_phe_time is not None and file_phe_time is not None:
        if not exists(file_phe_time):
            return result.NormalResultVo(True, paste(["Time file is not found, file=", file_phe_time], sep=""))
        try:

            phe_time = pd.read_csv(file_phe_time)

        except:
            print("! Can not open file(", file_phe_time, ")\n")
            return result.NormalResultVo(error=True, err_info=err_info)

        if verbose:
            print("* Individuals =", phe_time.index.size, "\n")
            print("* Times =", phe_time.columns.size, "\n")
            print("* Mean =", phe_time.mean(), "\n")
            print("* SD =", phe_time.mean(), "\n")
            print("  First 5 Items:\n")
            print(phe_time.head())

        idx_inter = np.intersect1d(phe_long.index.values, phe_time.index.values())
        if not (len(idx_inter) == phe_long.index.size and len(idx_inter) == phe_time.index.size):
            print("! phe_long don't have consistent IDs with phe_time.\n")
            return result.NormalResultVo(error=True, err_info=err_info)

    if verbose: print("Checking covariate file......\n")
    if verbose: print("* COV.FILE =", file_phe_cov, "\n")

    if file_phe_cov is not None:
        if not exists(file_phe_cov):
            return result.NormalResultVo(True, paste(["Covariate file is not found, file=", file_phe_cov], sep=""))
        try:
            phe_cov = pd.read_csv(file_phe_cov)
        except:
            print("! Can not open file(", file_phe_cov, ")\n")
            return result.NormalResultVo(True, err_info)
        if verbose:
            print("* Individuals =", phe_cov.index.size, "\n")
            print("* Covariates =", phe_cov.columns.size, "\n")
            print("* Mean =", phe_cov.mean(), "\n")
            print("* SD =", phe_cov.std(), "\n")
            print("  First 5 Items:\n")
            print(phe_cov.head())

        # y_ncov = phe_cov.columns.size
    # all.na < - which( is.na(rowSums(phe.long, na.rm = T)    ) )
    # if (length(all.na) > 0)
    #     {
    #         cat("!", length(all.na), "IDs dont have non-missing data.\n");
    #     cat("! First 5 IDs:\n");
    #     show(head(phe.long[all.na,], n=5));
    #     }
    return result.NormalResultVo(error=False, err_info=None)


# currently, no option items are used in this funtion
def fg_load_phenotype(file_phe_long, file_phe_cov=None, file_phe_time=None, curve_type="auto", covariance_type="auto",
                      file_plot_pdf=None, intercept=True, **kwargs):
    default_options = {"max_optim_failure": 100, "min_optim_success": 20, "R2_loop": 5, "verbose": False}
    if not kwargs:
        kwargs = default_options
    else:
        options0 = default_options
        options0.update(kwargs)
        kwargs = options0

    # if (missing(file_phe_long))
    #    os.exit("! file_phe_long must be assigned with the valid file name.")

    if curve_type is None: curve_type = "auto"
    if covariance_type is None: covariance_type = "auto"
    if file_phe_cov is None:  file_phe_cov = None
    if file_phe_time is None: file_phe_time = None
    if intercept is None: intercept = True

    if kwargs["verbose"]:
        print("[ Loading Phenotype Data]\n")
        print("Checking the parameters ......\n")
        print("* Phenotypic Data File = ", file_phe_long, "\n")
        print("* Covariate Data File = ", file_phe_cov, "\n")
        print("* Time Data File = ", file_phe_time, "\n")
        print("* Curve Type = ", curve_type, "\n")
        print("* Covariance Type = ", covariance_type, "\n")
        print("* Intercept = ", intercept, "\n")

        # print( "Checking the optional items......\n")
        # show_options( options)
        # pheObj
    ret = fgObj.FgPheObj()
    ret.intercept = intercept
    ret.options = kwargs
    ret.params = {"file_phe_long": file_phe_long,
                  "file_phe_cov": file_phe_cov,
                  "file_phe_time": file_phe_time,
                  "curve_typ": curve_type,
                  "covariance_type": covariance_type}

    r_chk = check_pheno_file(file_phe_long, file_phe_cov, file_phe_time, kwargs["verbose"])
    if r_chk.error:
        sys.exit(r_chk.err_info)
        ##### pheY file_phe_long 文件
    ret.pheY = pd.read_csv(file_phe_long)
    if file_phe_time is not None:
        ret.pheT = pd.read_csv(file_phe_time)
    else:
        row = ret.pheY.index.size
        pheT_data = np.tile(np.range(0, row), row).reshape(row, row, order='f')
        ret.pheT = pd.DataFrame(pheT_data, index=ret.pheY.index.values)

    if file_phe_cov is not None:
        ret.pheX = pd.read_csv(file_phe_cov)
    else:
        ret.pheX = None

        ## check the id consistent
    ret.ids = ret.pheY.index.values

    # if (! is.None(ret.pheX) & & is.data.frame(ret.pheX)) ret.pheX = as.matrix(ret.pheX)
    # if (! is.None(ret.pheY) & & is.data.frame(ret.pheY)) ret.pheY = as.matrix(ret.pheY)
    # if (! is.None(ret.pheT) & & is.data.frame(ret.pheT)) ret.pheT = as.matrix(ret.pheT)

    r_est = fg_dat_est(ret, curve_type, covariance_type, file_plot_pdf, **kwargs)
    r_est = None
    if r_est.error:
        sys.exit(r_est.err_info)
    else:
        ret.obj_curve = r_est.obj_curve
        ret.obj_covar = r_est.obj_covar
        ret.est_curve = r_est.est_curve
        ret.est_covar = r_est.est_covar
        ret.summary_curve = r_est.summary_curve
        ret.summary_covariance = r_est.summary_covar

    return ret
