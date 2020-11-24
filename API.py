from Simulation.simulate import fg_simulate


def siumulate(curve_type, covariance_type, n_obs, n_snp, time_points, par0=None, par1=None, par2=None,
              par_covar=None, par_X=None, sig_pos=None, phe_missing=0.03, snp_missing=0.03, plink_format=False,
              file_prefix=None):
    # TODO  param check

    obj = fg_simulate(curve_type, covariance_type, n_obs, n_snp, time_points, par0, par1, par2, par_covar, par_X,
                      phe_missing, snp_missing, sig_pos, plink_format, file_prefix)
    return obj


def check_param():
    pass
#     ## check 'curveType'
#     if math.isnan(curveType):
#         sys.exit("No curve type is specified in the paramater 'curveType'.")
#     else:
#         if isinstance(curveType, str) and len(curveType) > 1:
#             sys.exit("curveType should be a string indicating curve type.")
#
#     ## check 'covarianceType'
#     if math.isnan(covarianceType):
#         print("No covariance type is specified in the paramater 'covarianceType'.")
#     else:
#         if isinstance(covarianceType, str) and len(covarianceType) > 1:
#             sys.exit("'covarianceType' should be a string indicating covariance type.");
#
#     ## check 'n_obs'
#     if math.isnan(n_obs):
#         sys.exit("Individual count is not specified in the paramater 'n_obs'.")
#     else:
#         if not isinstance(n_obs, int):
#             sys.exit("Not integer in parameter 'n_obs'.")
#
#     ## check 'n_snp'
#     if math.isnan(n_snp):
#         sys.exit("SNP count is not specified in the paramater 'n_snp'.")
#     else:
#         if not isinstance(n_snp, int):
#             sys.exit("Not integer in parameter 'n_snp'.")
#
#     ## check 'time_points'
#     if math.isnan(time_points):
#         sys.exit("Time points are not specified in the paramater 'time_points'.")
#     else:
#         if not isinstance(time_points, int) or not isinstance(time_points, list):
#             sys.exit("Not integer or list in parameter 'time_points'.")
#         #
#         if len(time_points) == 1:
#             time_points = np.arrange(0, time_points + 1, 1)
#
#     ## check 'par_X'
#     if math.isnan(par_X):
#         # TODO deal par_X type
#         if np.dtype(par_X):
#             sys.exit("Covariate paramater should be numeric values.")
#
#     ## check 'file.prefix'
#     if math.isnan(file.prefix):
#         if not isinstance(file_prefix, str) or len(file_prefix) > 1:
#             sys.exit("'file.prefix' should be a string.")
#
#     if plink_format and math.isnan(file.prefix):
#         sys.exit("'file.prefix' is NULL")
#
#     ## check 'curveType'
#     if curveType.__class__.__name__ == "BaseCurve":
#         fg_curve = curveType
#     else:
#         fg_curve = base.getCurve(curveType)
#     # FG_CURVE GET_SIMU_PARAM
#     simu_param = base.get_simu_param(fg_curve, time_points)
#     simu_len = np.shape(simu_param)[1]
#
#     ## check 'par0'
#     if math.isnan(par0):
#         par0 = simu_param[0, :]
#
#         # TODO also par0 dtype
#         # else:
#         # 	if (!( all( is.numeric(par0)) & & length(par0) == simu_len) )
#         #     sys.exit(paste("Curve paramater should be", simu_len, "numeric values."))
#
#         ## check 'par1'
#
#     if math.isnan(par1):
#         par0 = simu_param[1, :]
#     # else
#     # {
#     # if (!( all( is.numeric(par1)) & & length(par1) == simu_len) )
#     # sys.exit(paste("Curve paramater should be", simu_len, "numeric values."));
#     # }
#
#     ## check 'par2'
#     if math.isnan(par2):
#         par2 = simu_param[2, :]
#     # else
#     # {
#     # if (!( all( is.numeric(par2)) & & length(par2) == simu_len) )
#     # sys.exit(paste("Curve paramater should be", simu_len, "numeric values."));
#     # }
#
#     ## check 'covarianceType'
#
#
# if (
#
#
# class(covarianceType)== "fgwas.covar")
# fg_covar < - covarianceType
# else
# fg_covar < - fg.getCovariance( covarianceType );
#
# simu_covar < - get_simu_param(fg_covar, time_points);
# simu_len < - NROW(simu_covar);
#
# ## check 'par.covar'
# if (missing(par.covar) | | is.null(par.covar) )
# par.covar < - simu_covar
# else
# {
#
#
# if (!( all( is.numeric(par.covar)) & & length(par.covar) == simu_len) )
# sys.exit(paste("Covariance paramater should be", simu_len, "numeric numbers."));
# }
#
# ## check 'sig.pos'
# if (missing(sig.pos) | | is.null(sig.pos))
# {
# sig.pos < - round(runif(1, n_snp * 0.25, n_snp * 0.75));
# cat(" * A significant SNP is randomly specified to location(", sig.pos, ")\n");
# }
# else
# {
# if (!( is.numeric( sig.pos) & & length(sig.pos)) )
# sys.exit("Not integer in parameter 'n_snp'.");
# }
# # optional items: verbose=F
# #                forece_split=T
#
# fg_load_plink < -function(file_plink_bed, file_plink_bim, file_plink_fam,
#                           plink_command = None, chr = None, options = list(verbose=F))
# {
# ## dont load plink data, just plink data object
# obj < - fg_load_plink(file_plink_bed, file_plink_bim, file_plink_fam, plink_command, chr, options);
#
# return (obj);
# }
#
# # optional items: verbose=F
# fg_load_simple < - function(file_simple_snp, options=list(verbose=F))
# {
# obj < - fg_load_simple(file_simple_snp, options)
# return (obj);
# }
#
# # optional items: verbose=F
# #                max_optim_failure= 100
# #                min_optim_success= 20,
# #                R2_loop = 5,
# fg_load_phenotype < -function(file_phe_long,
#                               file_phe_cov = None, file_phe_time = None, curve_type = None, covariance_type = None, file_plot_pdf = None, intercept = TRUE, options = list(
#     verbose=F))
# {
# if (! is_None(curve_type) & & tolower(curve_type) != "auto" & & is_None( fg_getCurve( curve_type ) ))
# stop("No curve object specified by the parameter 'curve_type'_ ")
#
# if (! is_None(covariance_type) & & tolower(
#         covariance_type) != "auto" & & is_None( fg_getCovariance( covariance_type ) ))
# stop("No covariance object specified by the parameter 'covariance_type'_ ")
#
# obj < - fg_load_phenotype(file_phe_long, file_phe_cov, file_phe_time, curve_type, covariance_type, file_plot_pdf,
#                           intercept, options);
# return (obj);
# }
#
# # optional items: verbose=F
# #                min_optim_failure= 100
# #                min_optim_success= 20,
# #                R2_loop = 5,
# fg_data_estimate < -function(obj_phe,
#                              curve_type = "auto", covariance_type = "auto", file_plot_pdf = None, options = list(
#     verbose=F) )
# {
# if (! is_None(curve_type) & & tolower(curve_type) != "auto" & & is_None( fg_getCurve( curve_type ) ))
# stop("No curve object specified by the parameter 'curve_type'_ ")
#
# if (! is_None(covariance_type) & & tolower(
#         covariance_type) != "auto" & & is_None( fg_getCovariance( covariance_type ) ))
# stop("No covariance object specified by the parameter 'covariance_type'_ ")
#
# obj < - fg_dat_est(obj_phe, curve_type, covariance_type, file_plot_pdf, options);
# return (obj);
# }
#
# # optional items: verbose=FALSE
# #                ncores=1,
# #                use_snowfall=TRUE
# #                max_optim_failure=20 for fgwas and optim-fgwas
# #                min_optim_success=2 for fgwas and optim-fgwas
# #                use_gradient=T for fgwas and optim-fgwas
# #                piecewise=1000  for fgwas and optim-fgwas
# #                degree=3 for mixed
# fg_snpscan < -function(fgwas_gen_obj, fgwas_phe_obj, method="optim-fgwas",
#                        curve_type = None, covariance_type = None, snp_sub = None, options = list(verbose=F))
# {
# if (! toupper(method) % in% toupper(c("gls", "mixed", "fast", "fgwas", "optim-fgwas")))
# stop("the parameter 'method' has 5 optional values: 'gls', 'mixed', 'fast', 'fgwas', 'optim-fgwas'_ ")
#
# if (! is_None(curve_type) & & tolower(curve_type) != "auto" & & is_None( fg_getCurve( curve_type ) ))
# stop("No curve object specified by the parameter 'curve_type'_ ")
#
# if (! is_None(covariance_type) & & tolower(
#         covariance_type) != "auto" & & is_None( fg_getCovariance( covariance_type ) ))
# stop("No covariance object specified by the parameter 'covariance_type'_ ")
#
# ret < - fg_snpscan(fgwas_gen_obj, fgwas_phe_obj, method=method,
#                    curve_type = curve_type, covariance_type = covariance_type, permutation = None, snp_sub = snp_sub, options = options );
# return (ret);
# }
#
#
# fg_select_sigsnp < - function(fgwas_scan_obj, sig_level = 0_05, pv_adjust = "bonferroni", options = list() )
# {
# ret < - fg_select_sigsnp(fgwas_scan_obj, sig_level, pv_adjust, options=list());
# return (ret);
# }
#
# summary_fgwas_gen_obj < - function(object, ___)
# {
# stopifnot(
#
#
# class(object)== "fgwas_gen_obj");
# summary_fgwas_gen_obj(object)
# }
#
#
# summary_fgwas_phe_obj < - function(object, ___)
# {
# stopifnot( class (object) == "fgwas_phe_obj");
# summary_fgwas_phe_obj(object)
# }
#
# summary_fgwas_scan_obj < - function(object, ___)
# {
# stopifnot( class (object) == "fgwas_scan_obj");
# summary_fgwas_scan_obj(object)
# }
#
# print_fgwas_gen_obj < - function(x, ___, useS4 = FALSE)
# {
# stopifnot( class (x) == "fgwas_gen_obj");
# print_fgwas_gen_obj(x  )
# }
#
# print_fgwas_phe_obj < - function(x, ___, useS4 = FALSE)
# {
# stopifnot( class (x) == "fgwas_phe_obj");
# print_fgwas_phe_obj(x )
# }
#
# print_fgwas_scan_obj < - function(x, ___, useS4 = FALSE )
# {
# stopifnot( class (x) == "fgwas_scan_obj");
# print_fgwas_scan_obj(x )
# }
#
# plot_fgwas_scan_obj < - function(x, y, ___, file_pdf=None, sig_level=0_05)
# {
# stopifnot( class (x) == "fgwas_scan_obj");
# plot_fgwas_scan_obj(x, y, file_pdf, sig_level, ___);
# }
#
# plot_fgwas_phe_obj < - function(x, y, ___, curve_fitting=T, file_pdf=None)
# {
# stopifnot( class (x) == "fgwas_phe_obj");
# plot_fgwas_phe_obj(x, file_pdf, curve_fitting, ___)
# }
#
# plot_fgwas_curve < - function( object, snp_sub, file_pdf=None, draw_rawdata=TRUE, draw_meanvector=TRUE, ___)
# {
# stopifnot( class (object) == "fgwas_scan_obj");
#
# if ( is_None(object$ret_fast) & & is_None(object$ret_fgwas))
# stop("No curve information in the result, only for FAST and fGWAS model");
#
# plot_fgwas_curve( object, snp_sub, file_pdf, draw_rawdata, draw_meanvector, ___);
# }
#
# # Inner function, not public
# profile_fgwas_curve < - function( object, snp_sub )
# {
# stopifnot( class (object) == "fgwas_scan_obj");
#
# if ( is_None(object$ret_fast) & & is_None(object$ret_fgwas))
# stop("No curve information in the result, only for FAST and fGWAS model");
#
# r < - profile_fgwas_curve( object, snp_sub );
#
#
# return (r);
# }
#
# fg_get_pca < -function(object, plink_path)
# {
# if (
#
#
#     class(object$reader) == "fg_dm_simple")
#
#
# return (fg_simple_getPCA(object, plink_path));
#
# if (
#
#
#     class(object$reader) == "fg_dm_plink")
#
#
# return (fg_plink_getPCA(object, plink_path));
# }
#
# fg_get_snp < - function(object, snp_names, options=list())
# {
# snp_idx < - object$reader$get_snpindex(snp_names);
# if (any( is_na(snp_idx)))
# stop("No snp in genome data", snp_names[which( is_na(snp_idx))], "\n");
#
# snp < - object$reader$get_snpmat(snp_idx);
# colnames(snp$snpmat) < - snp_names;
# return (snp);
# }
#
# fg_adjust_inflation < - function(object)
# {
# stopifnot(
#
#
# class(object)== "fgwas_scan_obj");
#
# if ( is_None(object$ret_fast) & & is_None(object$ret_fgwas) & & is_None(object$ret_gls))
# stop("No curve information in the result, only for FAST and fGWAS model");
#
# object < - adjust_fgwas_genomic_inflation( object );
#
#
# return (object);
# }
#
# fg_qqplot < -function(object, png_file, title="", width=480)
# {
# stopifnot(
#
#
# class(object)== "fgwas_scan_obj");
#
# fg_qqplot( object, png_file, title, width );
#
# }
#
#
