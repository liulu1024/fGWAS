# import os
# import  Simulation.FgObj as FgObj
# from pandas_plink  import  read_plink
# import  FgDmObject as FgDmObj
# from pyplink import  PyPlink
# import  numpy as np
# ## PLINK format
# #
# ## BIM file:
# # Chromosome
# # Marker ID
# # Genetic distance
# # Physical position
# # Allele 1
# # Allele 2
# #
# #
# #
# # Example of a BIM file of the binary PLINK format:
# #
# # 21 rs11511647 0 26765 A T
# # X      rs3883674 0 32380 C G
# # X     rs12218882 0 48172 T T
# # 9     rs10904045 0 48426 A T
# # 9     rs10751931 0 49949 C T
# # 8     rs11252127 0 52087 A C
# # 10 rs12775203 0 52277 A A
# # 8     rs12255619 0 52481 G T
# #
# ##
# ## FAM files
# #
# # Family ID
# # Sample ID
# # Paternal ID
# # Maternal ID
# # Sex (1=male 2=female other=unknown)
# # Affection (0=unknown 1=unaffected 2=affected)
# #
# #
# # FAM1 NA06985 0 0 1 1
# # FAM1 NA06991 0 0 1 1
# # 0  NA06993 0 0 1 1
# # 0  NA06994 0 0 1 1
# # 0  NA07000 0 0 2 1
# # 0  NA07019 0 0 1 1
# # 0  NA07022 0 0 2 1
# # 0  NA07029 0 0 1 1
# # FAM2 NA07056 0 0 0 2
# # FAM2 NA07345 0 0 1 1
#
#
# def fg_load_plink(file_plink_bed=None, file_plink_bim=None, file_plink_fam=None,command=None, chr=None, **kargs):
#         # 传入参数修改
#      # default_options ='force.split':True, 'verbose':False
#      # default_options[names(options)] =options
#      # options =default_options
#
#      if (file_plink_bed is None) or  (file_plink_bim is None) or (file_plink_fam is None):
#         os.exit("! file_plink_bed, file_plink_bim, file_plink_fam must be assigned with the valid values.")
#
#      force_split = kargs['force_split']
#      verbose = kargs['verbose']
#      if not isinstance(force_split,bool):
#         os.exit("! The parameter of force.split should be a logical value(True or False).")
#
#      if verbose:
#        print( "[ Loading PLINK ]")
#        print( "Checking the parameters ......")
#
#        print("* PLINK BED File = "+file_plink_bed)
#        print("* PLINK BIM File = "+file_plink_bim)
#        print("* PLINK FAM File = "+file_plink_fam)
#        print("* Chromosome     = "+ chr)
#        print("* PLINK Command = "+command)
#        print("* Force Split by PLINK Command = "+force_split)
#
#
#
#      chk = check_plink_file(file_plink_bed, file_plink_bim, file_plink_fam, verbose)
#      if chk.error:
#         os.exit("! Failed to check PLINK data files.")
#
#      obj_gen = FgObj.FgGenObj()
#      obj_gen.reader =  create_plink_obj ( file_plink_bed, file_plink_bim, file_plink_fam,command, chromosome=chr, verbose=verbose)
#
#      obj_gen.n_snp  =obj_gen.reader.n_snp
#      obj_gen.n_ind_total =obj_gen.reader.n_ind_total
#      obj_gen.n_ind_used  =len(obj_gen.reader.ids_used)
#
#      obj_gen.options = kargs
#      obj_gen.params  = (file_plink_bed,file_plink_bim, file_plink_fam,command )
#
#
#      return obj_gen
#
# def check_plink_file(file_plink_bed, file_plink_bim, file_plink_fam, verbose=False):
#
#      if verbose:
#
#        print("Checking PLINK file......")
#        print("* BED file ="+file_plink_bed)
#        print("* BIM file ="+file_plink_bim)
#        print("* FAM file ="+file_plink_fam)
#
#
#      bigdata  =False
#      error =False
#     #TODO read plink
#      snp_mat =None
#      try :
#         read_plink( file_plink_bed,  file_plink_bim, file_plink_fam)
#      except Exception:
#
#           if os.path.isfile(file_plink_bed) and os.path.isfile(file_plink_bim) and os.path.isfile(file_plink_fam ):
#            bigdata =True
#           else:
#            error =True
#
#
#      if not error:
#
#       # tb_fam =read.table(file_plink_fam, header=F)
#       # tb_bim =read.table(file_plink_bim, header=F)
#       # n_idv =NROW(tb.fam)
#       # n_snp =NROW(tb.bim)
#       tb_fam = None
#       tb_bim = None
#       n_idv = None
#       n_snp = None
#
#       if verbose:
#
#         print("* Individuals ="+n_idv+"SNPs="+n_snp)
#         print("* PLINK testing successfully.")
#
#       return {'error':False, 'bigdata':bigdata, 'family':tb_fam[:,2]}
#
#      else:
#       return {'error':True}
#
#
#
#
# def create_plink_obj(file_plink_bed=None, file_plink_bim=None, file_plink_fam=None, command=None, chromosome=None, verbose=False):
#     # TODO check if None
#      # if ( missing( file_plink_bed) || missing(file_plink_bim) || missing(file_plink_fam) )
#      #  stop("! file_plink_bed, file_plink_bim and file_plink_fam must be assigned with the valid values.")
#
#      # if(is.None(plink.command)) plink.command ="plink"
#
#      objref =FgDmObj.FgDmPlink(
#       file_bed    = file_plink_bed,  ##
#       file_bim    = file_plink_bim,  ##
#       file_fam    = file_plink_fam,  ##
#       command     = '' if command is None else command,   ##
#       chromosome        = -1 if chromosome is None else chromosome,
#       snp_block_size     = 1000,
#       snp_data           = [])
#
#      # t =try(system( paste( plink.command, "--noweb", sep=" "), ignore.stdout=True, ignore.stderr=True ))
#      # if(class(t)=="try-error")
#      #
#      #   print("! No PLINK command can be found in your environment( plink.command=",plink.command, ")")
#      #  return(None)
#
#
#      objref =plink_checkFiles( objref, verbose )
#
#      return objref
#
#
# def plink_checkFiles(refobj, verbose=False):
#
#      if verbose:
#
#        print( "Checking PLINK files ......")
#        print("* PLINK BED File = "+refobj.file_plink_bed)
#        print("* PLINK BIM File = "+refobj.file_plink_bim)
#        print("* PLINK FAM File = "+refobj.file_plink_fam)
#        print("* PLINK Command = "+refobj.plink.command)
#
#
#      if os.path.isfile(refobj.file_plink_bed) and os.path.isfile(refobj.file_plink_bim) and os.path.isfile(refobj.file_plink_fam)
#
#          with PyPlink("plink_file_prefix") as bed:
#              # Getting the BIM and FAM
#              bim = bed.get_bim()
#              fam = bed.get_fam()
#
#
#          refobj.snp_data.plink_fam =fam
#          fam.columns=["FID", "IID", "PID", "MID", "SEX", "PHE"]
#
#          refobj.snp_data.plink_bim =bim
#          refobj.snpdata.plink_bim=["CHR", "SNP", "CM", "POS", "A1", "A2"]
#
#          if not refobj.chromosome == -1:
#
#            snp_idx =np.isin(refobj.snp_data.plink_bim.loc[:,0], refobj.chromosome)
#            refobj.snp_data.plink_bim =refobj.snp_data.plink_bim.loc[snp_idx,:]
#
#
#          refobj.n.snp = refobj.snp_data.plink_bim.shape[0]
#          refobj.n_ind_total  = refobj.snp_data.plink_fam.shape[0]
#          refobj.n_ind.used   =refobj.n_ind_total
#          refobj.ids.used     = (refobj.snp_data.plink_fam.loc[:,1]).values.astype(str)
#
#          return refobj
#
#      else:
#       os.exit("PLInK files can not be found!")
#
#
#
#
#
#
#
#
# def plink_getSnpMat(objref, snp_idx=None, impute=False, allel=False ):
#        if  snp_idx is None:
#             snp_idx = np.arange(0,objref.n_snp)
#
#         # in decreasing order
#        snp_idx.new =snp_idx.sort()
#
#        def check_local_snpmat(snp_k):
#           if objref.snp_data.local_idx is None:
#              return False
#
#           return np.isin(snp_k,objref.snp_data.local_idx)
#
#
#        def load_local_snpmat(snp_k):
#
#         select_subjects =objref.snpdata.plink.fam[match(objref.ids.used, objref.snpdata.plink.fam[, 2] ), c(1, 2)]
#         loading.idx =seq(snp.k,
#                             ifelse(snp.k + objref.snp.blocksize > objref.n.snp, objref.n.snp, snp.k + objref.snp.blocksize), 1);
#         plink.obj =plink.cmd.select(
#             objref.plink.command, objref.file.plink.bed, objref.file.plink.bim, objref.file.plink.fam,
#                                                                                        select.snps = objref.snpdata.plink.bim[
#                                                                                                                         loading.idx, 2],
#                                                                                                                     select.subjects = select.subjects )
#
#         if (
#
#
#             class(plink.obj)== "try-error")
#             stop("Error in PLINK data extraction.");
#
#             objref.snpdata.local.idx =loading.idx;
#             objref.snpdata.local.fam =plink.obj.fam;
#             objref.snpdata.local.map =plink.obj.map;
#             objref.snpdata.local.snpmat =plink.obj.genotypes;
#
#
#             get_local_snpmat < -function(snp.k)
#
#             k.idx =which(snp.k == objref.snpdata.local.idx);
#
#
#         return (as.numeric(objref.snpdata.local.snpmat[, k.idx]) - 1 );
#
#
#         snpmat =matrix(NA, nrow=objref.n.ind.used, ncol = NROW(snp.idx.new));
#
#         for (i in 1:NROW(snp.idx.new))
#
#         if (!check_local_snpmat(snp_idx.new[i]) )
#         load_local_snpmat(snp.idx.new[i])
#
#         snpmat[, i] =get_local_snpmat(snp.idx.new[i])
#         snpmat[snpmat[, i] == -1, i] =NA;
#         if (mean(snpmat[, i], na.rm=T) / 2 > 0.5 )
#         snpmat[, i] =2 - snpmat[, i];
#
#
#         snpmat =snpmat[, order(snp.idx.ord), drop = F];
#
#         colnames(snpmat) = objref.snpdata.plink.bim[snp.idx, 2];
#         rownames(snpmat) = objref.ids.used;
#
#         if (impute)
#
#         snpmat.imp =impute_simple_snp(snpmat);
#         return (list(snpmat=snpmat.imp.snpmat, NMISS=snpmat.imp.NMISS, MAF=snpmat.imp.MAF))
#
#
#         else
#
#         NMISS =unlist(apply(snpmat, 2, function(snp.k)
#       length(which( is.na(snp.k)));}));
#         return (list(snpmat=snpmat, NMISS=NMISS, MAF=colMeans(snpmat, na.rm=T) / 2));
#
#
#
#
# def convert_simpe_to_plink(snp_info, snp_mat, snp_file_base):
#  # snp,mat : 0/1/2/NA
#  # PLINK raw data: 1/2/3==> AA,AB,BB, 0==>NA
#  snp_mat= (snp_mat+1).T
#  snp.mat[where(snp_mat)is.na(snp.mat)] =0
#
#  sub.name =colnames(snp.mat)
#  snp.name =rownames(snp.mat)
#
#  ###snps
#  dim.snps =dim(snp.mat)
#
#  snps =as.raw( as.matrix(snp.mat ) )
#  snps =array(snps, dim=dim.snps)
#  colnames(snps) =sub.name
#  rownames(snps) =snp.name
#  class(snps) ="SnpMatrix"
#
#  r =write.plink( file.base=snp.file.base, snp.major = F, snps=t(snps),
#       id=sub.name,
#       father=rep(0,dim.snps[2]),
#       mother=rep(0,dim.snps[2]),
#       sex=rep(0,dim.snps[2]),
#       phenotype = rep(-9,dim.snps[2]),
#    chromosome = as.character(snp.info[,1]),
#    genetic.distance = as.numeric(snp.info[,3]),
#    position= as.numeric(snp.info[,4]),
#    allele.1 = as.character(snp.info[,5]),,
#    allele.2 = as.character(snp.info[,6]),
#    na.code=0)
#
#   print(" Genotype files have been converted into PLINK binary format(bed/bim/fam)")
#
#  return(list(file_plink_bed = paste(snp.file.base, ".bed", sep=""),
#          file_plink_bim = paste(snp.file.base, ".bim", sep=""),
#          file_plink_fam = paste(snp.file.base, ".fam", sep="")))
#
# #
# # #get Identity-By-State matrix using PLINK command
# # fg_plink_getPCA( objref, plink.path )
# #
# #  plink.out.pca =tempfile( fileext=".pca" )
# #  t0 =plink_command( plink.path,
# #          c ( "--bed ", objref.reader.file_plink_bed,
# #              "--bim ", objref.reader.file_plink_bim,
# #              "--fam ", objref.reader.file_plink_fam,
# #              "--pca --out ", plink.out.pca)  )
# #
# #  tb =try( read.table(paste( plink.out.pca, "eigenvec", sep=".") ))
# #  if (class(tb)=="try-error")
# #
# #   show(t0)
# #   stop("Failed to call PLINK.")
# #
# #
# #     colnames(tb) =c("FID","IID", paste0("PCA", 1:(NCOL(tb)-2) ) )
# #  unlink( plink.out.pca )
# #
# #     return(tb)
# #
