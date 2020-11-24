# import os
# import  Simulation.FgObj as FgObj
# from pandas_plink  import  read_plink
# import  FgDmObject as FgDmObj
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
# def plink_checkFiles( refobj, verbose=False):
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
#      # if( all( file.exists( refobj.file_plink_bed,  refobj.file_plink_bim, refobj.file_plink_fam ) ) )
#      #
#      #      refobj.snpdata.plink.fam =read.table(refobj.file_plink_fam, header=F, stringsAsFactors=F)
#      #  colnames(refobj.snpdata.plink.fam) =c("FID", "IID", "PID", "MID", "SEX", "PHE")
#      #
#      #  refobj.snpdata.plink.bim =read.table(refobj.file_plink_bim, header=F, stringsAsFactors=F)
#      #  colnames(refobj.snpdata.plink.bim) =c("CHR", "SNP", "CM", "POS", "A1", "A2")
#      #
#      #  if(refobj.chromosome != -1)
#      #
#      #   snp.idx =which(refobj.snpdata.plink.bim[,1] %in% refobj.chromosome)
#      #   refobj.snpdata.plink.bim =refobj.snpdata.plink.bim[snp.idx, ]
#      #
#      #
#      #  refobj.n.snp        =NROW(refobj.snpdata.plink.bim)
#      #  refobj.n_ind_total  =NROW(refobj.snpdata.plink.fam)
#      #  refobj.n.ind.used   =refobj.n_ind_total
#      #  refobj.ids.used     =as.character(refobj.snpdata.plink.fam[,2])
#
#       #return refobj
#      #
#      # else:
#      #  # stop("PLInK files can not be found!")
#      #  pass
#
#
#      (objref, snp_idx=None, impute=False, allel=False )
#
#      if snp
#          snp.idx =1:objref.n.snp
#
#      snp.idx.ord =order( snp.idx, decreasing=F)
#      snp.idx.new =snp.idx[ snp.idx.ord ]
#
#      check_local_snpmat(snp.k)
#
#       if(is.None(objref.snpdata.local.idx))
#        return(False)
#
#       return(snp.k %in% objref.snpdata.local.idx)
#
#
#      load_local_snpmat(snp.k)
#
#       select.subjects =objref.snpdata.plink.fam[ match( objref.ids.used, objref.snpdata.plink.fam[,2] ), c(1,2)]
#       loading.idx =seq( snp.k, ifelse(snp.k + objref.snp.blocksize > objref.n.snp, objref.n.snp, snp.k + objref.snp.blocksize), 1)
#       plink.obj =plink.cmd.select( objref.plink.command, objref.file_plink_bed, objref.file_plink_bim, objref.file_plink_fam,
#           select.snps = objref.snpdata.plink.bim[loading.idx, 2],
#           select.subjects = select.subjects )
#
#       if(class(plink.obj)=="try-error")
#        stop("Error in PLINK data extraction.")
#
#       objref.snpdata.local.idx   =loading.idx
#       objref.snpdata.local.fam   =plink.obj.fam
#       objref.snpdata.local.map   =plink.obj.map
#       objref.snpdata.local.snpmat=plink.obj.genotypes
#
#
#      get_local_snpmat(snp.k)
#
#       k.idx =which(snp.k==objref.snpdata.local.idx)
#       return(as.numeric(objref.snpdata.local.snpmat[,k.idx]) -1 )
#
#
#      snpmat =matrix(NA, nrow = objref.n.ind.used, ncol=NROW(snp.idx.new))
#
#      for(i in 1:NROW(snp.idx.new))
#
#       if ( !check_local_snpmat( snp.idx.new[i] ) )
#        load_local_snpmat( snp.idx.new[i] )
#
#       snpmat[,i] =get_local_snpmat( snp.idx.new[i] )
#       snpmat[ snpmat[,i]==-1 , i ]  =NA
#       if (mean(snpmat[,i], na.rm=T)/2>0.5 )
#        snpmat[,i] =2 - snpmat[,i]
#
#
#      snpmat =snpmat[, order(snp.idx.ord), drop=F]
#
#      colnames(snpmat) = objref.snpdata.plink.bim[snp.idx,2]
#      rownames(snpmat) = objref.ids.used
#
#      if( impute )
#
#       snpmat.imp =impute_simple_snp(snpmat)
#       return(list(snpmat = snpmat.imp.snpmat, NMISS=snpmat.imp.NMISS, MAF=snpmat.imp.MAF))
#
#
#      else
#
#       NMISS =unlist(apply(snpmat, 2, function(snp.k) length(which(is.na(snp.k)))))
#       return(list(snpmat = snpmat, NMISS=NMISS, MAF=colMeans(snpmat, na.rm=T)/2))
#
#
#
#     plink.get.snpinfo(objref, snp.idx)
#
#      # remove genetic distance at 3rd postion
#
#      if(is.None(snp.idx))
#       snp.idx =1:objref.n.snp
#
#      return( objref.snpdata.plink.bim[snp.idx,c(2,1,4,5,6), drop=F])
#
#
#     plink.get.individuals(objref)
#
#      return( objref.snpdata.plink.fam[,2])
#
#
#     plink.select.individuals(objref, ids.used)
#
#      objref.ids.used =as.character(ids.used)
#
#
#     plink.get.used.individuals(objref)
#
#      return(objref.ids.used)
#
#
#     plink.get.snpindex =function( objref, snp.names )
#
#      return( match(snp.names, objref.snpdata.plink.bim[,2]) )
#
#
#     plink.shrink (objref)
#
#      objref.snpdata.plink.snpmat =None
#      objref.snpdata.plink.fam =None
#      objref.snpdata.plink.map =None
#
#
#     test_plink_func(plink)
#
#      n.snp =NCOL( plink.genotypes )
#
#       print("...Get 100 SNPs")
#      r.snp =get_plink_subsnp(plink, 1:100)
#
#       print("...Get 1000 SNPs")
#      r.snp =get_plink_subsnp(plink, 1:1000)
#
#       print("...Get 10000 SNPs")
#      r.snp =get_plink_subsnp(plink, 1:10000)
#
#       print("...Get 50000 SNPs")
#      r.snp =get_plink_subsnp(plink, 1:50000)
#
#       print("...Get 100000 SNPs")
#      r.snp =get_plink_subsnp( plink, sample(n.snp)[1:100000] )
#
#
#     get_plink_subsnp(snp.mat, snp.set.idx)
#
#      s.mat =as( snp.mat.genotypes[, snp.set.idx, drop=F ], "numeric")
#      snp.imp <-c()
#      snp.maf =c()
#      snp.names =c()
#
#      f.impute(s.mat.i )
#
#       s.miss =which( is.na(s.mat.i) )
#
#       if (length(s.miss)>0)
#
#        n.AA =length( which( s.mat.i == 0 ) )
#        n.Aa =length( which( s.mat.i == 1 ) )
#        n.aa =length( which( s.mat.i == 2 ) )
#        n.s  =n.AA + n.Aa + n.aa
#
#        r.miss =runif( length(s.miss) )
#        r.snp  =rep(2, length(s.miss))
#        r.snp [ r.miss <= n.AA/n.s ]<-0
#        r.snp [ r.miss <= (n.AA + n.Aa)/n.s ]<-1
#        s.mat.i[s.miss] =r.snp
#
#
#       if (mean(s.mat.i)/2>0.5) s.mat.i =2 - s.mat.i
#
#       return( s.mat.i )
#
#
#
#      snp.imp =t( apply(s.mat, 2, f.impute) )
#      snp.maf =rowMeans(snp.imp) /2
#
#      map =snp.mat.map[snp.set.idx, ,drop=F]
#      rownames(snp.imp) <-rownames( map )
#
#      return(list(maf=snp.maf, snp=snp.imp, info=map[,c(2,1,4)]) )
#
#
#     impute_simple_snp =function( snpmat )
#
#      f_imputed =function( snp )
#
#       s.miss =which( is.na( snp ) )
#       if ( length(s.miss)>0 )
#
#        n.AA =length(which( snp == 0 ) )
#        n.Aa =length(which( snp == 1 ) )
#        n.aa =length(which( snp == 2 ) )
#        n.s  =n.AA + n.Aa + n.aa
#
#        r.miss =runif( length(s.miss) )
#        r.snp  =rep(2, length( s.miss ))
#        r.snp [ r.miss <= n.AA/n.s ] =0
#        r.snp [ r.miss <= (n.AA + n.Aa)/n.s ] <-1
#        snp [ s.miss ] =r.snp
#
#
#       if (mean(  snp )/2>0.5)
#        snp =2 - snp
#
#       return(list(snp=snp, NMISS=length(s.miss), MAF=mean(  snp )/2 ))
#
#
#      total_miss =length(which(is.na(snpmat)))
#
#      r.imp  =apply(snpmat, 2, f_imputed)
#      snpmat =do.call("cbind", lapply(1:NCOL(snpmat), function(i)return(r.imp[[i]].snp)))
#      NMISS  =do.call("unlist", lapply(1:NCOL(snpmat), function(i)return(r.imp[[i]].NMISS)))
#      MAF    =do.call("unlist", lapply(1:NCOL(snpmat), function(i)return(r.imp[[i]].MAF)))
#
#
#       print("* Missing SNPs are imputed(", total_miss, "SNPs).")
#
#      return(list(snpmat=snpmat, NMISS=NMISS, MAF=MAF))
#
#
#     plink.cmd.load.chromosome(plink.command,file_plink_bed, file_plink_bim, file_plink_fam, chr)
#
#      tmp =tempfile()
#      str.cmd =paste( plink.command, "--noweb",
#          "--bed", file_plink_bed,
#          "--bim", file_plink_bim,
#          "--fam", file_plink_fam,
#          "--chr", chr,
#          "--make-bed",
#          "--out", tmp,
#          sep=" ")
#
#      t =try(system( str.cmd, ignore.stdout=True, ignore.stderr=True) )
#      if(class(t)=="try-error")
#
#        print("! Error in PLINK command.")
#       return(list(error=T, err.info="Error in PLINK command."))
#
#
#      plink =try( read.plink( paste(tmp, ".bed", sep=""),  paste(tmp, ".bim", sep=""), paste(tmp, ".fam", sep="")) )
#      if(class(plink) == "try-error")
#
#        print("! Package snpStats can not open PLINK data.")
#       return(list(error=T, err.info="Package snpStats can not open PLINK data."))
#
#
#      return(plink)
#
#
#
#     plink.cmd.select(plink.command, file_plink_bed, file_plink_bim, file_plink_fam, select.snps, select.subjects)
#
#      file.keep.snp =tempfile()
#      file.keep.subject =tempfile()
#      write.table(select.snps, file=file.keep.snp, quote=F, row.names=F, col.names=F)
#      write.table(select.subjects, file=file.keep.subject, quote=F, row.names=F, col.names=F)
#
#      tmp =tempfile()
#      str.cmd =paste( plink.command, "--noweb",
#          "--bed", file_plink_bed,
#          "--bim", file_plink_bim,
#          "--fam", file_plink_fam,
#          "--keep", file.keep.subject,
#          "--extract", file.keep.snp,
#          "--make-bed",
#          "--out", tmp,
#          sep=" ")
#
#      t =try(system( str.cmd, ignore.stdout=True, ignore.stderr=False) )
#      if(class(t)=="try-error")
#
#        print("! Error in PLINK command.")
#       return(list(error=T, err.info="Error in PLINK command."))
#
#
#      plink =try( read.plink( paste(tmp, ".bed", sep=""),  paste(tmp, ".bim", sep=""), paste(tmp, ".fam", sep="")) )
#      if(class(plink) == "try-error")
#
#        print("! Package snpStats can not open PLINK data.")
#       return(list(error=T, err.info="Package snpStats can not open PLINK data."))
#
#
#      unlink(c(file.keep.snp, file.keep.subject, paste(tmp, ".bed", sep=""),  paste(tmp, ".bim", sep=""), paste(tmp, ".fam", sep="") ))
#
#      return(plink)
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
