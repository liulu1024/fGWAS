# import os
# import  pandas as pd
# import  FgDmObject as fgDmObj
# def create_simple(file_simple_snp, *args):
#    if file_simple_snp is None:
#       os.exit("! file.simple.snp  must be assigned with the valid values.")
#
#    objref = fgDmObj.FgDmSimple(file_simple_snp,pd.DataFrame(),list())
#
#    # tb_gen = read.table( objref.file_simple_snp, header=T );
#    # objref.n.snp    = NROW(tb_gen)
#    # objref.n.ind.total    = NCOL(tb_gen)-5
#    # objref.n.ind.used    = NCOL(tb_gen)-5
#    # objref.ids.used    = colnames(tb_gen)[-c(1:5)];
#    #
#    # load_simple_snp(objref, options.verbose);
#
#    return objref
#
# def simple_shrink(objref):
#    return objref
#
# def  covert_geno(allel_info, allel_geno):
#       snpB    = as.character(unlist(allel.info[5]));
#       snpA    = as.character(unlist(allel.info[4]));
#       QQ2    = paste(snpB, snpB, sep="");
#       qq0    = paste(snpA, snpA, sep="");
#       Qq1    = c(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;
#
#       d.g    = rep( NA, length(allel.geno) );
#       d.g[which(allel.geno == QQ2)]   =2;
#       d.g[which(allel.geno == qq0)]   =0;
#       d.g[which(allel.geno == Qq1[1])]   =1;
#       d.g[which(allel.geno == Qq1[2])]   =1;
#
#       if( mean(d.g, na.rm=T)/2 > 0.5)
#          d.g    = 2 - d.g;
#
#       return d_g
#
# def simple_get_snpmat(objref, snp_idx=None, impute=False, allel=False):
#
#    allel_mat = None
#    allel_info = None
#    if snp_idx is None:
#       allel_mat =  objref.snpdata.snp_mat[:,objref.ids_used]
#       allel_info  = objref.snpdata.snp_info
#
#    else:
#       allel_mat  = objref.snpdata.snp_mat[ snp_idx, objref.ids_used]
#       allel_info = objref.snpdata.snp_info[ snp.idx, ]
#
#    ##BB/AB/AA/A./B.==>2/1/0/NA
#
#
#    ## if data element is AA/AB/BB
#    if(is.character(allel.mat[1,1]) || is.factor(allel.mat[1,1]) )
#    {
#       snpmat    = do.call( "cbind", lapply(1:NROW(allel.mat), function(i){
#       return( covert_geno ( allel.info[i,], as.character(unlist(allel.mat[i,])) ) );
#       }));
#    }
#    ## if data element is 0/1/2/NA
#    else
#    {
#       snpmat    = t(allel.mat);
#    }
#
#    rownames(snpmat)    = colnames(allel.mat);
#    colnames(snpmat)    = rownames(allel.mat);
#
#    if(!allel)
#    {
#       if(impute)
#       {
#          snpmat.imp    = impute_simple_snp( snpmat );
#          return(list(snpmat = snpmat.imp.snpmat, NMISS=snpmat.imp.NMISS, MAF=snpmat.imp.MAF))
#       }
#       else
#       {
#          NMISS    = unlist(apply(snpmat, 2, function(snp.i){length(which(is.na(snp.i)))}));
#          return(list(snpmat = snpmat, NMISS=NMISS, MAF=colMeans(snpmat, na.rm=T)/2));
#       }
#    }
#    else
#       return(list(snpmat = allel.mat, NMISS=NA, MAF=NA))
# }
#
# def simple.get.snpinfo(objref, snp.idx=None ):
#    if (is.None(snp.idx))
#       return(objref.snpdata.snp.info)
#    else
#       return(objref.snpdata.snp.info[snp.idx,,drop=F]);
#
# def simple_get_individuals(objref):
#    return( colnames(objref.snpdata.snp.mat) );
#
# def simple_select_individuals(objref, ids_used):
#    objref.ids.used    = ids.used;
#
#    load_simple_snp(objref);
#
#    ids.all    = colnames(objref.snpdata.snp.mat)
#    if(any(is.na(match(ids.used, ids.all))))
#       stop("Some IDs are not matached in the SNP data file");
#
#    objref.n.ind.used    = length(ids.used);
#    return(objref);
#
# def simple_get_used_individuals(objref):
#    return( objref.ids.used );
#
# def simple_get_snpindex( objref, snp_names):
#    return( match(snp.names, rownames(objref.snpdata.snp.info) ) );
#
#
# #scaffoldId, Loci, RefBase, AltBase, P1, P2, P3,....
# def load_simple_snp(objref, verbose=True):
#    if not objref.file_simple_snp =="":
#       tb_gen    = read.table( objref.file_simple_snp, header=T, stringsAsFactors=F, colClasses=c(SNP="character", CHR="character", RefAllele="character", AltAllele="character"))
#    else:
#       tb_gen    = objref.rawdata;
#
#    if(verbose)
#       cat("gen.table[", dim(tb_gen), "]\n");
#
#    snp_info    = tb_gen[, c(1:5)]
#    snp_mat     = tb_gen[, -c(1:5)]
#    rownames(snp_info)    = snp_info[,1];
#    objref.snpdata.snp.info    = snp_info
#
#    ids.idx    = match( objref.ids.used, colnames(snp_mat) );
#    objref.snpdata.snp.mat    = snp_mat[, ids.idx];
#
# def fg_load_simple( file_simple_snp, options ):
#    if ( missing(file.simple.snp) )
#       stop("! file.simple.snp must be assigned with the valid file name..");
#
#    if(options.verbose)
#    {
#       cat( "[ Loading Simple ] \n");
#       cat( " Checking the parameters ......\n");
#       cat("* SIMPLE SNP File = ",  file.simple.snp, "\n");
#    }
#
#    if(!file.exists( file.simple.snp))
#       stop("Failed to open Simple data files.");
#
#    params    = list( file.simple.snp = file.simple.snp );
#
#    obj.gen    = create.simple.obj ( file.simple.snp, options );
#
#    ret    = list( reader=obj.gen, options=options, params=params );
#    ret.n.snp    = obj.gen.n.snp
#    ret.n.ind.total    = obj.gen.n.ind.total
#    ret.n.ind.used    = obj.gen.n.ind.used
#    class(ret)    = "fgwas.gen.obj";
#
#    return(ret);
#
# def plink_command(plink.path, plink.parms):
#    t1    = try(system(paste(c(plink.path, plink.parms), collapse=" "), intern = TRUE))
#    ## show(t1);
#
#    return(t1)
#
# #get Identity-By-State matrix using PLINK command
# def fg_simple_getPCA(objref, plink_path):
#    snp_mat    = objref.reader.get_snpmat( None, impute=False, allel=False).snpmat;
#    snp_info    = objref.reader.get_snpinfo(None );
# colnames(snp.info)    = c("SNP","CHR","POS","A1","A2");
#
#    snp.mat    = snp.mat[, with(snp.info, order(CHR, POS))]
#    snp.info    = snp.info[ with(snp.info, order(CHR, POS)),]
#
#    #change chromosome string to number;
#    snp.info.CHR=as.numeric(factor(snp.info.CHR))
#
#    snp.file.base    = tempfile(fileext = ".plink")
#    r    = convert_simpe_to_plink( data.frame(snp.info[,c(2,1)], 0, snp.info[,c(3:5)]),  snp.mat, snp.file.base);
#
#    plink.out.pca    = paste(snp.file.base, "pca", sep=".")
#    t0    = plink_command( plink.path, c ( "--bfile ", snp.file.base, "--pca --out ", plink.out.pca)  ) ;
#
#    tb    = try(read.table(paste(plink.out.pca, "eigenvec", sep=".")));
#    if (class(tb)=="try-error")
#    {
#       show(t0);
#       stop("Failed to call PLINK.");
#    }
#
#    unlink( paste(snp.file.base, c("bim", "bed", "fam", "pca.eigenvec"), sep=".") );
#
#    rownames(tb)    = tb[,1];
#    tb    = tb[, -c(1,2)];
#        colnames(tb)    = paste0("PCA", 1:NCOL(tb) );
#
#    return(tb);
# }
