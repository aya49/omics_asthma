## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result")


# asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
eqtl_dir = paste0(stat_dir,"/eqtl"); dir.create(eqtl_dir, showWarnings=F)



## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("data.table")
libr("MatrixEQTL")
libr("foreach")
libr("doMC")
libr("stringr")
libr("qdapTools") # make dummy variables
libr("Matrix")



## options
no_cores = 15#detectCores()-3
registerDoMC(no_cores)

overwrite = F
writecsv = T #write results as csv on top of Rdata?

writetra = F#write all associations not only local ones

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis

id_col = "id"
class_col = "response"
categorical = T # is class column categorical?
interested_cols = c("response")
# interested_cols = c("age","bmi","sex","centre","batch","race","response") 
# interested_cont_cols = ""

caucasians_only = T

f1_bins = c("","01","12") # 01 make all 2s 1, 12 make all 0s 1
# split_f1_col = "time" #there are
useModels = c("modelLINEAR", "modelANOVA", "modelLINEAR_CROSS")

# plotting size
width = 800
height = 600




# parameters: data1 (genotype!!), data2 (rna!!),
#  
# # Only associations significant at this level will be saved
# pvalthres_cis = 2e-2; # numeric. Significance threshold for all/distant tests
# pvalthres_tra = 1e-2; # numeric. Same as pvOutputThreshold, but for local eQTLs.
# 
# # Model:
# # int OR
# # modelLINEAR to model the effect of the genotype as additive linear and test for its significance using t-statistic
# # modelANOVA to treat genotype as a categorical variables and use ANOVA model and test for its significance using F-test. The default number of ANOVA categories is 3. Set otherwise like this: options(MatrixEQTL.ANOVA.categories=4)
# # modelLINEAR_CROSS to add a new term to the model equal to the product of genotype and the last covariate; the significance of this term is then tested using t-statistic
# useModel = modelLINEAR
# 
# # Error covariance matrix
# # Set to numeric() for identity.
# errorCovariance = numeric();
# # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
# 
# # Distance for local gene-SNP pairs
# #  numeric. SNP-gene pairs within this distance are considered local. 
# #  The distance is measured from the nearest end of the gene. 
# #  SNPs within a gene are always considered local.
# cisDist = 1e6;
feat_types = list(c("genotype", "rnaseqgenes.pre", 2e-6, 1e-6, 1e6),
                  c("genotype", "rnaseqgenes.post", 2e-6, 1e-6, 1e6),
                  c("genotype", "rnaseqgenes.diff", 2e-6, 1e-6, 1e6),
                  c("genotype", "rnaseqisoforms.pre", 2e-6, 1e-6, 1e6),
                  c("genotype", "rnaseqisoforms.post", 2e-6, 1e-6, 1e6),
                  c("genotype", "rnaseqisoforms.diff", 2e-6, 1e-6, 1e6)
                  # c("rnaseqgenes.pre", "rnaseqgenes.post", 2e-2, 1e-2, "modelLINEAR", "NA", 1)
)







start = Sys.time()

meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))

for (feat_type in feat_types) {
  
  # name parameters
  feat_type1 = feat_type[1]
  feat_type2 = feat_type[2]
  pvalthres_cis = as.numeric(feat_type[3])
  pvalthres_tra = as.numeric(feat_type[4])
  # useModel = feat_type[5]
  # errorCovariance = numeric()
  cisDist = as.numeric(feat_type[5])
  
  # load genotype data
  f1_m0 = get(load(paste0(feat_dir,"/",feat_type1,".Rdata")))
  f1_meta_col0 = get(load(paste0(meta_col_dir,"-",str_split(feat_type1,"[.]")[[1]][1],".Rdata")))
  f1_meta_col0 = f1_meta_col0[match(colnames(f1_m0),f1_meta_col0[,id_col]),]
  # if (colnames(f1_m0)[1]%in%f1_meta_file0[,id_col]) f1_m0 = t(f1_m0)
  
  # load continuous data
  f2_m0 = get(load(paste0(feat_dir,"/",feat_type2,".Rdata")))
  f2_meta_col0 = get(load(paste0(meta_col_dir,"-",str_split(feat_type2,"[.]")[[1]][1],".Rdata")))
  f2_meta_col0 = f2_meta_col0[match(colnames(f2_m0),f2_meta_col0[,id_col]),]
  # if (colnames(f1_m0)[1]%in%f1_meta_file0[,id_col]) f1_m0 = t(f1_m0)
  
  # get col genotype & rna features
  f1_meta_col_ind = apply(f1_meta_col0[,c("dbSNP","chromosome","pos_phys")],1,function(x) !any(is.na(x)))
  f1_meta_col = f1_meta_col0[f1_meta_col_ind,]
  f1_good_col_na = f1_meta_col0$dbSNP #SNP ONLY
  
  f2_meta_col_ind = apply(f2_meta_col0[,c("symbol","chr","start","end")],1,function(x) !any(is.na(x)))
  f2_meta_col = f2_meta_col0[f2_meta_col_ind,]
  f2_good_col_na = f2_meta_col0$symbol #GENE ONLY
  
  # get row files/samples & covariate
  samples_to_include = intersect(rownames(f2_m0), rownames(f1_m0))
  meta_file = 
    meta_file0[!is.na(meta_file0[,class_col]) & 
                 meta_file0[,id_col] %in% samples_to_include &
                 ifelse(caucasians_only, grepl("Caucasian",meta_file0[,"race"]),rep(T,nrow(meta_file0))),]
  rownames(meta_file) = meta_file[,id_col]
  
  
  # trim matrices
  f1_m = f1_m0[meta_file[,id_col], f1_meta_col_ind]
  f2_m = f2_m0[meta_file[,id_col], f2_meta_col_ind]
  
  
  
  ## prepare meta_col
  #  data.frame with columns snpid (Snp_01), chr (1), pos (725123)
  f1_pos = data.frame(snpid=f1_meta_col$dbSNP, chr=paste0("chr",f1_meta_col$chromosome), pos=f1_meta_col$pos_phys)
  # levels(f1_pos$chr) = paste("chr", c(1:22, "X", "Y", "M"), sep="") #convert to bioconductor format
  #  data.frame with columns geneid (Gene_01), chr (1), left (721289), right (731289)
  cgene_col = ifelse("gene"%in%colnames(f2_meta_col),"gene","id")
  f2_pos = data.frame(geneid=f2_meta_col$symbol, chr=paste0("chr",f2_meta_col$chr), 
                      left=f2_meta_col$start, right=f2_meta_col$end)
  
  ## prepare meta_file
  cvrt00 = meta_file[,interested_cols,drop=F]
  cvrt_num_ind = apply(cvrt00,2, function(x)
    all(!is.na(as.numeric(x[!is.na(x) & x!=""]))) )
  cvrt0 = as.data.frame(lapply(append(which(cvrt_num_ind),which(!cvrt_num_ind)), function(x) {
    if (cvrt_num_ind[x]) return(as.numeric(cvrt00[,x]))
    if (length(unique(cvrt00[,x]))>2) return(mtabulate(cvrt00[,x]))
    return(as.numeric(factor(cvrt00[,x]))-1)
  }))
  rownames(cvrt0) = rownames(meta_file)
  # colnames(cvrt0)[1:sum(cvrt_num_ind)] = colnames(cvrt00)[cvrt_num_ind]
  # cvrt1 = model.matrix(~cvrt0)
  
  cvrt = SlicedData$new()
  cvrt$CreateFromMatrix(t(as.numeric(as.matrix(cvrt0))))
  
  for (f1_bin in f1_bins) {
    for (useModel_ in useModels) {
      
      # # Output file name
      # output_file_name_cis = tempfile();
      # output_file_name_tra = tempfile();
      eqtl_cis_dir = paste0(eqtl_dir,"/", 
                            feat_type1,ifelse(f1_bin!="",".",""),f1_bin,"-", feat_type2,
                            "_cov-", paste(interested_cols,collapse="."), 
                            "_cis", ".Rdata")
      eqtl_tra_dir = NULL; if (writetra) eqtl_tra_dir = gsub("cis_","tra_",eqtl_cis_dir)
      
      # to overwrite or not
      if1 = file.exists(eqtl_cis_dir)
      if2 = T; if (!is.null(eqtl_tra_dir)) if2 = file.exists(eqtl_tra_dir)
      if ((if1 & if2 & !overwrite) | 
          (f1_bin!="" & useModel_=="modelLINEAR_CROSS")) next()
      print(paste0(eqtl_cis_dir))
      
      # set model
      useModel = NULL
      if (useModel_=="modelLINEAR") useModel = modelLINEAR # genotype is assumed to have only additive effect on expression
      if (useModel_=="modelANOVA") useModel = modelANOVA # assume genotype to have both additive and dominant effects (ANOVA model). In this case genotype data set musts take at most 3 distinct values (i.e. 0/1/2/NA)
      if (useModel_=="modelLINEAR_CROSS") useModel = modelLINEAR_CROSS
      
      
      
      ## prepare m genotyping and rnaseq data as SlicedData for input into MatrixEQTL()
      # to make genotype binary or not
      if (f1_bin=="01") f2_m[f2_m==2] = 1
      if (f1_bin=="12") f2_m[f2_m==0] = 1
      
      #  Can be real-valued for linear models 
      #  and must take at most 3 distinct values for ANOVA 
      #  unless the number of ANOVA categories is set to a higher number (see useModel parameter).
      #  Must have matching columns
      colnames(f1_m) = f1_pos$snpid
      f1_sd = SlicedData$new()
      f1_sd$CreateFromMatrix(t(f1_m))
      
      colnames(f2_m) = f2_pos$geneid
      f2_sd = SlicedData$new()
      f2_sd$CreateFromMatrix(t(f2_m))
      # sd$ResliceCombined(sliceSize = 2L) # Slice it in pieces of 2 rows
      # length(sd) # Show the number of slices (equivalent function calls)
      # sd$IsCombined() # Is it all in one slice? (No)
      # colnames(sd) # Show the column names (equivalent function calls)
      # rownames(sd) # Show all row names (equivalent function calls)
      # print(sd[[2]]) # Print the second slice
      # sd$ColumnSubsample( c(1,3,4) ) # Reorder and subset columns
      # sd$RowReorder( c(3,1) ) # Reorder and subset rows  
      # sd$FindRow("row1") # Find the row with name "row1" (it is second in the first slice)
      # show(sd) # Show the detail of the object (one slice again)
      
      
      ## Run the analysis
      
      try({
        # eqtl_cis_dir = tempfile()
        me = Matrix_eQTL_main(
          snps = f1_sd,
          gene = f2_sd,
          cvrt = cvrt, # SlicedData object with additional covariates. Can be an empty SlicedData object in case of no covariates. The constant is always included in the model and would cause an error if included in cvrt. The order of columns must match those in snps and gene.
          
          output_file_name = eqtl_tra_dir, # significant associations (all significant associations if pvOutputThreshold=0 or only distant if pvOutputThreshold>0). If the file with this name exists, it is overwritten.
          output_file_name.cis = eqtl_cis_dir, #output_file_name_cis=tempfile(); significant local associations
          
          pvOutputThreshold = pvalthres_tra, 
          pvOutputThreshold.cis = pvalthres_cis,
          
          useModel = useModel, 
          # errorCovariance = errorCovariance, # numeric. The error covariance matrix. Use numeric() for homoskedastic independent errors.
          verbose = F, 
          
          snpspos = f1_pos,
          genepos = f2_pos, 
          cisDist = cisDist,
          pvalue.hist = "qqplot", # logical, numerical, or "qqplot" (faster if false); To record information for a histogram set pvalue.hist to the desired number of bins of equal size. Finally, pvalue.hist can also be set to a custom set of bin edges.
          min.pv.by.genesnp = T, # record the minimum p-value for each SNP and each gene in the returned object. The minimum p-values are recorded even if if they are above the corresponding thresholds of pvOutputThreshold and pvOutputThreshold.cis (faster if false)
          noFDRsaveMemory = F # save significant gene-SNP pairs directly to the output files, reduce memory footprint and skip FDR calculation. The eQTLs are not recorded
        ) 
        # unlink(output_file_name_tra);
        # unlink(eqtl_cis_dir);
        
        save(me, file=gsub("_cis","",eqtl_cis_dir))
        
        ## Results:
        cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
        cat('Detected local eQTLs:', '\n');
        write.csv(me$cis$eqtls, file=gsub(".Rdata",".csv",eqtl_cis_dir))
        # write.csv(me$trans$eqtls, file=gsub(".Rdata",".csv",gsub("cis","tra",eqtl_cis_dir)))
        cat('Detected distant eQTLs:', '\n');
        show(me$trans$eqtls)
        
        ## Plot the histogram of local and distant p-values
        png(gsub(".Rdata","_qq.png",eqtl_cis_dir), width=width, height=height)
        plot(me)
        graphics.off()
        
      })
    } # useModel
  } # f1_bin
} # feat_type


time_output(start)


















