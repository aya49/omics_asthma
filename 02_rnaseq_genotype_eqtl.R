## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614



## root directory
root = "~/projects/asthma"
setwd(root)

result0_dir = paste0(root, "/result")
result_dir = list.dirs(result0_dir, recursive=F)


# asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")


## output directory
stat_dir = paste0(result0_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
eqtl_cis_dir = paste0(stat_dir,"/eqtl")



## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("data.table")
libr("MatrixEQTL")
libr("foreach")
libr("doMC")
libr("stringr")
libr("Matrix")



## options
no_cores = 15#detectCores()-3
registerDoMC(no_cores)

overwrite = F
writecsv = T #write results as csv on top of Rdata?

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis

id_col = "fileName"
incommon_col = "sample" #meta_files columns gt and rs has in common
cid_col = "id"
class_col = "response"
categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""

gt_bins = c("","01","12") # 01 make all 2s 1, 12 make all 0s 1
split_rs_col = "time" #there are

# plotting size
width = 800
height = 600

# Only associations significant at this level will be saved
pvalthres_cis = 2e-2; # numeric. Significance threshold for all/distant tests
pvalthres_tra = 1e-2; # numeric. Same as pvOutputThreshold, but for local eQTLs.

# Model:
# int OR
# modelLINEAR to model the effect of the genotype as additive linear and test for its significance using t-statistic
# modelANOVA to treat genotype as a categorical variables and use ANOVA model and test for its significance using F-test. The default number of ANOVA categories is 3. Set otherwise like this: options(MatrixEQTL.ANOVA.categories=4)
# modelLINEAR_CROSS to add a new term to the model equal to the product of genotype and the last covariate; the significance of this term is then tested using t-statistic
useModel = modelLINEAR

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
#  numeric. SNP-gene pairs within this distance are considered local. 
#  The distance is measured from the nearest end of the gene. 
#  SNPs within a gene are always considered local.
cisDist = 1e6;









start = Sys.time()


gt_inds = grep("genotyp",result_dir)
rs_inds = grep("RNAseq",result_dir)

for (gt_ind in gt_inds) {
  for (rs_ind in rs_inds) {
    
    # load data
    rs_meta_file0 = get(load(paste0(meta_file_dir[rs_ind],".Rdata")))
    rs_meta_col0 = get(load(paste0(meta_col_dir[rs_ind],".Rdata")))
    rs_m0_names = list.files(feat_dir[rs_ind],full.names=F)
    rs_m0_name = rs_m0_names[which.min(nchar(rs_m0_names))]
    rs_m0 = as.matrix(get(load(paste0(feat_dir[rs_ind],"/",rs_m0_name))))
    if (colnames(rs_m0)[1]%in%rs_meta_file0[,id_col]) rs_m0 = t(rs_m0)
    
    rs_good_col_inds = intersect(colnames(rs_m0), rs_meta_col0[!is.na(rs_meta_col0$symbol),cid_col])
    rs_m0 = rs_m0[,rs_good_col_inds]
    rs_meta_col0 = rs_meta_col0[match(rs_good_col_inds,rs_meta_col0[,cid_col]),]
    rs_meta_file0 = rs_meta_file0[match(rownames(rs_m0),rs_meta_file0[,id_col]),]
    
    gt_meta_file0 = get(load(paste0(meta_file_dir[gt_ind],".Rdata")))
    gt_meta_col0 = get(load(paste0(meta_col_dir[gt_ind],".Rdata")))
    gt_m0_names = list.files(feat_dir[gt_ind],full.names=F)
    gt_m0_name = gt_m0_names[which.min(nchar(gt_m0_names))]
    gt_m0 = as.matrix(get(load(paste0(feat_dir[gt_ind],"/",gt_m0_name))))
    if (colnames(gt_m0)[1]%in%gt_meta_file0[,id_col]) gt_m0 = t(gt_m0)
    
    gt_good_col_inds = colnames(gt_m0)[apply(gt_m0,2,function(x) min(table(x))>good_col & length(unique(x))>1)]
    gt_m0 = gt_m0[,gt_good_col_inds]
    gt_meta_col0 = gt_meta_col0[match(gt_good_col_inds,gt_meta_col0[,cid_col]),]
    gt_meta_file0 = gt_meta_file0[match(rownames(gt_m0),gt_meta_file0[,id_col]),]
    gt_good_row_inds = ! ((duplicated(gt_meta_file0$sample, fromLast=T) | duplicated(gt_meta_file0$sample)) & gt_meta_file0$kit=="Minikit")
    gt_m0 = gt_m0[gt_good_row_inds,]
    gt_meta_file0 = gt_meta_file0[gt_good_row_inds,]
    
    for (gt_bin in gt_bins) {
      for (rs_split in unique(rs_meta_file0[,split_rs_col])) {
        
        # # Output file name
        # output_file_name_cis = tempfile();
        # output_file_name_tra = tempfile();
        eqtl_cis_dir = paste0(eqtl_cis_dir,".tra_",gsub(".Rdata","",gt_m0_name),".",gt_ind,"_",gsub(".Rdata","",rs_m0_name),".","split-",split_rs_col,".",rs_split,".Rdata")
        eqtl_tra_dir = gsub("tra_","cis_",eqtl_cis_dir)
        
        if (file.exists(eqtl_cis_dir) & file.exists(eqtl_tra_dir) & !overwrite) next()
        print(paste0(eqtl_cis_dir,"\n",eqtl_tra_dir))
        
        ## trim matrices
        rs_meta_file = rs_meta_file0[rs_meta_file0[,split_rs_col]==rs_split,]
        
        samples_to_include = intersect(gt_meta_file0[,incommon_col],rs_meta_file[,incommon_col])
        rs_meta_file = rs_meta_file[rs_meta_file[,incommon_col]%in%samples_to_include,]
        gt_meta_file = gt_meta_file0[match(rs_meta_file[,incommon_col],gt_meta_file0[,incommon_col]),]
        
        rs_meta_col = rs_meta_col0
        rs_m = rs_m0[rs_meta_file[,id_col],]
        
        gt_meta_file = gt_meta_file0
        gt_meta_col = gt_meta_col0
        gt_m = gt_m0
        if (gt_bin=="01") gt_m[gt_m==2] = 1
        if (gt_bin=="12") gt_m[gt_m==0] = 1
        
        
        # prepare SNP and gene positions
        #  data.frame with columns snpid (Snp_01), chr (1), pos (725123)
        gt_pos = data.frame(snpid=gt_meta_col$dbSNP, chr=paste0("chr",gt_meta_col$chromosome), pos=gt_meta_col$pos_phys)
        # levels(gt_pos$chr) = paste("chr", c(1:22, "X", "Y", "M"), sep="") #convert to bioconductor format
        #  data.frame with columns geneid (Gene_01), chr (1), left (721289), right (731289)
        cgene_col = ifelse("gene"%in%colnames(rs_meta_col),"gene","id")
        rs_pos = data.frame(geneid=rs_meta_col$symbol, chr=paste0("chr",rs_meta_col$chr), 
                            left=rs_meta_col$start, right=rs_meta_col$end)
        
        
        # prepare genotyping and rnaseq data as SlicedData for input into MatrixEQTL()
        #  Can be real-valued for linear models 
        #  and must take at most 3 distinct values for ANOVA 
        #  unless the number of ANOVA categories is set to a higher number (see useModel parameter).
        #  Must have matching columns
        colnames(rs_m) = rs_meta_col$symbol
        rownames(rs_m) = rs_meta_file[,incommon_col]
        rs_sd = SlicedData$new()
        rs_sd$CreateFromMatrix(t(rs_m))
        
        colnames(gt_m) = gt_meta_col$dbSNP
        rownames(gt_m) = gt_meta_file[,incommon_col]
        gt_sd = SlicedData$new()
        gt_sd$CreateFromMatrix(t(gt_m[match(rownames(rs_m),rownames(gt_m)),]))
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
        
        me = Matrix_eQTL_main(
          snps = gt_sd,
          gene = rs_sd,
          # cvrt = cvrt, # SlicedData object with additional covariates. Can be an empty SlicedData object in case of no covariates. The constant is always included in the model and would cause an error if included in cvrt. The order of columns must match those in snps and gene.
          
          output_file_name = eqtl_cis_dir, # significant associations (all significant associations if pvOutputThreshold=0 or only distant if pvOutputThreshold>0). If the file with this name exists, it is overwritten.
          output_file_name.cis = eqtl_tra_dir, #output_file_name_cis=tempfile(); significant local associations
          
          pvOutputThreshold = pvalthres_tra, 
          pvOutputThreshold.cis = pvalthres_cis,
          
          useModel = useModel, 
          # errorCovariance = errorCovariance, # numeric. The error covariance matrix. Use numeric() for homoskedastic independent errors.
          verbose = T, 
          
          snpspos = gt_pos,
          genepos = rs_pos, 
          cisDist = cisDist,
          pvalue.hist = T, # logical, numerical, or "qqplot" (faster if false); To record information for a histogram set pvalue.hist to the desired number of bins of equal size. Finally, pvalue.hist can also be set to a custom set of bin edges.
          min.pv.by.genesnp = T, # record the minimum p-value for each SNP and each gene in the returned object. The minimum p-values are recorded even if if they are above the corresponding thresholds of pvOutputThreshold and pvOutputThreshold.cis (faster if false)
          noFDRsaveMemory = F # save significant gene-SNP pairs directly to the output files, reduce memory footprint and skip FDR calculation. The eQTLs are not recorded
        ) 
        # unlink(output_file_name_tra);
        # unlink(output_file_name_cis);
        
        ## Results:
        cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
        cat('Detected local eQTLs:', '\n');
        show(me$cis$eqtls)
        cat('Detected distant eQTLs:', '\n');
        show(me$trans$eqtls)
        
        ## Plot the histogram of local and distant p-values
        png(gsub(".Rdata",".csv",gsub("tra_","",eqtl_cis_dir)), width=width, height=height)
        plot(me)
        graphics.off()
        
      } #rs_split
    } # gt_bin
  } # rs_ind
} # gt_ind


time_out(start)


















