## input: features & meta_file
## output: eQTL etc
## aya43@sfu.ca
## created 20180614


# Cell Specific eQTL Analysis without Sorting Cells
# y ~ g*c # y=rnaseq, g=dna, c=cells
# NOTE: Yi = β0 + Xiβ1 + Ziβ2 + Wiβ3 + XiZiβ4 + XiWiβ5 + ZiWiβ6 + XiZiWiβ7 + ei
# Y ~ X + Z + W + X:Z + X:W + Z:W + X:Z:W
# Y ~ X * Z * W
# Y ~ (X + Z + W)^3


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/_func-classifiers.R"))
libr(append(pkgs(),c("qdaptools","MatrixEQTL"))) #parallelized

# no_cores = 1#detectCores()-2 #number of cores to use for parallel processing
# registerDoMC(no_cores)


## options
overwrite = T # overwrite if existing analysis eexists?
writecsv = T # write results as csv on top of Rdata?

writetra = F # write all associations not only local ones

categorical = T # is class column categorical?
# interested_cols = c("response")
interested_cols = c("response","sex")
# interested_cont_cols = ""
scale_cont = F #if values are continuous, scale?
caucasians_only = T

coln = list(c("dbSNP","chromosome","pos_phys"), c("symbol","chr","start","end"))
colnum = c("pos_phys","chr","start","end")
colstr = c("dbSNP","symbol")

# split_f1_col = "time" #there are
useModels = c("modelLINEAR", "modelANOVA", "modelLINEAR_CROSS")

# plotting size
width = 800
height = 600

# eqtl cut-offs
pvalthres_cis = .005
pvalthres_tra = .0005
# useModel = feat_type_set[5]
# errorCovariance = numeric()
cisDist = 1e6 #10000#


# parameters: data1 (dna!!), data2 (rna!!),
#
# # Only associations significant at this level will be saved
# pvalthres_cis = 2e-2; # numeric. Significance threshold for all/distant tests
# pvalthres_tra = 1e-2; # numeric. Same as pvOutputThreshold, but for local eQTLs.
#
# # Model:
# # int OR
# # modelLINEAR to model the effect of the dna as additive linear and test for its significance using t-statistic
# # modelANOVA to treat dna as a categorical variables and use ANOVA model and test for its significance using F-test. The default number of ANOVA categories is 3. Set results/enrichrwise like this: options(MatrixEQTL.ANOVA.categories=4)
# # modelLINEAR_CROSS to add a new term to the model equal to the product of dna and the last covariate; the significance of this term is then tested using t-statistic
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
feat_type_sets = feat_types_annots


## do eQTL ---------------------------------------------

start = Sys.time()
if (caucasians_only) meta_file0 = meta_file0[grepl("caucasian",meta_file0[,"race"]),]

# foreach (feat_type_set = feat_type_sets) %dopar% {
for (feat_type_set in feat_type_sets) {  
  start1 = Sys.time()
  cat("\n",feat_type_set, " ", sep=" ")
  
  # get dna & gene/isoforms data
  source(paste0(root,"/code/_func-asthma_mset0-load.R"))
  
  class_coli = 1
  for (file_ind_n in names(file_inds)) {
    # for (f1_bin in f1_bins) {
      # get row files/samples & covariate
      meta_file = meta_file0[meta_file0[,id_col] %in% Reduce(intersect,lapply(m_set0,rownames)),, drop=F]
      rownames(meta_file) = meta_file[,id_col]
      if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiment
      
      # trim matrices
      source(paste0(root,"/code/_func-asthma_mset-trim.R"))
      
      # get row files/samples & covariate
      meta_file = meta_file[match(rownames(m_set[[1]]),meta_file[,id_col]),]
      if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiment
      
      cvrt00 = meta_file[,interested_cols,drop=F]
      cvrt_num_ind = apply(cvrt00,2, function(x) all(!is.na(as.numeric(x[!is.na(x) & x!=""]))) )
      cvrt0 = as.data.frame(lapply(append(which(cvrt_num_ind),which(!cvrt_num_ind)), function(x) {
        if (cvrt_num_ind[x]) return(as.numeric(cvrt00[,x]))
        if (length(unique(cvrt00[,x]))>2) return(mtabulate(cvrt00[,x]))
        return(as.numeric(factor(cvrt00[,x]))-1)
      })); rownames(cvrt0) = rownames(meta_file)
      
      # get col feature meta
      meta_col_set = NULL
      for (x in 1:length(feat_name)) {
        file_name = paste0(meta_col_dir,feat_name[x],".Rdata")
        if (!file.exists(file_name)) next()
        meta_col_setx = get(load(file_name))
        if (!all(coln[[x]]%in%colnames(meta_col_setx))) next()
        meta_col_set[[feat_type_set[x]]] = 
          meta_col_setx[match(colnames(m_set[[feat_type_set[x]]]), meta_col_setx[,id_col]),]
      }
      cole = feat_type_set%in%names(meta_col_set)
      
      for (useModel_ in useModels) {
        # # Output file name
        # output_file_name_cis = tempfile();
        # output_file_name_tra = tempfile();
        eqtl_cis_dir = paste0(eqtl_dir,"/", feat_type_set[1], 
                              "-", feat_type_set[2], "-", file_ind_n, "Xall",
                              "_class-", paste(interested_cols,collapse="."),"_", 
                              paste0(names(table(meta_file[,class_col])),
                                     table(meta_file[,class_col]), collapse="v"),
                              "_", useModel_,"_cisdist-",cisDist, 
                              ifelse(caucasians_only, "_cauc",""), ".Rdata")
        # eqtl_cis_dir = paste0(eqtl_dir,"/", feat_type_set[1], ifelse(f1_bin!="",".",""), f1_bin,
        #                       "-", feat_type_set[2], "-", file_ind_n, "Xall",
        #                       "_class-", paste(interested_cols,collapse="."),"_", 
        #                       paste0(names(table(meta_file[,class_col])),
        #                              table(meta_file[,class_col]), collapse="v"),
        #                       "_", useModel_,"_cisdist-",cisDist, 
        #                       ifelse(caucasians_only, "_cauc",""), ".Rdata")
        eqtl_tra_dir = NULL; if (writetra) eqtl_tra_dir = gsub(".Rdata","tra_.Rdata",eqtl_cis_dir)
        
        # set model
        useModel = NULL
        if (useModel_=="modelLINEAR") useModel = modelLINEAR # dna is assumed to have only additive effect on expression
        if (useModel_=="modelANOVA") useModel = modelANOVA # assume dna to have both additive and dominant effects (ANOVA model). In this case dna data set musts take at most 3 distinct values (i.e. 0/1/2/NA); NOTE: can't do if dna has more than three categories
        if (useModel_=="modelLINEAR_CROSS") useModel = modelLINEAR_CROSS
        
        print(paste0(eqtl_cis_dir))
        if1 = file.exists(eqtl_cis_dir)
        if2 = T; if (!is.null(eqtl_tra_dir)) if2 = file.exists(eqtl_tra_dir)
        try({
          ## prepare m genotyping and rnaseq data as SlicedData for input into MatrixEQTL()
          cvrt = SlicedData$new()
          cvrt$CreateFromMatrix(t(as.data.frame(cvrt0)))
          
          # reformat col features
          f1_pos = data.frame(snpid=colnames(m_set[[1]]), chr=paste0("chr",1), pos=1)
          if (cole[1]) { x = names(m_set)[1]
          #  data.frame with columns snpid (Snp_01), chr (1), pos (725123)
          f1_pos = data.frame(snpid=as.character(meta_col_set[[x]]$dbSNP),
                              chr=paste0("chr",meta_col_set[[x]]$chromosome),
                              pos=as.numeric(meta_col_set[[x]]$pos_phys), stringsAsFactors=F)
          # if (sum(is.na(f1_pos$snpid))>0)
          f1_pos$snpid[is.na(f1_pos$snpid)] = colnames(m_set[[1]])[is.na(f1_pos$snpid)]
          f1_pos$chr[grepl("NA",f1_pos$chr)] = paste0("chr",1)
          f1_pos$pos[is.na(f1_pos$pos)] = 1
          }
          f1_pos$snpid[is.na(f1_pos$snpid)] = colnames(m_set[[1]])[is.na(f1_pos$snpid)]
          colnames(m_set[[1]]) = f1_pos$snpid
          
          f2_pos = data.frame(geneid=colnames(m_set[[2]]), chr=paste0("chr",1), left=1, right=1)
          if (cole[2]) { x = names(m_set)[2]
          #  data.frame with columns geneid (Gene_01), chr (1), left (721289), right (731289)
          f2_pos = data.frame(geneid=as.character(meta_col_set[[x]]$symbol),
                              chr=paste0("chr",meta_col_set[[x]]$chr),
                              left=as.numeric(meta_col_set[[x]]$start),
                              right=as.numeric(meta_col_set[[x]]$end), stringsAsFactors=F)
          f2_pos$geneid[is.na(f2_pos$geneid)] = colnames(m_set[[2]])[is.na(f2_pos$geneid)]
          f2_pos$chr[grepl("NA",f2_pos$chr)] = paste0("chr",1)
          f2_pos$left[is.na(f2_pos$left)] = f2_pos$right[is.na(f2_pos$right)] = 1
          }
          f2_pos$snpid[is.na(f2_pos$geneid)] = colnames(m_set[[2]])[is.na(f2_pos$geneid)]
          colnames(m_set[[2]]) = f2_pos$geneid
          
          #  Can be real-valued for linear models
          #  and must take at most 3 distinct values for ANOVA
          #  unless the number of ANOVA categories is set to a higher number (see useModel parameter).
          #  Must have matching columns
          f1_sd = SlicedData$new()
          f1_sd$CreateFromMatrix(t(m_set[[1]]))
          
          f2_sd = SlicedData$new()
          f2_sd$CreateFromMatrix(t(m_set[[2]]))
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
          
          # to overwrite or not
          
          if (!(if1 & if2 & !overwrite)) {
            
            
            
            
            
            # ## TESTING
            # i = 1
            # x = factor(f1_m[,i])
            # y = f2_m[,i]
            # lm_featresp = lm(x ~ y*meta_file$response)
            #
            #
            # library(MASS)
            # data(crabs)
            # library("reshape2")
            # melt_crabs <- melt(crabs,id.var=c("sp","sex","index"))
            #
            # ggplot(aes(x=interaction(x, y), y=y)) +
            #   geom_boxplot(aes(fill=y), alpha=.5) +
            #   geom_line(aes(group=y), alpha = 0.5, colour="darkgrey") +
            #   facet_grid(meta_file$response~.,scales="free_x") +
            #   scale_x_discrete(labels="")
            #
            # ggplot(melt_crabs, aes(x = interaction(sex, variable), y = value)) +
            #   geom_boxplot(aes(fill = sex), alpha = 0.5) +
            #   geom_line(aes(group = interaction(index, variable)),
            #             alpha = 0.5, colour = "darkgrey") +
            #   facet_grid(sp~variable,scales="free_x") +
            #   scale_x_discrete(labels = "")
            #
            #
            #
            # hsb2 <- read.csv("https://stats.idre.ucla.edu/stat/data/hsb2.csv")
            #
            # # creating the factor variable
            # lmt = lm(write ~ factor(race)*factor(female), data = hsb2)
            # ggplot(hsb2, aes(x=factor(race), y=write, fill=factor(female))) +
            #   geom_boxplot(position=position_dodge(1)) +
            #   stat_summary(fun.y=mean, geom="line")
            # # + geom_jitter(width=0.2)
            # summary(lmt)
            # # contrasts, including treatment, Helmert, sum and poly
            # hsb2 <- within(hsb2, {
            #   race.ct <- C(race.f, treatment)
            #   print(attributes(race.ct))
            # })
            # summary(lm(write ~ race.ct, data = hsb2))
            # hsb2 <- within(hsb2, {
            #   race.ch <- C(race.f, helmert)
            #   print(attributes(race.ch))
            # })
            # summary(lm(write ~ race.ch, data = hsb2))
            # hsb2 <- within(hsb2, {
            #   race.ch1 <- C(race.f, helmert, 3)
            #   print(attributes(race.ch1))
            # })
            # summary(lm(write ~ race.ch1, data = hsb2))
            #
            
            
            ## Run the analysis
            me = NULL
            filename = tempfile()
            me = Matrix_eQTL_main(
              snps = f1_sd,
              gene = f2_sd,
              cvrt = cvrt, # SlicedData object with additional covariates. Can be an empty SlicedData object in case of no covariates. The constant is always included in the model and would cause an error if included in cvrt. The order of columns must match those in snps and gene.
              
              output_file_name = filename, # significant associations (all significant associations if pvOutputThreshold=0 or only distant if pvOutputThreshold>0). If the file with this name exists, it is overwritten.
              output_file_name.cis = filename, #output_file_name_cis=tempfile(); significant local associations
              
              pvOutputThreshold = pvalthres_tra,
              pvOutputThreshold.cis = pvalthres_cis,
              
              useModel = useModel,
              # errorCovariance = errorCovariance, # numeric. The error covariance matrix. Use numeric() for homoskedastic independent errors.
              verbose = FALSE,
              
              snpspos = f1_pos,
              genepos = f2_pos,
              cisDist = ifelse(all(cole), cisDist, Inf),
              pvalue.hist = "qqplot", # logical, numerical, or "qqplot" (faster if false); To record information for a histogram set pvalue.hist to the desired number of bins of equal size. Finally, pvalue.hist can also be set to a custom set of bin edges.
              min.pv.by.genesnp = T, # record the minimum p-value for each SNP and each gene in the returned object. The minimum p-values are recorded even if if they are above the corresponding thresholds of pvOutputThreshold and pvOutputThreshold.cis (faster if false)
              noFDRsaveMemory = F # save significant gene-SNP pairs directly to the output files, reduce memory footprint and skip FDR calculation. The eQTLs are not recorded
            )
            unlink(filename)
            
            if (is.null(me)) next()
            time_output(start1, message=paste0("sig=",nrow(me$cis$eqtl)))
            if (me$all$neqtls==0) next()
            
            try ({
              me$trans$eqtl = cbind(me$trans$eqtl, rep("trans",nrow(me$trans$eqtl)))
              me$cis$eqtl = cbind(me$cis$eqtl, rep("cis",nrow(me$cis$eqtl)))
              colnames(me$trans$eqtl)[ncol(me$trans$eqtl)] = colnames(me$cis$eqtl)[ncol(me$cis$eqtl)] = "cis_trans"
            })
            
            if (all(f1_pos$chr=="chr1") | !all(f2_pos$chr=="chr1")) {
              mecist = me$trans$eqtl
            } else {
              mecist = rbind(me$cis$eqtl, me$trans$eqtl)
              # save(me$cis$eqtl, file=gsub(".Rdata","_cis.Rdata",eqtl_cis_dir))
              # save(me$trans$eqtl, file=gsub(".Rdata","_trans.Rdata",eqtl_cis_dir))
            }
            if (nrow(mecist)>0) {
              mecist = mecist[!duplicated(paste0(mecist$snps,"_",mecist$gene)),, drop=F]
              save(mecist, file=eqtl_cis_dir)
              rown = meta_file[,id_col]
              write.csv(rown,file=gsub(".Rdata","_id.csv",eqtl_cis_dir),row.names=F)
            }
            
            
            png(file=gsub(".Rdata","_qq.png",eqtl_cis_dir))
            plot(me)
            graphics.off()
            
            ## Results -----------------------------------------
            cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
            # cat('Detected local eQTLs:', '\n');
            # # write.csv(me$trans$eqtls, file=gsub(".Rdata",".csv",gsub("cis","tra",eqtl_cis_dir)))
            # cat('Detected distant eQTLs:', '\n');
            # show(me$trans$eqtls)
          }
        })
      } # useModel
    # } # f1_bin
  } #file_inds
} # feat_type_set


time_output(start)


