## input: features & meta_file
## output: eQTL
## aya43@sfu.ca
## created 20180614


# Cell Specific eQTL Analysis without Sorting Cells
# y ~ g*c # y=rnaseq, g=genotype, c=cells
# NOTE: Yi = β0 + Xiβ1 + Ziβ2 + Wiβ3 + XiZiβ4 + XiWiβ5 + ZiWiβ6 + XiZiWiβ7 + ei
# Y ~ X + Z + W + X:Z + X:W + Z:W + X:Z:W
# Y ~ X * Z * W
# Y ~ (X + Z + W)^3


## logistics
root = "~/projects/asthma"; commandArgs <- function(...) root  # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/_func-classifiers.R"))
libr(append(pkgs(),c("qdaptools","MatrixEQTL")))

no_cores = detectCores()-2 #number of cores to use for parallel processing
registerDoMC(no_cores)


## options
overwrite = T # overwrite if existing analysis eexists?
writecsv = T # write results as csv on top of Rdata?

writetra = F # write all associations not only local ones

categorical = T # is class column categorical?
# interested_cols = c("response")
interested_cols = c("sex","response")
# interested_cont_cols = ""

caucasians_only = T

f1_bins = c("","01","12") # 01 make all 2s 1, 12 make all 0s 1
# split_f1_col = "time" #there are
useModels = c("modelLINEAR", "modelANOVA", "modelLINEAR_CROSS")

# plotting size
width = 800
height = 600

# eqtl cut-offs
pvalthres_cis = .005
pvalthres_tra = 1e-6
# useModel = feat_type[5]
# errorCovariance = numeric()
cisDist = 1e6


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
feat_types = feat_types_annots


## do eQTL ---------------------------------------------

start = Sys.time()


foreach (feat_type = feat_types) %dopar% {
  start1 = Sys.time()
  cat(feat_type, " ", sep=" ")

  # name parameters
  feat_type1 = feat_type[1]
  feat_type2 = feat_type[2]

  # get genotype & gene/isoforms data
  f1_m0 = get(load(paste0(feat_dir,"/",feat_type1,".Rdata")))
  f2_m0 = get(load(paste0(feat_dir,"/",feat_type2,".Rdata")))
  f1_m0 = f1_m0[rownames(f1_m0)%in%meta_file0[,id_col],]
  f2_m0 = f2_m0[rownames(f2_m0)%in%meta_file0[,id_col],]

  for (file_ind_n in names(file_inds)) {
    for (f1_bin in f1_bins) {

      if (feat_type1!="genotype" & f1_bin!="") next()

      file_ind = file_inds[[file_ind_n]]
      if (file_ind_n=="all") file_ind = meta_file0[,id_col]

      # get row files/samples & covariate
      samples_to_include = intersect(rownames(f2_m0), rownames(f1_m0))
      meta_file =
        meta_file0[
          # !is.na(meta_file0[,class_col]) &
            meta_file0[,id_col] %in% samples_to_include &
            ifelse(caucasians_only, grepl("caucasian",meta_file0[,"race"]), rep(T,nrow(meta_file0))) &
            meta_file0[,id_col] %in% file_ind
        ,,drop=F]
      if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiments
      rownames(meta_file) = meta_file[,id_col]

      cvrt00 = meta_file[,interested_cols,drop=F]
      cvrt_num_ind = apply(cvrt00,2, function(x)
        all(!is.na(as.numeric(x[!is.na(x) & x!=""]))) )
      cvrt0 = as.data.frame(lapply(append(which(cvrt_num_ind),which(!cvrt_num_ind)), function(x) {
        if (cvrt_num_ind[x]) return(as.numeric(cvrt00[,x]))
        if (length(unique(cvrt00[,x]))>2) return(mtabulate(cvrt00[,x]))
        return(as.numeric(factor(cvrt00[,x]))-1)
      }))
      rownames(cvrt0) = rownames(meta_file)

      f1_m = f1_m0[as.character(meta_file[,id_col]),, drop=F]
      f2_m = f2_m0[as.character(meta_file[,id_col]),, drop=F]

      # get col features
      f1_meta_col_ind = f1_meta_col_ind0 = apply(f1_m0,2,function(x) sum(!is.na(x))>(good_na*length(x)))
      f2_meta_col_ind = f2_meta_col_ind0 = apply(f2_m0,2,function(x) sum(!is.na(x))>(good_na*length(x)))
      if (sum(f1_meta_col_ind)==0 | sum(f2_meta_col_ind)==0) next()

      if (f1_bin=="01") f1_m[f1_m==2] = 1
      if (f1_bin=="12") f1_m[f1_m==0] = 1

      if (feat_type1=="genotype") f1_meta_col_ind = f1_meta_col_ind0 &
        apply(f1_m,2,function(x) { tx = table(x); return(min(tx)>good_col & length(tx)>1) })

      f1_cole = file.exists(paste0(meta_col_dir,"-",feat_type1,".Rdata"))
      if (f1_cole) {
        f1_meta_col0 = get(load(paste0(meta_col_dir,"-",feat_type1,".Rdata")))
        f1_cole = f1_cole & all(c("dbSNP","chromosome","pos_phys")%in%colnames(f1_meta_col0))
        if (f1_cole) {
          f1_meta_col0 = f1_meta_col0[match(colnames(f1_m0),f1_meta_col0[,id_col]),]
          f1_meta_col_ind = f1_meta_col_ind & apply(f1_meta_col0[,c("dbSNP","chromosome","pos_phys")],1,function(x) !any(is.na(x)))
          f1_meta_col = f1_meta_col0[f1_meta_col_ind,]
        }
      }

      f2_cole = file.exists(paste0(meta_col_dir,"-",feat_type2,".Rdata"))
      if (f2_cole) {
        f2_meta_col0 = get(load(paste0(meta_col_dir,"-",feat_type2,".Rdata")))
        f2_cole = f2_cole & all(c("symbol","chr","start","end")%in%colnames(f2_meta_col0))
        if (f2_cole & sum(is.na(f2_meta_col0))>(.5*nrow(f2_meta_col0))) {
          f2_cole = F
          colnames(f1_m) = paste0(f2_meta_col0$symbol[match(colnames(f1_m),f2_meta_col0[,id_col])],"_",colnames(f1_m))
        }
        if (f2_cole) {
          f2_meta_col0 = f2_meta_col0[match(colnames(f2_m0),f2_meta_col0[,id_col]),]
          f2_meta_col_ind = f2_meta_col_ind & apply(f2_meta_col0[,c("symbol","chr","start","end")],1,function(x) !any(is.na(x)))
          f2_meta_col = f2_meta_col0[f2_meta_col_ind,]
          # if (colnames(f1_m0)[1]%in%f1_meta_file0[,id_col]) f1_m0 = t(f1_m0)
        }
      }

      # trim matrices
      f1_m = f1_m[, f1_meta_col_ind, drop=F]
      f2_m = f2_m[, f2_meta_col_ind, drop=F]
      if (sum(f1_meta_col_ind)==0 | sum(f2_meta_col_ind)==0) next()
      s
      for (useModel_ in useModels) {
        # # Output file name
        # output_file_name_cis = tempfile();
        # output_file_name_tra = tempfile();
        eqtl_cis_dir = paste0(eqtl_dir,"/", feat_type1,
                              ifelse(f1_bin!="",".",""),f1_bin,"-", feat_type2,
                              "-", file_ind_n, "Xall",
                              "_class-", paste(interested_cols,collapse="."),
                              "_", useModel_,"_cisdist-",cisDist, ".Rdata")
        eqtl_tra_dir = NULL; if (writetra) eqtl_tra_dir = gsub(".Rdata","tra_.Rdata",eqtl_cis_dir)

        # set model
        useModel = NULL
        if (useModel_=="modelLINEAR") useModel = modelLINEAR # genotype is assumed to have only additive effect on expression
        if (useModel_=="modelANOVA") useModel = modelANOVA # assume genotype to have both additive and dominant effects (ANOVA model). In this case genotype data set musts take at most 3 distinct values (i.e. 0/1/2/NA)
        if (useModel_=="modelLINEAR_CROSS") useModel = modelLINEAR_CROSS


        print(paste0(eqtl_cis_dir))
        if1 = file.exists(eqtl_cis_dir)
        if2 = T; if (!is.null(eqtl_tra_dir)) if2 = file.exists(eqtl_tra_dir)
        try({
          ## prepare m genotyping and rnaseq data as SlicedData for input into MatrixEQTL()

          cvrt = SlicedData$new()
          cvrt$CreateFromMatrix(t(as.numeric(as.matrix(cvrt0))))

          # reformat col features
          f1_pos = data.frame(snpid=colnames(f1_m), chr=paste0("chr",1), pos=1, stringsAsFactors=F)
          if (f1_cole) {
            #  data.frame with columns snpid (Snp_01), chr (1), pos (725123)
            f1_pos = data.frame(snpid=as.character(f1_meta_col$dbSNP),
                                chr=paste0("chr",f1_meta_col$chromosome),
                                pos=as.numeric(f1_meta_col$pos_phys),
                                stringsAsFactors=F)
          }

          f2_pos = data.frame(geneid=colnames(f2_m), chr=paste0("chr",1), left=1, right=1, stringsAsFactors=F)
          if (f2_cole) {
            #  data.frame with columns geneid (Gene_01), chr (1), left (721289), right (731289)
            f2_pos = data.frame(geneid=as.character(f2_meta_col$symbol),
                                chr=paste0("chr",f2_meta_col$chr),
                                left=as.numeric(f2_meta_col$start),
                                right=as.numeric(f2_meta_col$end),
                                stringsAsFactors=F)
          }

          colnames(f1_m) = f1_pos$snpid
          colnames(f2_m) = f2_pos$geneid

          #  Can be real-valued for linear models
          #  and must take at most 3 distinct values for ANOVA
          #  unless the number of ANOVA categories is set to a higher number (see useModel parameter).
          #  Must have matching columns
          f1_sd = SlicedData$new()
          f1_sd$CreateFromMatrix(t(f1_m))

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
              cisDist = ifelse(f1_cole & f2_cole, cisDist, Inf),
              pvalue.hist = "qqplot", # logical, numerical, or "qqplot" (faster if false); To record information for a histogram set pvalue.hist to the desired number of bins of equal size. Finally, pvalue.hist can also be set to a custom set of bin edges.
              min.pv.by.genesnp = T, # record the minimum p-value for each SNP and each gene in the returned object. The minimum p-values are recorded even if if they are above the corresponding thresholds of pvOutputThreshold and pvOutputThreshold.cis (faster if false)
              noFDRsaveMemory = F # save significant gene-SNP pairs directly to the output files, reduce memory footprint and skip FDR calculation. The eQTLs are not recorded
            )
            unlink(filename)

            if (is.null(me)) next()
            time_output(start1, message=paste0("sig=",nrow(me$cis$eqtl)))
            if (me$all$neqtls==0) next()

            if (all(f1_pos$chr=="chr1") | !all(f2_pos$chr=="chr1")) {
              save(me$trans$eqtl, file=eqtl_cis_dir)
            } else {
              mecist = rbind(cbind(rep("cis",nrow(me$cis$eqtl)), me$cis$eqtl),
                             cbind(rep("trans",nrow(me$trans$eqtl)), me$trans$eqtl))
              # save(me$cis$eqtl, file=gsub(".Rdata","_cis.Rdata",eqtl_cis_dir))
              # save(me$trans$eqtl, file=gsub(".Rdata","_trans.Rdata",eqtl_cis_dir))
              if (nrow(mecist)>0) {
                mecist = mecist[!duplicated(paste0(mecist$snps,"_",mecist$gene)),, drop=F]
                colnames(mecist)[1] = "cis_trans"
                save(mecist, file=eqtl_cis_dir)
              }
            }

            png(file=gsub(".Rdata","_qq.png",eqtl_cis_dir))
            plot(me)
            graphics.off()

            ## Results -----------------------------------------
            cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
            cat('Detected local eQTLs:', '\n');
            # write.csv(me$trans$eqtls, file=gsub(".Rdata",".csv",gsub("cis","tra",eqtl_cis_dir)))
            cat('Detected distant eQTLs:', '\n');
            # show(me$trans$eqtls)
          }
        })
      } # useModel
    } # f1_bin
  } #file_inds
} # feat_type


time_output(start)

