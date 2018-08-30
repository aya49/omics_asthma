## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/visualizationFunctions.R"))
libr(append(pkgs(),c("pROC", "rafalib", "ROCR", "CellCODE", "GGally", "mixOmics")))

# no_cores = 5#detectCores()-3
# registerDoMC(no_cores)


## options
overwrite = F

pthres = .01

# plot
height = 400
width = 500

ncomp = 2

scale_cont = F #if values are continuous, scale?
caucasians_only = T



dnasigonly = T # if dna; only keep significant snps


# if there is a cell feat, there must be a rna feat with gene as colnames
# feat with same name will be binded together e.g. rna="rnapcgenes.pre" and rna="rnaelements.pre"
feat_type_sets1 = list(unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna="rnapcgenes.pre",
                                   rna="rnaelements.pre",
                                   dna="dna")),
                       unlist(list(dna="dna")),
                       unlist(list(rna.pc="rnapcgenes.pre")),
                       unlist(list(rna.elements="rnaelements.pre")),
                       unlist(list(rna="rnaseqstarens.pre")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                                   metab="metab.pre",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.post",
                                   cell="cellseqgenes.post",
                                   metab="metab.post",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.diff",
                                   cell="cellseqgenes.diff",
                                   metab="metab.diff",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   cell="cell.pre",
                                   metab="metab.pre",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.post",
                                   cell="cell.post",
                                   metab="metab.post",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.diff",
                                   cell="cell.diff",
                                   metab="metab.diff",
                                   dna="dna")))

feat_type_sets = append(feat_types_annots,feat_type_sets1)

celldiff = get(load(meta_col_rnacells_dir))
if (caucasians_only) meta_file0 = meta_file0[grepl("caucasian",meta_file0[,"race"]),]

for (feat_type_set in feat_type_sets) {
  start1 = Sys.time()
  cat("\n", unlist(feat_type_set), " ", sep=" ")
  
  # get dna & gene/isoforms data
  source(paste0(root,"/code/_func-asthma_mset0-load.R"))
  
  class_coli = 1
  for (file_ind_n in names(file_inds)) {
    # for (f1_bin in f1_bins) {
      # get row files/samples & covariate
      # get row files/samples & covariate
      # get row files/samples & covariate
      meta_file = meta_file0[meta_file0[,id_col] %in% Reduce(intersect,lapply(m_set0,rownames)),, drop=F]
      rownames(meta_file) = meta_file[,id_col]
      if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiment
      
      # trim matrices
      source(paste0(root,"/code/_func-asthma_mset-trim.R"))
      
      # get row files/samples & covariate
      meta_file = meta_file[match(rownames(m_set[[1]]),meta_file[,id_col]),]
      if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiment
      
      # prepare meta data for samples (their classes)
      Y = factor(meta_file[,class_col], levels = c("ER", "DR"))
      names(Y) = meta_file[,id_col]
      
      # get col feature meta
      meta_col_set = NULL
      for (x in 1:length(feat_name)) {
        file_name = paste0(meta_col_dir,feat_name[x],".Rdata")
        if (!file.exists(file_name)) next()
        meta_col_setx = get(load(file_name))
        meta_col_set[[feat_type_set[x]]] = 
          meta_col_setx[match(colnames(m_set[[feat_type_set[x]]]), meta_col_setx[,id_col]),]
      }
      cole = feat_type_set%in%names(meta_col_set)
      
      # prepare actual data
      X0 = list()
      for (x in names(m_set)) {
        mx = m_set[[x]]
        if (grepl("dna",x) & dnasigonly) { #trim based on p value
          gwf = list.files(gwas_dir, pattern=".Rdata", full.names=T)
          gwas_g = get(load(gwf[grepl(paste0(x,"-",file_ind_n,"Xall"),gwf)]))
          gwas_gl = gwas_g[,grepl("none",colnames(gwas_g))]
          mx = mx[,names(gwas_gl)[gwas_gl<pthres]]
          if (!is.null(dim(meta_col_set[[x]]))) meta_col_set[[x]] = meta_col_set[[x]][match(colnames(mx),meta_col_set[[x]][,id_col]),]
        }
        if (grepl("dna",x)) colnames(mx)[!is.na(meta_col_set[[x]][,"dbSNP"])] = meta_col_set[[x]][!is.na(meta_col_set[[x]][,"dbSNP"]),"dbSNP"]
        if (grepl("rna",x)) colnames(mx)[!is.na(meta_col_set[[x]][,"symbol"])] = meta_col_set[[x]][!is.na(meta_col_set[[x]][,"symbol"]),"symbol"]
        if (grepl("rna",x) | grepl("dna",x)) {
          dupind = duplicated(colnames(mx), fromLast=F) | duplicated(colnames(mx), fromLast=T)
          colnames(mx)[dupind] = paste0(colnames(mx)[dupind], "_", meta_col_set[[x]][dupind,id_col])
        }
        X0[[x]] = mx
      }

      # 0 has all the original feat names
      X = X0
      meta_col_set0 = meta_col_set 
      
      
      #get rid of rna genes that went into inferring cell count if cell count is used
      if (sum(grepl("cellseqgenes",names(X)))>0 & sum(grepl("rna",names(X)))>0)
        X[grepl("rna",names(X))] = lapply(X[grepl("rna",names(X))], function(X) x[,!colnames(x)%in%celldiff])
      
      #bind feat with same name together
      dupfeat = feat_type_set[feat_type_set%in%names(X)]
      dupfeat = dupfeat[duplicated(names(dupfeat)) | duplicated(names(dupfeat), fromLast=T)]
      if (length(dupfeat)>0) {
        for (duf in unique(names(dupfeat))) {
          dupfeati = dupfeat[names(dupfeat)%in%duf]
          Xrna = Reduce("cbind",X[dupfeati])
          X = X[!names(X)%in%dupfeati]
          X[[duf]] = Xrna
          
          meta_col_setrna = meta_col_set[names(meta_col_set)%in%dupfeati]
          rnacols = Reduce(intersect, lapply(meta_col_setrna, colnames))
          meta_col_setrna = Reduce(rbind, lapply(meta_col_setrna, function(x) x[,rnacols]))
          meta_col_set = meta_col_set[!names(meta_col_set)%in%dupfeati]
          meta_col_set[[duf]] = meta_col_setrna
        }
      }
      
      # prepare number of pls-da factors you want for each data
      design = matrix(0, nrow = length(X), ncol = length(X))
      
      
      
      
      for (tune_ in c("","tune")) {
        
        # prepare file name
        pname = paste0(blocksplsda_dir, "/", paste(feat_type_set,collapse="-"), "-",  file_ind_n, "Xall", "_class-", paste(interested_cols,collapse="."),"_", paste0(names(table(meta_file[,class_col])),table(meta_file[,class_col]), collapse="v"), "_pthres-", pthres,ifelse(tune_=="","","_"),tune_)
        # pname = paste0(blocksplsda_dir, "/", paste(feat_type_set,collapse="-"), ifelse(f1_bin=="","","."), f1_bin, "-",  file_ind_n, "Xall", "_class-", paste(interested_cols,collapse="."),"_", paste0(names(table(meta_file[,class_col])),table(meta_file[,class_col]), collapse="v"), "_pthres-", pthres,ifelse(tune_=="","","_"),tune_)
        
        if (overwrite & file.exists(paste0(pname,".Rdata"))) next()
        
        keepX = lapply(names(X), function(x) {
          if (x=="cell") return(rep(2, ncomp))
          if (grepl("rna",x) | grepl("dna",x)) return(rep(10, ncomp))
          if (x=="metab") return(rep(2, ncomp))
        })
        names(keepX) = names(X)
        
        if (tune_=="tune") {
          # definition of the keepX value to be tested for each block mRNA miRNA and protein
          # names of test.keepX must match the names of 'data'
          test.keepX = lapply(names(X), function(x) seq(5,20,5))
          names(test.keepX) = names(X)
          
          # the following may take some time to run, note that for through tuning
          # nrepeat should be > 1
          if (length(X)==1) {
            tune = tune.splsda(X = X[[1]], Y = Y, ncomp = ncomp, 
                               test.keepX = test.keepX[[1]], nrepeat = 1,
                               cpus = no_cores, validation = "Mfold")
            
          } else {
            tune = tune.block.splsda(X = X, Y = Y, ncomp = ncomp, 
                                     test.keepX = test.keepX, design = design, nrepeat = 1,
                                     cpus = no_cores, validation = "Mfold")
          }
          
          tune$choice.keepX.constraint # NULL as constraint = F per default
          tune$choice.keepX
          
          pdf(paste0(pname,"_tune.pdf"))
          par(mfrow=c(2,1))
          plot(tune$error.rate[, "comp1"][order(tune$error.rate[, ncol(tune$error.rate)])] ~ factor(rownames(tune$error.rate)), col = 1, ylim = c(0.2, 1), main=as.character(which.min(tune$error.rate[,ncol(tune$error.rate)])))
          points(tune$error.rate[,ncol(tune$error.rate)][order(tune$error.rate[,ncol(tune$error.rate)])] ~ factor(rownames(tune$error.rate)), col = 2)
          dev.off()
          
          # use tuned keepX
          keepX = tune$choice.keepX
        }
        
        ## block.splsda = horizontal integration PLS-DA model with a specified number of components per block either by Y or by its position indY in the list of blocks X
        # result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
        #                       keepX = keepX, design = design,
        #                       mode = "regression", bias = T)
        if (length(X)==1) {
          if (typeof(keepX)=="list") keepX = keepX[[1]]
          result = splsda(X = X[[1]], Y = Y, ncomp = ncomp, 
                          keepX = keepX, #design = design,
                          mode = "regression")
        } else {
          result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                                keepX = keepX, design = design,
                                mode = "regression") # PLS regression ("regression"), PLS canonical analysis ("canonical"), redundancy analysis ("invariant") and the classical PLS algorithm ("classic")
        }
        
        # get variables that have >0 weights
        save(result, file=paste0(pname,".Rdata"))
        feat1 = lapply(result$loadings, function(x) apply(x, 2, function(i) names(i)[which(i != 0)]))
        
        pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"plots.pdf"), width = 7)
        try ({
          # plotDiablo(result, ncomp = 1, groupOrder = c("DR", "ER"))
          # plotDiablo3(result, ncomp = 1, groupOrder = c("DR", "ER"))
          plotIndiv_diablo(result, ncomp = 1, groupOrder = c("DR", "ER"))
        })
        dev.off()
        
        pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"arrows.pdf"), width = 7)
        try ({
          plotArrow(result)
        })
        dev.off()
        
        pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"circos.pdf"), width = 7)
        try ({
          circosPlot(result, showIntraLinks = F, cutoff=.5)#, corThreshold = 0.5, cex.label = 0.5)
          # circosPlot_diabloModif(result, corThreshold = 0.5, cex.label = 0.5, showIntraLinks = F) #grep for rna, name in X list must have rna for transcriptome!
        })
        dev.off()
        
        pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"heatmap.pdf"), width = 7)
        try ({
          heatmap_diablo(result, margins = c(2, 12))
          # cim(result, margins = c(2, 12))
          # cimDiablo(result, margins = c(2, 12))
        })
        dev.off()
        
        pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"var.pdf"), width = 7)
        try ({
          plotVar(result)
        })
        dev.off()
        
        pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"loadings.pdf"), width = 7)
        try ({
          plotLoadings(result)
        })
        dev.off()
        
        pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"network.pdf"), width = 7)
        try ({
          network(result) # network: similarity between features is obtained by calculating the sum of the correlations between the original variables and each of the latent components of the model
        })
        dev.off()

        
        # performance of splsda
        for (test in c("loo","Mfold")) {
          cv = perf(result, validation = test)
          cv$WeightedPredict.error.rate
          # cv$error.rate
          
          png(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"predictscore-",test,".png"), height=(length(X))*height, width=2*width)
          par(mfrow=c(length(X)+1,2))
          for (feat in c(names(X),"avg")) {
            
            # plot AUC
            if (length(X)>1) {
              if (feat=="avg") {
                predictScores = Reduce("+", lapply(cv$predict$nrep1, function(i) i[[ncomp]]))/length(X)
              } else {
                predictScores = cv$predict$nrep1[[feat]][[ncomp]] #average prediction
              }
            } else {
              if (feat=="avg") next()
              predictScores = cv$predict[[ncomp]][,,1]
            }
            
            roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[,"DR"], plot=F, percent = F, na.rm =T, direction = "<")   ## had the wrong direction!!
            roc.res4 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
            roc.res4$Specificity = 100 - as.numeric(roc.res4$Specificity)
            roc.res4$Sensitivity = as.numeric(roc.res4$Sensitivity)
            # roc.score$auc
            # plot(roc.score, main=paste0(feat," auc=",round(roc.score$auc,4)))
            
            plot(roc.res4$Sensitivity ~ roc.res4$Specificity, type = "o", pch = 19,
                 xlab = "100-Specificity", ylab = "Sensitivity", col = "#1F78B4",
                 main=paste0(feat," auc=",round(roc.score$auc,4)))
            abline(a = 0, b = 1, col = "gray", lty = 2)
            text(x = 80, y = 20, labels = paste0("AUC = ", roc.score$auc), col = "#1F78B4")
            
            if (feat=="avg") next()
            
            
            # save features
            if (grepl("rna",feat)) {
              hk.known = meta_col_set[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col_set[[feat]][,"symbol"]),]
            } else if (feat=="dna") {
              hk.known = meta_col_set[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col_set[[feat]][,"dbSNP"]),]
            } else {
              hk.known = as.character(feat1[[feat]])
            }
            write.csv(hk.known, paste0(pname, "_",tune_,ifelse(tune_=="","","_"), feat, "_loading.csv"))
            
            
            
            # plot(X[[1]][, "dmap_HSC3"]~ X[[3]][, "xLeucine"])
            # abline(lm(X[[1]][, "dmap_HSC3"]~ X[[3]][, "xLeucine"]))
            
            
            ## performance on each dataset separately; do plsda again #############
            # Cells
            result.cc = splsda(X = X[[feat]], Y = Y, keepX = c(2, 2))
            cv.cc = perf(result.cc, validation = "loo")
            predictScores = cv.cc$predict[[length(cv.cc$predict)]]
            roc.score = roc(response=as.character(Y), predictor=predictScores[,"DR",1], plot=F, percent=T, na.rm=T)
            roc.res3 = data.frame(Specificity=rev(roc.score$specificities), Sensitivity=rev(roc.score$sensitivities))
            roc.res3$Specificity = 100 - as.numeric(roc.res3$Specificity)
            roc.res3$Sensitivity = as.numeric(roc.res3$Sensitivity)
            # plot(roc.score, main=paste0(feat," (plsda) auc=",round(roc.score$auc,4)))
            plot(roc.res3$Sensitivity ~ roc.res3$Specificity, type = "o", pch = 19,
                 xlab = "100-Specificity", ylab = "Sensitivity", col = "#1F78B4",
                 main=paste0(feat," (splsda) auc=",round(roc.score$auc,4)))
            abline(a = 0, b = 1, col = "gray", lty = 2)
            text(x = 80, y = 20, labels = paste0("AUC = ", roc.score$auc/100), col = "#1F78B4")
            
            
            # save perf
            biomarkerGenescv = 
              unlist(lapply(strsplit(names(unlist(
                ifelse (length(X)==1, cv$features$stable, cv$features$stable$nrep1[[feat]])
              )), "[.]"), function(i) i[2]))
            write.csv(biomarkerGenescv, paste0(pname, "_",tune_,ifelse(tune_=="","","_"),feat,"_perf.csv"))
            
          } #feat
          dev.off()
          
        } #test
      } #tune_
    # } #f1_bin
  } #file_ind
} #feat_type_set
