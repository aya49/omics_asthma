## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720


## logistics
root = "~/projects/asthma"; commandArgs <- function(...) root  # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/visualizationFunctions.R"))
libr(append(pkgs(),c("pROC", "rafalib", "ROCR", "CellCODE", "GGally", "mixOmics")))

no_cores = 5#detectCores()-3
registerDoMC(no_cores)


## options
overwrite = F

pthres = .01

# plot
height = 400
width = 500

ncomp = 2

# if there is a cell feat, there must be a rna feat with gene as colnames
feat_types1 = list(list(rna="rnaseqstarens.pre",
                       rna.pc="rnapcgenes.pre",
                       rna.elements="rnaelements.pre",
                       gene="genotype"),
                  list(gene="genotype"),
                  list(rna.pc="rnapcgenes.pre"),
                  list(rna.elements="rnaelements.pre"),
                  list(rna="rnaseqstarens.pre"),
                  list(rna="rnaseqstarens.pre",
                       rna.pc="rnapcgenes.pre",
                       rna.elements="rnaelements.pre",
                       cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                       metab="metab.pre",
                       gene="genotype"),
                  list(rna="rnaseqstarens.post",
                       cell="cellseqgenes.post",
                       metab="metab.post",
                       gene="genotype"),
                  list(rna="rnaseqstarens.diff",
                       cell="cellseqgenes.diff",
                       metab="metab.diff",
                       gene="genotype"),
                  list(rna="rnaseqstarens.pre",
                       rna.pc="rnapcgenes.pre",
                       rna.elements="rnaelements.pre",
                       cell="cell.pre",
                       metab="metab.pre",
                       gene="genotype"),
                  list(rna="rnaseqstarens.post",
                       cell="cell.post",
                       metab="metab.post",
                       gene="genotype"),
                  list(rna="rnaseqstarens.diff",
                       cell="cell.diff",
                       metab="metab.diff",
                       gene="genotype"))

feat_types = append(feat_types_annots,feat_types1)

celldiff = get(load(meta_col_rnacells_dir))

for (feat_type in feat_types) {
  m00 = lapply(feat_type, function(x) get(load(paste0(feat_dir,"/",x,".Rdata"))))
  m00_names = Reduce("intersect",lapply(m00,rownames))
  
  for (file_ind_n in names(file_inds)) {
    # prepare meta_file indices
    file_ind = file_inds[[file_ind_n]]
    if (file_ind_n!="all" & all(rownames(m00_names)%in%file_ind)) next()
    if (file_ind_n=="all") file_ind = m00_names
    
    m0_names = intersect(m00_names,file_ind)
    m0 = lapply(m00, function(x) x[m0_names,] )
    
    meta_col00 = lapply(feat_type, function(x) {
      mcname = paste0(meta_col_dir,"-", str_split(x,"[.]")[[1]][1],".Rdata")
      if (file.exists(mcname)) return(get(load(mcname)))
      return(NA)
    })
    names(m0) = names(meta_col00) = names(feat_type)
    
    for (bin in c("12","01","")) {
      
      
      # prepare meta data for samples (their classes)
      meta_file = meta_file0[match(m0_names,meta_file0[,id_col]),]
      Y = factor(meta_file[,class_col], levels = c("ER", "DR"))
      names(Y) = meta_file[,id_col]
      
      # prepare actual data
      X0 = list()
      meta_col0 = list()
      for (x in names(m0)) {
        mx = m0[[x]]
        if (x=="gene") { #trim based on p value
          gwf = list.files(gwas_dir, pattern=".Rdata", full.names=T)
          gwfi = gwf[grepl("goodpplXall",gwf) & grepl("genotype",gwf)]
          gwas_g = get(load(gwfi))
          gwas_gl = gwas_g[,grepl("none",colnames(gwas_g))]
          mx = mx[,names(gwas_gl)[gwas_gl<pthres]]
          
          if (bin=="01") mx[mx==2] = 1
          if (bin=="12") mx[mx==0] = 1
        }
        mx = mx[apply(mx,1,function(y) any(!is.na(y))), 
                apply(mx,2,function(y) sum(!is.na(y))>(good_na*length(y)))]
        if (x=="gene") mx = mx[,apply(mx,2,function(x) min(table(x))>good_col & length(table(x))>1)]
        if (!is.null(dim(meta_col00[[x]]))) meta_col0[[x]] = meta_col00[[x]][match(colnames(mx),meta_col00[[x]][,id_col]),]
        if (x=="gene") colnames(mx) = meta_col0[[x]][,"dbSNP"]
        if (grepl("rna",x)) colnames(mx) = meta_col0[[x]][,"symbol"]
        if (grepl("rna",x) | x=="gene") {
          dupind = duplicated(colnames(mx), fromLast=F) | duplicated(colnames(mx), fromLast=T)
          colnames(mx)[dupind] = paste0(colnames(mx)[dupind], "_", meta_col0[[x]][dupind,id_col])
        }
        X0[[x]] = mx
      }

      X = X0
      meta_col = meta_col0
      
      
      #get rid of rna genes that went into inferring cell count if cell count is used
      if (sum(grepl("cellseqgenes",names(X)))>0 & sum(grepl("rna",names(X)))>0)
        X[grepl("rna",names(X))] = lapply(X[grepl("rna",names(X))], function(X) x[,!colnames(x)%in%celldiff])
      
      #bind rna type data sets together
      if (sum(grepl("rna",names(meta_col0)))>1) {
        Xrna = Reduce("cbind",X[grepl("rna",names(X))])
        X = X[!grepl("rna",names(X))]
        X$rna = Xrna
        
        meta_colrna = Reduce("rbind",meta_colrna[grepl("rna",names(meta_col))])
        meta_col = meta_colrna[!grepl("rna",names(meta_col))]
        
      }
      
      # prepare number of pls-da factors you want for each data
      design = matrix(0, nrow = length(X), ncol = length(X))
      
      
      
      
      for (tune_ in c("","tune")) {
        
        # prepare file name
        pname = paste0(blocksplsda_dir, "/", paste(feat_type,collapse="-"), ifelse(bin=="","","."), bin, "-",  file_ind_n, "Xall", "_pthres-", pthres,ifelse(tune_=="","","_"),tune_)
        if (overwrite & file.exists(paste0(pname,".Rdata"))) next()
        
        keepX = lapply(names(X), function(x) {
          if (x=="cell") return(rep(2, ncomp))
          if (grepl("rna",x) | x=="gene") return(rep(10, ncomp))
          if (x=="metab") return(rep(5, ncomp))
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
              hk.known = meta_col0[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col0[[feat]][,"symbol"]),]
            } else if (feat=="gene") {
              hk.known = meta_col0[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col0[[feat]][,"dbSNP"]),]
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
    } #bin
  } #file_ind
} #feat_type
