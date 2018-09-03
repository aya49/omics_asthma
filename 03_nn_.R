## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result")


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_dna_dir = paste0(feat_dir,"/snp-file-dna")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
eqtl_dir = paste0(stat_dir,"/eqtl"); dir.create(eqtl_dir, showWarnings=F)
gwas_dir = paste0(stat_dir,"/gwas"); dir.create(gwas_dir,showWarnings=F)

super_dir = paste0(result_dir,"/supervised"); dir.create(super_dir, showWarnings=F)
blockplsda_dir = paste0(super_dir,"/blockplsda"); dir.create(blockplsda_dir,showWarnings=F)



source("src/_func.R")
source("src/visualizationFunctions.R")

libr(c("biomaRt", "mixOmics",
       "limma", "GGally", "ggplot2",
       "reshape2", "ROCR", "CellCODE", # devtools::install_github("mchikina/CellCODE")
       "pROC", "rafalib", "dplyr","stringr"))




## options
no_cores = 5#detectCores()-3
registerDoMC(no_cores)

id_col = "id"
class_col = "response"

pthres = .01

overwrite = F

height = 400
width = 500

good_na = .75 #proportion of na more than this, then delete the column in matrix
good_col = 3

ncomp = 2



# if there is a cell feat, there must be a rna feat with gene as colnames
feat_types = list(list(rna="rnaseqstarens.pre",
                       gene="dna"),
                  list(gene="dna"),
                  list(rna="rnaseqstarens.pre"),
                  list(rna="rnaseqstarens.pre",
                       cell="cellseqgenes.pre",
                       metab="metab.pre",
                       gene="dna"),
                  list(rna="rnaseqstarens.post",
                       cell="cellseqgenes.post",
                       metab="metab.post",
                       gene="dna"),
                  list(rna="rnaseqstarens.diff",
                       cell="cellseqgenes.diff",
                       metab="metab.diff",
                       gene="dna"))


meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))

for (feat_type in feat_types) {
  m0 = lapply(feat_type, function(x) get(load(paste0(feat_dir,"/",x,".Rdata"))))
  m0_names = Reduce("intersect",lapply(m0,rownames))
  m0 = lapply(m0, function(x) x[m0_names,] )
  
  meta_col0 = lapply(feat_type, function(x) {
    mcname = paste0(meta_col_dir,"-", str_split(x,"[.]")[[1]][1],".Rdata")
    if (file.exists(mcname)) return(get(load(mcname)))
    return(NA)
  })
  names(m0) = names(meta_col0) = names(feat_type)
  
  for (bin in c("12","01","")) {
    # prepare file name
    pname = paste0(blockplsda_dir, "/", paste(feat_type,collapse="-"), ifelse(bin=="","","."),bin, "_pthres-", pthres)
    if (overwrite & file.exists(paste0(pname,".Rdata"))) next()
    
    # prepare meta data for samples (their classes)
    meta_file = meta_file0[match(m0_names,meta_file0[,id_col]),]
    Y = factor(meta_file[,class_col], levels = c("ER", "DR"))
    names(Y) = meta_file[,id_col]
    
    # prepare actual data
    X = m = lapply(names(m0), function(x) {
      mx = m0[[x]]
      if (x=="gene") { #trim based on p value
        gwf = list.files(gwas_dir, pattern=".Rdata", full.names=T)
        gwfi = gwf[grepl("goodpplXall",gwf) & grepl("dna",gwf)]
        gwas_g = get(load(gwfi))
        gwas_gl = gwas_g[,grepl("none",colnames(gwas_g))]
        mx = mx[,names(gwas_gl)[gwas_gl<pthres]]
        
        if (bin=="01") mx[mx==2] = 1
        if (bin=="12") mx[mx==0] = 1
      }
      mx = mx[apply(mx,1,function(y) any(!is.na(y))), 
              apply(mx,2,function(y) sum(!is.na(y))>(good_na*length(y)))]
      if (x=="gene") mx = mx[,apply(mx,2,function(x) min(table(x))>good_col & length(table(x))>1)]
      if (!is.null(dim(meta_col0[[x]]))) meta_col = meta_col0[[x]][match(colnames(mx),meta_col0[[x]][,id_col]),]
      if (x=="gene") colnames(mx) = meta_col[,"dbSNP"]
      if (x=="rna") colnames(mx) = meta_col[,"symbol"]
      if (x=="rna" | x=="gene") {
        dupind = duplicated(colnames(mx), fromLast=F) | duplicated(colnames(mx), fromLast=T)
        colnames(mx)[dupind] = paste0(colnames(mx)[dupind], "_", meta_col[dupind,id_col])
      }
      return(mx)
    })
    names(X) = names(m) = names(m0)
    # prepare number of pls-da factors you want for each data
    design = matrix(0, nrow = length(X), ncol = length(X))
    
    
    
    
    for (tune_ in c("","tune")) {
      
      keepX = lapply(names(X), function(x) {
        if (x=="cell") return(rep(2, ncomp))
        if (x=="rna" | x=="gene") return(rep(10, ncomp))
        if (x=="metab") return(rep(5, ncomp))
      })
      names(keepX) = names(m0)
      
      if (tune_=="tune") {
        # definition of the keepX value to be tested for each block mRNA miRNA and protein
        # names of test.keepX must match the names of 'data'
        test.keepX = lapply(names(X), function(x) seq(5,20,5))
        names(test.keepX) = names(X)
        
        # the following may take some time to run, note that for through tuning
        # nrepeat should be > 1
        tune = tune.block.splsda(X = X, Y = Y,
                                 ncomp = ncomp, test.keepX = test.keepX, design = design, nrepeat = 1,
                                 cpus = no_cores, validation = "Mfold")
        
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
      result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                            keepX = keepX, design = design,
                            mode = "regression") # PLS regression ("regression"), PLS canonical analysis ("canonical"), redundancy analysis ("invariant") and the classical PLS algorithm ("classic")
      # get variables that have >0 weights
      save(result, file=paste0(pname,ifelse(tune_=="","","_"),tune_,".Rdata"))
      feat1 = lapply(result$loadings, function(x) apply(x, 2, function(i) names(i)[which(i != 0)]))
      
      pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"sampleplot.pdf"), width = 9, height = 5)
      try ({
        plotDiablo(result, ncomp = 1, groupOrder = c("DR", "ER"))
        plotDiablo3(result, ncomp = 1, groupOrder = c("DR", "ER"))
      })
      dev.off()
      
      pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"circosplot.pdf"), width = 7)
      try ({
        circosPlot(result, showIntraLinks = F, cutoff=.5)#, corThreshold = 0.5, cex.label = 0.5)
        circosPlot_diabloModif(result, corThreshold = 0.5, cex.label = 0.5, showIntraLinks = F)
      })
      dev.off()
      
      pdf(paste0(pname,"_",tune_,ifelse(tune_=="","","_"),"heatplot.pdf"), width = 7)
      heatmap_diablo(result, margins = c(2, 12))
      # cimDiablo(result, margins = c(2, 12))
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
          if (feat=="avg") {
            predictScores = Reduce("+", lapply(cv$predict$nrep1, function(i) i[[ncomp]]))/length(X)
          } else {
            predictScores = cv$predict$nrep1[[feat]][[ncomp]] #average prediction
          }
          
          roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot=F, percent = F, na.rm =T, direction = "<")   ## had the wrong direction!!
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
          if (feat=="rna") {
            hk.known = meta_col0[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col0[[feat]][,"symbol"]),]
          } else if (feat=="gene") {
            hk.known = meta_col0[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col0[[feat]][,"dbSNP"]),]
          } else {
            hk.known = as.character(feat1[[feat]])
          }
          write.csv(hk.known, paste0(pname, "_",tune_,ifelse(tune_=="","","_"), feat, "_loading.csv"))
          
          
          # save perf
          biomarkerGenescv = unlist(lapply(strsplit(names(unlist(cv$features$stable$nrep1[[feat]])), "\\."), function(i) i[2]))
          write.csv(biomarkerGenescv, paste0(pname, "_",tune_,ifelse(tune_=="","","_"),feat,"_perf.csv"))
          
          
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
          
        } #feat
        dev.off()
        
      } #test
    } #tune_
  } #bin
} #feat_type
