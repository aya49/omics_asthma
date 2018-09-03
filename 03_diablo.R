## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/src/_dirs.R"))
source(paste0(root, "/src/_func.R"))
source(paste0(root, "/src/_func-asthma.R"))
source(paste0(root, "/src/visualizationFunctions.R"))
libr(append(pkgs(),c("pROC", "rafalib", "ROCR", "CellCODE", "GGally", "mixOmics"))) # other plots

# no_cores = 5#detectCores()-3
# registerDoMC(no_cores)


## options
overwrite = T

pthres = .01

# plot
height = 400
width = 500

ncomp = 3 #number of components to use

nfeat_short = 5 #number of features in each component
nfeat_long = 10

scale_cont = F #if values are continuous, scale?

celldiff = get(load(meta_col_rnacells_dir))

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
                       
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   dna="dna01")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna="rnapcgenes.pre",
                                   rna="rnaelements.pre",
                                   dna="dna01")),
                       
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   dna="dna12")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna="rnapcgenes.pre",
                                   rna="rnaelements.pre",
                                   dna="dna12")),
                       
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                                   metab="metab.pre",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna="rnapcgenes.pre",
                                   rna="rnaelements.pre",
                                   cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                                   metab="metab.pre",
                                   dna="dna")),
                       
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                                   metab="metab.pre",
                                   dna="dna01")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna="rnapcgenes.pre",
                                   rna="rnaelements.pre",
                                   cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                                   metab="metab.pre",
                                   dna="dna01")),
                       
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                                   metab="metab.pre",
                                   dna="dna12")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna="rnapcgenes.pre",
                                   rna="rnaelements.pre",
                                   cell="cellseqgenes.pre", # get rid of rnaseq stuff that went into predicting cells
                                   metab="metab.pre",
                                   dna="dna12")),
                       
                       unlist(list(rna="rnaseqstarens.post",
                                   cell="cellseqgenes.post",
                                   metab="metab.post",
                                   dna="dna")),
                       unlist(list(rna="rnaseqstarens.post",
                                   cell="cellseqgenes.post",
                                   metab="metab.post",
                                   dna="dna01")),
                       unlist(list(rna="rnaseqstarens.post",
                                   cell="cellseqgenes.post",
                                   metab="metab.post",
                                   dna="dna12")),
                       
                       unlist(list(rna="rnaseqstarens.diff",
                                   cell="cellseqgenes.diff",
                                   metab="metab.diff",
                                   dna="dna")),
                       
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   cell="cell.pre",
                                   metab="metab.pre",
                                   dna="dna01")),
                       unlist(list(rna="rnaseqstarens.pre",
                                   rna.pc="rnapcgenes.pre",
                                   rna.elements="rnaelements.pre",
                                   cell="cell.pre",
                                   metab="metab.pre",
                                   dna="dna12")),
                       
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

feat_type_sets2 = 
  lapply(feat_types_annot, function(x) { names(x) = feats[sapply(feats,function(y) grepl(y,x))]; x})
feat_type_sets = append(append(feat_type_sets1,feat_type_sets2),feat_types_annots)

feat_types = feat_types_annot

# load and trim all features
source(paste0(root,"/src/_func-asthma_mset0-load.R"))

start = Sys.time()

for (feat_type_set in feat_type_sets) { try ({
  start1 = Sys.time()
  cat("\n",feat_type_set, " ", sep=" ")
  if (!all(feat_type_set%in%names(m0_inds))) next()
  
  feat_name = gsub("[0-9]|[.]pre|[.]post|[.]diff","",feat_type_set)
  
  class_coli = 1
  col_ind_n = "all"
  for (file_ind_n in names(file_inds)) {
    if (!all(sapply(m0_inds[feat_type_set], function(x) any(names(x)%in%file_ind_n)))) next()
    
    ## trim matrices
    source(paste0(root,"/src/_func-asthma_mset-trim.R"))
    if (length(m_set)<1) { cat(" skipped, no features left "); next() }
    if (length(unique(meta_file[,class_col]))<2) { cat(" skipped, no variety in class "); next() } 
    
    # prepare meta data for samples (their classes)
    Y = factor(meta_file[,class_col], levels=c("ER", "DR"))
    names(Y) = meta_file[,id_col]
    
    # get col feature meta
    meta_col_set = NULL
    for (x in 1:length(feat_name)) {
      if (is.null(m_col0s[[feat_name[x]]])) next()
      meta_col_setx = m_col0s[[feat_name[x]]]
      meta_col_set[[feat_type_set[x]]] = 
        meta_col_setx[match(colnames(m_set[[feat_type_set[x]]]), meta_col_setx[,id_col]),]
    }
    cole = feat_type_set%in%names(meta_col_set)
    
    # prepare actual data
    X0 = list()
    for (xi in 1:length(m_set)) {
      x = names(m_set)[xi]
      mx = m_set[[x]]
      if (grepl("dna",x) & dnasigonly) { #trim based on p value
        gwf = list.files(gwas_dir, pattern=".Rdata", full.names=T)
        pn = gwf[grepl(paste0(x,"-",file_ind_n,"Xall"),gwf)]
        gwas_g = get(load(pn))
        gwas_gl = gwas_g[match(colnames(mx),rownames(gwas_g)),grepl("none",colnames(gwas_g))]
        names(gwas_gl) = rownames(gwas_g)[match(colnames(mx),rownames(gwas_g))]
        mx = mx[,names(gwas_gl)[gwas_gl<pthres & !is.na(gwas_gl)]]
        if (!is.null(dim(meta_col_set[[x]]))) meta_col_set[[x]] = meta_col_set[[x]][match(colnames(mx),meta_col_set[[x]][,id_col]),]
      }
      if (grepl("dna",x)) { nonars = !is.na(meta_col_set[[x]][,"dbSNP"])
      colnames(mx)[nonars] = paste0(colnames(mx)[nonars], "_", meta_col_set[[x]][nonars,"dbSNP"])  } 
      if (grepl("rna",x)) { nonasy = !is.na(meta_col_set[[x]][,"symbol"])
      colnames(mx)[nonasy] = paste0(colnames(mx)[nonars], "_", meta_col_set[[x]][nonasy,"symbol"]) } 
      # if (grepl("rna",x) | grepl("dna",x)) {
      #   dupind = duplicated(colnames(mx), fromLast=F) | duplicated(colnames(mx), fromLast=T)
      #   colnames(mx)[dupind] = paste0(colnames(mx)[dupind], "_", meta_col_set[[x]][dupind,id_col])
      # }
      X0[[x]] = mx
    }
    
    # 0 has all the original feat
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
    
    
    ## blocksplsda -------------------------------------
    
    for (constrains in c(0,.5,1)) { # how constrainsed is each feature type with each other)
      for (tune_ in c("","tune")) {
        
        # prepare file name
        pname = paste0(blocksplsda_dir, "/", paste(feat_type_set,collapse="-"), "-",  file_ind_n, "Xall", 
                       "_class-", paste(class_col,collapse="."),"_", 
                       paste0(names(table(meta_file[,class_col])),table(meta_file[,class_col]), collapse="v"),
                       "_pthres-", pthres, ifelse(tune_=="","","_"),tune_,
                       "_constrains-",constrains)
        # pname = paste0(blocksplsda_dir, "/", paste(feat_type_set,collapse="-"), ifelse(f1_bin=="","","."), f1_bin, "-",  file_ind_n, "Xall", "_class-", paste(interested_cols,collapse="."),"_", paste0(names(table(meta_file[,class_col])),table(meta_file[,class_col]), collapse="v"), "_pthres-", pthres,ifelse(tune_=="","","_"),tune_)
        
        if (overwrite & file.exists(paste0(pname,".Rdata"))) next()
        dir.create(pname,showWarnings=F)
        
        
    # prepare how constrainsed each feature is with each other
    design = matrix(constrains, nrow = length(X), ncol = length(X))
    diag(design) = 0
    
        rown = Reduce(intersect,lapply(X,rownames))
        write.csv(rown,file=paste0(pname,"_id.csv"),row.names=F)
        
        # prepare number of pls-da factors you want for each data
        keepX = lapply(names(X), function(x) {
          if (grepl("cell|metab",x)) return(rep(nfeat_short, ncomp))
          if (grepl("rna|dna",x)) return(rep(nfeat_long, ncomp))
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
            tune = tune.splsda(X=X[[1]], Y=Y, ncomp=ncomp, 
                               test.keepX=test.keepX[[1]], nrepeat=1,
                               # cpus=no_cores,
                               validation="Mfold")
            
          } else {
            tune = tune.block.splsda(X=X, Y=Y, ncomp=ncomp, 
                                     test.keepX=test.keepX, design=design, nrepeat=1,
                                     # cpus=no_cores, 
                                     validation="Mfold")
          }
          
          tune$choice.keepX.constrainst # NULL as constrainst = F per default
          tune$choice.keepX
          
          pdf(paste0(pname,"/tune.pdf"))
          par(mfrow=c(2,1))
          plot(tune$error.rate[, "comp1"][order(tune$error.rate[, ncol(tune$error.rate)])] ~ factor(rownames(tune$error.rate)), col=1, ylim=c(0.2, 1), main=as.character(which.min(tune$error.rate[,ncol(tune$error.rate)])))
          points(tune$error.rate[,ncol(tune$error.rate)][order(tune$error.rate[,ncol(tune$error.rate)])] ~ factor(rownames(tune$error.rate)), col=2)
          dev.off()
          
          # use tuned keepX
          keepX = tune$choice.keepX
        }
        
        ## block.splsda = horizontal integration PLS-DA model with a specified number of components per block either by Y or by its position indY in the list of blocks X
        # result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
        #                       keepX = keepX, design = design,
        #                       mode = "regression", bias = T)
        if (length(X)==1) {
          if (typeof(keepX)=="list") keepX=keepX[[1]]
          result = splsda(X=X[[1]], Y=Y, ncomp=ncomp, 
                          keepX=keepX, #design=design,
                          mode="regression")
        } else {
          result = block.splsda(X=X, Y=Y, ncomp=ncomp, 
                                keepX=keepX, design=design,
                                mode="regression") # PLS regression ("regression"), PLS canonical analysis ("canonical"), redundancy analysis ("invariant") and the classical PLS algorithm ("classic")
        }
        
        # get variables that have >0 weights
        save(result, file=paste0(pname,".Rdata"))
        feat1 = lapply(result$loadings, function(x) apply(x, 2, function(i) names(i)[which(i != 0)]))
        
        ## blocksplsda plots -------------------------------------------
        pdf(paste0(pname,"/plots.pdf"), width=7)
        try ({
          # plotDiablo(result, ncomp=1, groupOrder=c("DR", "ER"))
          # plotDiablo3(result, ncomp=1, groupOrder=c("DR", "ER"))
          plotIndiv_diablo(result, ncomp=1, groupOrder=c("DR", "ER"))
        })
        dev.off()
        
        pdf(paste0(pname,"/arrows.pdf"), width=7)
        try ({
          plotArrow(result)
        })
        dev.off()
        
        pdf(paste0(pname,"/circos.pdf"), width=7)
        try ({
          circosPlot(result, showIntraLinks=F, cutoff=.5)#, corThreshold=0.5, cex.label=0.5)
          # circosPlot_diabloModif(result, corThreshold=0.5, cex.label=0.5, showIntraLinks=F) #grep for rna, name in X list must have rna for transcriptome!
        })
        dev.off()
        
        pdf(paste0(pname,"/heatmap.pdf"), width=7)
        try ({
          heatmap_diablo(result, margins=c(2, 12))
          # cim(result, margins=c(2, 12))
          # cimDiablo(result, margins=c(2, 12))
        })
        dev.off()
        
        pdf(paste0(pname,"/var.pdf"), width=7)
        try ({
          plotVar(result)
        })
        dev.off()
        
        pdf(paste0(pname,"/loadings.pdf"), width=7)
        try ({
          plotLoadings(result)
        })
        dev.off()
        
        pdf(paste0(pname,"/network.pdf"), width=7)
        try ({
          network(result) # network: similarity between features is obtained by calculating the sum of the correlations between the original variables and each of the latent components of the model
        })
        dev.off()
        
        weight = result$weights
        write.csv(weight, paste0(pname, "/weights.csv"))
        
        ev = result$explained_variance$Y
        write.csv(ev, paste0(pname, "/ev_class.csv"))
        
        # performance of splsda
        for (test in c("loo","Mfold")) {
          cv = perf(result, validation=test)
          cv$WeightedPredict.error.rate
          # cv$error.rate
          
          png(paste0(pname,"/predictscore-",test,".png"), height=(length(X))*height, width=width)
          par(mfrow=c(length(X)+1,1))
          for (feat in c(names(X),"avg")) {
            
            # plot AUC
            if (length(X)>1) {
              if (feat=="avg") {
                predictScores = Reduce("+", lapply(cv$predict$nrep1, function(i) i[[length(i)]]))/length(X)
              } else {
                predictScores = cv$predict$nrep1[[feat]][[length(cv$predict$nrep1[[feat]])]] #average prediction
              }
            } else {
              if (feat=="avg") next()
              predictScores = cv$predict[[length(cv$predict)]][,,1]
            }
            
            roc.score = roc(response=factor(as.character(Y), levels=c("ER", "DR")), predictor=predictScores[,"DR"], plot=F, percent=F, na.rm=T, direction="<")   ## had the wrong direction!!
            roc.res4=data.frame(Specificity=rev(roc.score$specificities), Sensitivity=rev(roc.score$sensitivities))
            roc.res4$Specificity = 1 - as.numeric(roc.res4$Specificity)
            roc.res4$Sensitivity = as.numeric(roc.res4$Sensitivity)
            # roc.score$auc
            # plot(roc.score, main=paste0(feat," auc=",round(roc.score$auc,4)))
            
            plot(roc.res4$Sensitivity ~ roc.res4$Specificity, type="o", pch=19,
                 xlab="1-Specificity", ylab="Sensitivity", col="#1F78B4",
                 main=paste0(feat," auc=",round(roc.score$auc,4)))
            abline(a=0, b=1, col="gray", lty=2)
            text(x=80, y=20, labels=paste0("AUC = ", roc.score$auc), col="#1F78B4")
            
            if (feat=="avg") next()
            
            
            # save features
            if (grepl("rna",feat)) {
              hk.known = meta_col_set[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col_set[[feat]][,"symbol"]),]
            } else if (feat=="dna") {
              hk.known = meta_col_set[[feat]][match(sapply(strsplit(as.character(feat1[[feat]]), "[.]|_"), function(i) i[1]), meta_col_set[[feat]][,"dbSNP"]),]
            } else {
              hk.known = as.character(feat1[[feat]])
            }
            # write.csv(hk.known, paste0(pname, "/loading_", feat, ".csv"))
            
            loading = result$loadings[[feat]]
            loading = loading[apply(loading,1,function(x) any(x)>1),]
            if (nrow(loading)>0) write.csv(loading, paste0(pname, "/loading_", feat, ".csv"))
            
            comp = result$variates[[feat]]
            write.csv(comp, paste0(pname, "/component_", feat, ".csv"))
            
            ev = result$explained_variance[[feat]]
            write.csv(ev, paste0(pname, "/ev_", feat, ".csv"))
            
            # save perf
            biomarkerGenescv = 
              unlist(lapply(strsplit(names(unlist(
                ifelse (length(X)==1, cv$features$stable, cv$features$stable$nrep1[[feat]])
              )), "[.]"), function(i) i[2]))
            write.csv(biomarkerGenescv, paste0(pname, "/perf_",feat,".csv"))
            
          } #feat
          dev.off()
          
        } #test
        
      } #tune_
    } #constrains
  } #file_ind
}) } #feat_type_set

time_output(start)
