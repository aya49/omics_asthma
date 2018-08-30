## input: features & meta_file
## output: association study p values
## aya43@sfu.ca
## created 20180524



## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
libr(pkgs())

no_cores = detectCores()-3 #number of cores to use for parallel processing
registerDoMC(no_cores)


## options
overwrite = T
writecsv = T #write results as csv on top of Rdata?

pthres = .01
padjust = p.adjust.methods
padjust = padjust[!padjust%in%c("none")]

scale_cont = F

tests = c("chi2", "lmbayes")
tests_cont = c("lmbayes")
tests_cate = c("chi2")

printmaxrows = 200 #if total number of features under this value, just print everything on csv

categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""

cid_col = "probe"

#plot size
width = 800
height = 400


## features and indices
feat_types = feat_types_annot



## calculate p values ---------------------------------------------
start = Sys.time()

for (feat_type in sort(feat_types,decreasing=T)) { start1 = Sys.time(); try({
  cat("\n", feat_type, sep="")
  
  # load m0 using: meta_file0, feat_dir, feat_type, col_inds0 --> m0, col_inds
  source(paste0(root, "/code/_func-asthma_m0-load.R")) 
  
  class_coli = 1 # response
  # for (class_coli in 1:length(class_col)) {
  for (file_ind_n in names(file_inds)) {
    for (col_ind_n in names(col_inds)) {
      # for (f1_bin in f1_bins) {
      
      # prep m, meta_file, meta_col: m0 class_coli, class_cols, file_ind_n, file_inds
      source(paste0(root, "/code/_func-asthma_m-trim.R")); if (nextm) next()
      if (scale_cont) m[,m_cont_col] = scale(as.numeric(as.matrix(m[,m_cont_col])))
      meta_file = meta_file0[match(rownames(m),meta_file0[,id_col]),]
      if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiment
      
      # meta_col = meta_col0[match(colnames(m),meta_col0[,cid_col]),]
      class_names = levels(factor(meta_file[,class_col]))
      class = as.numeric(factor(meta_file[,class_col]))
      coni = class%in%grep(control,class_names)
      
      # filename; 
      pname0 = paste0(gwas_dir,"/",feat_type,"-",file_ind_n,"X",col_ind_n,"_class-",class_col,"_", paste0(names(table(meta_file[,class_col])),table(meta_file[,class_col]), collapse="v"))
      # pname0 = paste0(gwas_dir,"/",feat_type,ifelse(f1_bin!="",".",""), f1_bin,"-",file_ind_n,"X",col_ind_n,"_class-",class_col,"_", paste0(names(table(meta_file[,class_col])),table(meta_file[,class_col]), collapse="v"))
      
      
      loop_ind = loop_ind_f(1:ncol(m),no_cores)
      design_resp <- model.matrix(~factor(meta_file[,class_col], unique(meta_file[,class_col])))
      is_cont = sum(!m_cont_col)<(ncol(m)*good_catecont)
      is_cate = sum(m_cont_col)<(ncol(m)*good_catecont)
      
      # get differential expression
      fold11 = NULL
      if (is_cont) {
        fold11 = unlist(foreach(ii = loop_ind) %dopar% { 
          foldii = apply(m[,ii,drop=F],2,function(x) {
            ind = x!=-1 & !is.na(x)
            experx = mean(x[ind & !coni]); if (experx==0) experx = 1
            contrx = mean(x[ind &  coni]); if (contrx==0) contrx = 1
            if (sum(x[ind]<0)>0) return( experx - contrx )
            return(log10( experx / contrx ))
          })
        })
      }
      
      for (test in tests) {
        # filename; overwrite?
        pname = paste0(pname0,"_test-",test)
        # if (!file.exists(paste0(pname,".Rdata"))) break
        # next
        
        if ((test%in%tests_cont & is_cate) |
            (test%in%tests_cate & is_cont)) next()
        if (!overwrite & file.exists(paste0(pname,".Rdata"))) next()
        
        ## p value calculation; if Xall already exists, just use that and recalculate adjustment
        cpname = paste0(gsub(paste0("X",col_ind_n),"Xall",pname),".Rdata")
        if (col_ind_n!="all" & file.exists(cpname)) {
          pvalt_temp0 = get(load(cpname))
          pvalt_temp = pvalt_temp0[rownames(pvalt_temp0)%in%colnames(m),,drop=F]
          rownames(pvalt_temp) = rownames(pvalt_temp0)[rownames(pvalt_temp0)%in%colnames(m)]
          pval11 = pvalt_temp[,!grepl(paste(padjust,collapse="|"),colnames(pvalt_temp)),drop=F]
        } else {
          
          if (test=="lmbayes") {
            mfit = t(as.matrix(apply(m,2,as.numeric)))
            fit = NULL
            try({ fit <- eBayes(lmFit(mfit, design_resp)) })
            if (is.null(fit)) next()
            top <- topTable(fit, coef=2, adjust.method="none", n=nrow(fit), sort.by="none")
            pval11 = data.frame(lmbayes_p_none=top$P.Value, logoddsDE=top$B)
            rownames(pval11) = colnames(m)
            
            # col_datasetLabels <- rep(NA, length(datasetLabels))
            # col_datasetLabels[datasetLabels == "UCSC genes"] <- "#66C2A5"
            # col_datasetLabels[datasetLabels == "UCSC gene-isoforms"] <- "#FC8D62"
            # col_datasetLabels[datasetLabels == "Ensembl"] <- "#8DA0CB"
            # col_datasetLabels[datasetLabels == "Trinity"] <- "#E78AC3"
            png(file=paste0(pname,"_logfc.vavg.png"), width=width, height=height*2)
            par(mfrow=c(1,2))
            plot(top$logFC ~ top$AveExpr, pch = 19) #col = col_datasetLabels, 
            abline(h = c(-0.15, 0.15), lty = 2)
            # points(top$logFC ~ top$AveExpr, col = 1, pch = 21)
            # legend("topright", c("UCSC genes", "UCSC gene-isoforms", "Ensembl", "Trinity"),
            #        col = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"), pch = 19, bty = "n")
            graphics.off()
            
            
            # df <- data.frame(FC = top$logFC,
            #                  AveExpr = top$AveExpr)
            # , datasets = datasetLabels)
            # colPalette <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854")
            # png(paste0(pname, "_avg.v.fc.png"), width = 6, height = 4)
            # ggplot(df, aes(x = AveExpr, y = FC)) + #color=datasets
            #   geom_point(size = 3) +
            #   geom_hline(yintercept = -0.15, lty=2) + geom_hline(yintercept = 0.1, lty=2) +
            #   customTheme(sizeStripFont=15, xAngle=0, hjust = 0.5, vjust = 0.5, xSize=10, ySize=15, xAxisSize=15, yAxisSize=15) +
            #   ylab(expression("log"[2]~"fold-change (DR minus ER)")) + xlab(expression("Average Expression log"[2]~"cpm")) +
            #   # scale_color_manual(values=colPalette, name = "Datasets") +
            #   annotate("text",x=12.5,y=0,label="House-keepers", size = 3) +
            #   theme(legend.position = c(0.83, 0.8))
            # graphics.off()
            
          }
          
          if (test=="chi2") {
            m_unique = as.character(sort(unique(as.numeric(m))))
            pval11 = foreach(ii = loop_ind, .combine="cbind") %dopar% { 
              pvalii = apply(m[,ii,drop=F],2,function(x) {
                ind = !is.na(x)
                chi2 = chisq.test(table(x[ind],class[ind]))
                p = chi2$p.value
                ctable_ = chi2$observed
                ctable = NULL
                if (identical(dim(ctable_), c(length(class_names),length(m_unique)))) {
                  ctable = ctable_
                } else {
                  for (colno in m_unique) {
                    if (colno%in%colnames(ctable_)) {
                      ctable = cbind(ctable, ctable_[,colno])
                    } else {
                      ctable = cbind(ctable, rep(0,length(class_names)))
                    }
                  }
                }
                return(append(p, as.vector(ctable)))
              })
            }
            pval11 = as.data.frame(t(pval11))
            colnames(pval11) = append("chi2_p_none", apply(expand.grid(class_names, m_unique), 1, paste, collapse="_") )
            rownames(pval11) = colnames(m)
          }
          
        }
        pvalt = pval11[,grepl("_none$", colnames(pval11))]
        
        ## p value adjust
        pval1 = foreach(padj = padjust, .combine="cbind") %dopar% {
          # pval[[paste0(pvaln,"_",padj)]] = p.adjust(pvalt,method=padj)
          return(p.adjust(pvalt,method=padj))
        }
        colnames(pval1) = paste0(test,"_p_",padjust)
        pv_table = as.matrix(cbind(pval11, pval1))
        rownames(pv_table) = rownames(pval11)
        pv_table = as.data.frame(pv_table[
          apply(pv_table[,grepl("_p_",colnames(pv_table))],1,function(x) any(x<1)),
          apply(pv_table,2,function(x) any(x<1)), drop=F])
        
        if (is_cont & !is.null(fold11)) pv_table$log10de = fold11[rownames(pv_table)]
        
        rown = rownames(m)
        write.csv(rown,file=paste0(pname,"_id.csv"),row.names=F)
        save(pv_table,file=paste0(pname,".Rdata"))
        
        # } #bins 
      } #test
    } #col_ind
  } #row_ind
  # } #class_col
}); time_output(start1) } #feat_type

time_output(start)








