## input: features & meta_file
## output: p vaue features testing significance of ER DR correlation
## aya43@sfu.ca
## created 20180524



## logistics
root = "~/projects/asthma"; commandArgs <- function(...) root  # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/_func-classifiers.R"))
libr(append(pkgs(),c("qqman","traseR","FunciSNP","SeqGSEA","DOSE")))

no_cores = detectCores()-2 #number of cores to use for parallel processing
registerDoMC(no_cores)

data(taSNP) #subset of below
data(taSNPLD)

taSNP0 = as.data.frame(taSNP)
write.csv(taSNP0, paste0(gwas_dir, "/taSNP.csv"))
taSNPLD0 = as.data.frame(taSNPLD)
write.csv(taSNPLD0, paste0(gwas_dir, "/taSNPLD.csv"))


## options

overwrite = T
writecsv = T #write results as csv on top of Rdata?

pthres = .01
padjust = p.adjust.methods
padjust = padjust[!padjust%in%c("none")]
val_max1t = 5e-8 #comment out if just want default, large large p value threshold guideline value

tests = c("chi2", "lmbayes")
tests_cont = c("lmbayes")
tests_cate = c("chi2")

cont_col = 15 #if a column has more than this unique values, it's continuous
printmaxrows = 200 #if total number of features under this value, just print everything on csv

compare_classes = c("goodppl","flipperdr") #do a crossover of features from these classes

categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""

cid_col = "probe"

#plot size
width = 800
height = 400

#plot number of top p values
reg_no = 16


## features and indices
feat_types = feat_types_annot




## calculate p values ---------------------------------------------
start = Sys.time()

for (feat_type in sort(feat_types,decreasing=T)) {
  start1 = Sys.time()
  try({
    cat("\n", feat_type, sep="")
    
    m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
    m0 = m0[rownames(m0)%in%meta_file0[,id_col],]
    # if (sum(colnames(m0)%in%meta_file0[,id_col])==ncol(m0)) m0 = t(m0)
    
    col_inds_i = sapply(names(col_inds0), function(x) grepl(x,feat_type))
    if (sum(col_inds_i)==0) {
      col_inds = list(all=c(""))
    } else {
      col_inds = col_inds0[col_inds_i]
    }
    
    
    for (class_coli in 1:length(class_cols)) {
      for (file_ind_n in names(file_inds)) {
        for (col_ind_n in names(col_inds)) {
          
          class_col = class_cols[class_coli]
          control = controls[class_coli]
          
          file_ind = file_inds[[file_ind_n]]
          if (file_ind_n!="all" & all(rownames(m0)%in%file_ind)) next()
          if (file_ind_n=="all") file_ind = rownames(m0)
          # file_ind_n = paste0("-",file_ind_n)
          
          ## ensure there is a class label for each sample
          # file_ind = file_ind[!is.na(meta_file0[match(file_ind,meta_file0[,id_col]),class_col])]
          
          col_ind = col_inds[[col_ind_n]]
          if (file_ind_n!="all" & all(colnames(m0)%in%col_ind)) next()
          if (col_ind_n=="all") col_ind = colnames(m0)
          # col_ind_n = paste0(".",col_ind_n)
          
          m = as.matrix(m0)[rownames(m0)%in%file_ind, colnames(m0)%in%col_ind]
          class(m) = "numeric"
          m = m[apply(m,1,function(x) any(!is.na(x))), 
                apply(m,2,function(x) sum(!is.na(x))>(good_na*length(x)))]
          if (grepl("genotype",feat_type)) m = m[,apply(m,2,function(x) min(table(x))>good_col & length(table(x))>1)]
          if (any(dim(m)==0) | sum(!is.na(m))<sum(is.na(m))) {
            cat(" (empty or too many NA's skipped!)")
            next()
          }
          
          m_del_col = apply(m,2,function(x) length(unique(x[!is.na(x)]))>1)
          m = m[,m_del_col]
          
          m_cont_col = apply(m,2,function(x) length(unique(x))>cont_col)
          if (scale_cont) m[,m_cont_col] = scale(as.numeric(as.matrix(m[,m_cont_col])))
          
          if (sum(m_del_col)==0) next()
          m[is.na(m)] = -1
          
          cat(" (",file_ind_n, " x ",col_ind_n,"; ",nrow(m)," x ",ncol(m),") ", sep="")
          
          meta_file = meta_file0[match(rownames(m),meta_file0[,id_col]),]
          if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiments
          
          
          # meta_col = meta_col0[match(colnames(m),meta_col0[,cid_col]),]
          class_names = levels(factor(meta_file[,class_col]))
          
          class = as.numeric(factor(meta_file[,class_col]))
          # class = meta_file[,class_col]
          
          # for (i in 1:ncol(m)) {
          #   x = m[,i]
          #   ind = !is.na(x)
          #   a = chisq.test(class[ind],as.numeric(x)[ind])$p.value
          #   cat(i," ")
          # }
          
          for (test in tests) {
            pname = paste0(gwas_dir,"/",feat_type,"-",file_ind_n,"X",col_ind_n,"_class-",class_col,"_test-",test)
            
            if ((test%in%tests_cont & sum(m_cont_col)<(ncol(m)*good_catecont)) |
                (test%in%tests_cate & sum(!m_cont_col)<(ncol(m)*good_catecont))) next()
            if (!overwrite & file.exists(paste0(pname,".Rdata"))) next()
            
            ## p value calculation
            cpname = paste0(gsub(paste0("X",col_ind_n),"Xall",pname),".Rdata")
            if (col_ind_n!="all" & file.exists(cpname)) {
              pvalt_temp0 = get(load(cpname))
              pvalt_temp = pvalt_temp0[rownames(pvalt_temp0)%in%colnames(m),,drop=F]
              rownames(pvalt_temp) = rownames(pvalt_temp0)[rownames(pvalt_temp0)%in%colnames(m)]
              pval11 = pvalt_temp[,!grepl(paste(padjust,collapse="|"),colnames(pvalt_temp)),drop=F]
            } else {
              loop_ind = loop_ind_f(1:ncol(m),no_cores)
              
              if (test=="chi2") {
                m_unique = as.character(sort(unique(as.numeric(m))))
                pval11 = foreach(ii = loop_ind, .combine="cbind") %dopar% { 
                  pvalii = apply(m[,ii],2,function(x) {
                    ind = !is.na(x)
                    chi2 = chisq.test(class[ind],x[ind])
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
                  # for (xi in ii) {
                  #   x = m[,xi]
                  #   ind = !is.na(x)
                  #   a = ifelse(length(unique(x[ind]))==1, 1, chisq.test(class[ind],x[ind])$p.value)
                  # }
                  # return(pvalii)
                }
                pval11 = as.data.frame(t(pval11))
                colnames(pval11) = append("chi2_p_none", apply(expand.grid(class_names, m_unique), 1, paste, collapse="_") )
                rownames(pval11) = colnames(m)
              }
              
              if (test=="lmbayes") {
                design_resp <- model.matrix(~factor(meta_file[,class_col], unique(meta_file[,class_col])))
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
            pv_table = pv_table[apply(pv_table[,grepl("_p_",colnames(pv_table))],1,function(x) any(x<1)),
                                apply(pv_table,2,function(x) any(x<1)), drop=F]
            
            save(pv_table,file=paste0(pname,".Rdata"))
            
            
            # pval1$chi2 = apply(m,2,function(x) { 
            #   ind = !is.na(x)
            #   return(ifelse(length(unique(x[ind]))==1, 1, chisq.test(class[ind],x[ind])$p.value)) 
            # })
            
            # pval1$mi = apply(m,2,function(x) {
            #   ind = !is.na(x)
            #   return(ifelse(length(unique(x[ind]))==1, 1, mcnemar.test(factor(class[ind]),factor(x[ind]))))
            # })
            
          } #test
        } #col_ind
      } #row_ind
    } #class_col
  })
  time_output(start1)
} #feat_type

time_output(start)








#### list out grid for each gwas ###

## plot / enrichment --------------------------------------
start = Sys.time()

# enrichment data
data(CEU) #system("plink -file plink.file -freq -out chr"); to retrieve all SNPs withcorresponding MAF (minor allele frequency) from the CEU vcf files downloaded


pval_paths = sort(gsub(".Rdata","",list.files(gwas_dir, pattern=".Rdata", full.names=T)))
pval_filenames = fileNames(pval_paths)
feat_types = sapply(str_split(pval_filenames,"[-]"), function(x) x[1])

for (pi in 1:length(pval_paths)) { try({
  
  pval_path = pval_paths[pi]
  feat_type = feat_types[pi]
  
  cat("\n",pval_path)
  start1 = Sys.time()
  
  pval0 = get(load(paste0(pval_path,".Rdata")))
  
  m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
  if (sum(colnames(m0)%in%meta_file0[,id_col])==ncol(m0)) m0 = t(m0)
  
  pval_filename = str_split(pval_filenames[pi],"[-]")[[1]]
  class_col = gsub("class-","", str_extract(pval_filenames[pi], "class-[a-zA-Z]+"))
  file_ind_n = gsub("[-]|X","", str_extract(pval_filenames[pi], "-[a-z]+X"))
  m = m0[,rownames(pval0)]
  if (!file_ind_n=="all") m = m[rownames(m)%in%file_inds[[file_ind_n]],]
  meta_file = meta_file0[match(rownames(m), meta_file0[,id_col]),]
  class = as.numeric(factor(meta_file[,class_col]))
  class_names = levels(factor(meta_file[,class_col]))
  
  meta_col0 = meta_col = NULL
  meta_col_path = paste0(meta_col_dir,"-",strsplit(feat_type,"[.]")[[1]][1],".Rdata")
  mcf = file.exists(meta_col_path)
  if (mcf) {
    meta_col0 = as.data.frame(get(load(meta_col_path)))
    meta_col = meta_col0[match(rownames(pval0),meta_col0[,id_col]),]
    
    # symb = sapply(str_extract_all( meta_col0$`Associated Gene`, "ENST[0-9]+" ), function(x) paste(grch38_dt$symbol[match(x,grch38_dt$enstxp)],"_"))
  } 
  
  
  ## save as csv
  pv_table = pval0[apply(pval0,1,function(x) any(x<1)),apply(pval0,2,function(x) !all(x==1))]
  if (all(dim(pv_table)>0)) {
    # if (is.null(dim(pv_table))) pv_table = matrix(pv_table,nrow=1)
    if (mcf) pv_table = cbind(pv_table,meta_col[match(rownames(pv_table),meta_col[,id_col]),])
    if (grepl("genotype",feat_type)) pv_table = cbind(taSNP0[match(pv_table$dbSNP, taSNP0$SNP_ID), c("Trait", "Trait_Class")], pv_table)
    pv_table = pv_table[order(pv_table[,grepl("none",colnames(pv_table))]),]
    
    pval000 = pv_table
    if (sum(pv_table[,grepl("none",colnames(pv_table))]<.05)>1 & nrow(pv_table)>printmaxrows) 
      pval000 = pv_table[pv_table[,grepl("none",colnames(pv_table))]<.05,,drop=F]
    pval000 = pval000[order(pval000[,grepl("none",colnames(pval000))]),]
    
    write.csv(pval000, file=paste0(pval_path,".csv"))
  }
  
  
  
  
  pvalt = pval0[,grepl("_none$", colnames(pval0))]
  pvalt1 = sort(pvalt, decreasing=F)
  nplots = min(reg_no, sum(pvalt<pthres))
  if(nplots>0) {
    ndim_r = ceiling(sqrt(nplots))
    
    png(file=paste0(pval_path,"_reg.png"), width=300*3, height=200*3)
    par(mfrow=c(ndim_r,ndim_r))
    for (si in 1:nplots) {
      boxplot(m[,si]~class, main=paste0(test," unadj p: ", round(pvalt1[si],4)),
              xlab=class_col, ylab=colnames(m)[si], xaxt="n")
      axis(1, at=1:length(unique(class[!is.na(class)])), labels=levels(factor(meta_file[,class_col])))
      points(jitter(class), m[,si], pch=16, col="blue")
    }
    graphics.off()
  }
  
  
  
  
  
  
  ## PLOT!
  prow = sum(grepl("_p_",colnames(pval0)))
  pcol = 1#length(pval)/prow
  
  # qq plot :( not working)
  png(paste0(pval_path,"_qq.png"), width=pcol*width, height=prow*height)
  par(mfcol=c(prow,pcol))
  for (pvaln in sort(colnames(pval0)[grepl("_p_",colnames(pval0))])) { try({
    pvalt = pval0[,pvaln]
    # if (sum(pvalt>val_max1t)>5) val_max1t = 5e-6
    pv_ind = pvalt<1
    if (sum(pv_ind)==0) next()
    
    qq(pvalt)
  }) }
  graphics.off()
  
  # manhattan plot
  if (mcf) {
    
    png(paste0(pval_path,"_manhattan.png"), width=pcol*width, height=prow*height)
    par(mfcol=c(prow,pcol))
    
    for (pvaln in sort(colnames(pval0)[grepl("_p_",colnames(pval0))])) { try({
      pvalt = pval0[,pvaln]
      # if (sum(pvalt>val_max1t)>5) val_max1t = 5e-6
      pv_ind = pvalt<1
      if (sum(pv_ind)==0) next()
      
      # chrom = meta_col$chromosome
      # chrom[chrom=="X"]="23"
      # chrom[chrom=="Y"]="24"
      # chrom[chrom=="XY"]="25"
      # chrom[chrom=="MT"]="26"
      # 
      # manhattan(data.frame(P=pvalt, BP=meta_col$pos_phys, CHR=as.numeric(chrom), SNP=meta_col$dbSNP))
      
      if (grepl("genotype",feat_type)) {
        pos = meta_col[,"pos_phys"]
        chr = meta_col[,"chromosome"]
      }
      if (grepl("rna",feat_type)) {
        pos = meta_col[,"start"]
        chr = meta_col[,"chr"]
      }
      plotind = !is.na(chr) & nchar(chr)<3
      
      try({
        manhattan_plot(val=-log(pvalt[plotind],10), pos=pos[plotind], chrom=chr[plotind], 
                       val_thres=-log(pthres,10), val_max2=-log(1e-4,10), val_max1=-log(val_max1t,10),
                       main=paste0(pvaln," sig=", sum(pvalt < pthres)),
                       lines=c(-log(.01,10),-log(.025,10),-log(.05,10),-log(.1,10)))
      })
      
    }) }
    graphics.off()
    # require(gridExtra)
    # grid.arrange(plots,ncol=pcol,nrow=prow)
    # plotsa = marrangeGrob(grob=plots,ncol=pcol,nrow=prow)
    # ggexport(plotsa, filename = paste0(gwas_dir,"/snp-file-genotype_norm-NA.png"))
    # ggsave(file=paste0(gwas_dir,"/snp-file-genotype_norm-NA.png"), plotsa)
    
  }
  
  
  
  
  # snp enrichment: traseR
  if (mcf & grepl("genotype",feat_type) & grepl("Xall",pval_path)) {
    
    png(paste0(pval_path,"_traseR.png"), width=width, height=height)
    # par(mfcol=c(prow,pcol))
    
    pvaln = sort(colnames(pval0)[grepl("none", colnames(pval0))])
    pvalt = pval0[,pvaln]
    # if (sum(pvalt>val_max1t)>5) val_max1t = 5e-6
    
    pv_ind = pvalt<pthres & meta_col[,"dbSNP"]%in%taSNPLD0$SNP_ID
    pv_ind1 = pvalt<pthres & meta_col[,"dbSNP"]%in%taSNP0$SNP_ID
    if (sum(pv_ind)==0) next()
    pvalt_pv = pvalt[pv_ind]
    pvalt_pv1 = pvalt[pv_ind1]
    meta_col_pv = meta_col[pv_ind,]
    meta_col_pv1 = meta_col[pv_ind1,]
    taSNPLD0_pv = taSNPLD0[match(meta_col_pv$dbSNP,taSNPLD0$SNP_ID),]
    taSNP0_pv = taSNP0[match(meta_col_pv1$dbSNP,taSNP0$SNP_ID),]
    
    # chr, start, end
    snps = taSNPLD0_pv[,c("seqnames","start","end")]#, pvalt_pv[taSNPLD0_ind])
    snps$end = snps$end + 1
    snps1 = taSNP0_pv[,c("seqnames","start","end")]
    snps1$end = snps1$end + 1
    colnames(snps)[c(1)] = colnames(snps1)[c(1)] = c("chr")
    
    plotContext(snpdb=taSNP, region=snps1)#, pvalue=pvalt[pv_ind1])
    mtext(side=1, text=paste0("sig snps in SNP db: ", nrow(taSNP0_pv), " | + lin diseq ", nrow(taSNPLD0_pv), "\n",
                              "sig snps w/ p <",pthres," : ", sum(pvalt<pthres), "\n",
                              "total snps: ", nrow(pval0)))
    
    # binomial
    x1 = traseR(snpdb=taSNP, region=snps1, snpdb.bg=CEU)
    x2 = traseR(snpdb=taSNPLD, region=snps, snpdb.bg=CEU)
    
    x1$tb1$Trait = sapply(str_split(x1$tb1$Trait,"[(]"), function(x) x[1])
    x2$tb1$Trait = sapply(str_split(x2$tb1$Trait,"[(]"), function(x) x[1])
    
    write.csv(x1$tb1, file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNP.csv"))
    write.csv(x2$tb1, file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNPLD.csv"))
    write.csv(x1$tb2, file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNP_class.csv"))
    write.csv(x2$tb2, file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNPLD_class.csv"))
    
    x1tb1 = x1$tb1; x1tb1$trait_type = "trait";       x1tb1$ta_type = "taSNP"
    x1tb2 = x1$tb2; x1tb2$trait_type = "trait_class"; x1tb2$ta_type = "taSNP"
    x2tb1 = x2$tb1; x2tb1$trait_type = "trait";       x2tb1$ta_type = "taSNPLD"
    x2tb2 = x2$tb2; x2tb2$trait_type = "trait_class"; x2tb2$ta_type = "taSNPLD"
    xdf1 = rbind.fill(list(x1tb1[x1tb1$taSNP.hits>0,],x2tb1[x2tb1$Trait%in%x1tb1$Trait[x1tb1$taSNP.hits>0],]), as.matrix)
    xdf2 = rbind.fill(list(x1tb2,x2tb2), as.matrix)
    
    p1a = ggplot(xdf1) + 
      geom_bar(aes(x=interaction(ta_type,Trait), y=taSNP.num, fill=ta_type), stat='identity') +
      ylab("bar = # of total SNPs in database") + ggtitle("traits with 1+ hits") +
      theme(axis.text.x=element_text(angle=90,hjust=1))
    p1b = ggplot(xdf1) + 
      geom_bar(aes(x=interaction(ta_type,Trait), y=taSNP.hits, fill=ta_type), stat='identity') +
      geom_point(aes(x=interaction(ta_type,Trait), y=-log(as.numeric(p.value,10)), fill=ta_type)) +
      ylab("point = -log10(p.value) & bar = # of hits") +
      theme(axis.text.x=element_text(angle=90,hjust=1))
    
    p2a = ggplot(xdf2) + 
      geom_bar(aes(x=interaction(ta_type,Trait_Class), y=taSNP.num, fill=ta_type), stat='identity') +
      ylab("bar = # of total SNPs in database") + 
      theme(axis.text.x=element_text(angle=90,hjust=1))
    p2b = ggplot(xdf2) + 
      geom_bar(aes(x=interaction(ta_type,Trait_Class), y=taSNP.hits, fill=ta_type), stat='identity') +
      geom_point(aes(x=interaction(ta_type,Trait_Class), y=-log(as.numeric(p.value,10)), fill=ta_type)) +
      ylab("point = -log10(p.value) & bar = # of hits") +
      theme(axis.text.x=element_text(angle=90,hjust=1))
    
    png(file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNP.png"), width=width, height=height*2)
    print(p1a)
    graphics.off()
    png(file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNP-total.png"), width=width, height=height*2)
    print(p1b)
    graphics.off()
    
    png(file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNP_class.png"), width=width, height=height*2)
    print(p2a)
    graphics.off()
    png(file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNP_class-total.png"), width=width, height=height*2)
    print(p2b)
    graphics.off()
  }
  
  
  
  #use another package to do SNP disease enrichment
  if (mcf & grepl("genotype",feat_type) & grepl("Xall",pval_path)) {
    pv_ind2 = pvalt<pthres
    if (sum(pv_ind)==0) next()
    pvalt_pv = pvalt[pv_ind]
    sig_snps = meta_col[meta_col[,id_col]%in%names(pvalt_pv),"dbSNP"]
    edsnps = enrichDGNv(sig_snps)
    if (nrow(engenes@result)>0) {
      edsnps_table = cbind(edsnps@result, t(sapply(edsnps@result$ID, function(x) c(length(edsnps@geneSets[[x]]), paste(edsnps@geneSets[[x]][edsnps@geneSets[[x]]%in%sig_snps], collapse="_")) )))
      
      write.csv(edsnps_table,file=paste0(pval_path,"_dose-snp.csv"))
      png(file=paste0(pval_path,"_dose-snp.png"), width=width, height=height)
      enrichMap(edsnps)
      graphics.off()
    }
    
    
    # # enrichment analysis using hypergeometric test
    # sea00 <- sear(input=sigGenes, type = 'mrna')
    # sea0 = filter(sea00, collection %in% c("TISSUES")) %>% 
    #   arrange(fdr)
    # sea0$geneset <- as.character(sea0$geneset)
    # sea0$geneset[sea0$subcollection == "CNS"] <- "CNS"
    # sea0$geneset <- factor(sea0$geneset)
    # sea <- sea0 %>% dplyr::group_by(subcollection, geneset) %>% 
    #   dplyr::summarise(fdr = mean(fdr))
    # sea$subcollection <- factor(sea$subcollection, 
    #                             levels = toupper(c("cns", "thymus", "tonsils", "spleen", "lymph", "bcells","myeloid", "nkcells","blood","rbcs","neutrophils", "tcells","th1cells","th2cells","tregs","cd8tcells","cd4tcells")))
    # subCol <- sea %>% dplyr::group_by(subcollection) %>% 
    #   dplyr::arrange(subcollection) %>% dplyr::select(subcollection, geneset, fdr) %>% as.data.frame() 
    # subCol$geneset <- factor(as.character(subCol$geneset), levels = unique(as.character(subCol$geneset)))
    
    # #dev.off(); dev.off();
    # cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
    #                brewer.pal(9, "Set1"))
    # pdf(paste0(data_dir, "/results/Figures/Figure2/Figure2c.pdf"), width = 8, height = 5.5)
    # ggplot(subCol, aes(x = geneset, y = -log10(fdr), fill = subcollection)) + 
    #   geom_bar(stat = "identity") +
    #   customTheme(sizeStripFont=10, xAngle=90, hjust=1, vjust=0.5, 
    #               xSize=10, ySize=10, xAxisSize=10, yAxisSize=10) +
    #   ylab(expression("-log"[10]~"(BH-FDR)")) + xlab("") + scale_fill_manual(values=cbPalette, name = "Tissue")
    # dev.off()
    
    
    
    
    
    
    # library(InterMineR)
    # 
    # # get HumanMine instance
    # im.human = initInterMine(listMines()["HumanMine"])
    # 
    # # retrieve data model for HumanMine
    # hsa_model = getModel(im.human)
    # 
    # # all targets of HumanMine enrichment widgets possess InterMine ids
    # subset(hsa_model, type %in% c("Gene", "Protein", "SNP") & child_name == "id")
    # 
    # data("PL_DiabetesGenes")
    # head(PL_DiabetesGenes, 3)
    # 
    # # get Gene.symbol
    # hsa_gene_names = as.character(PL_DiabetesGenes$Gene.symbol)
    # head(hsa_gene_names, 3)
    # 
    # # get Gene.primaryIdentifiers (ENTREZ Gene identifiers)
    # hsa_gene_entrez = as.character(PL_DiabetesGenes$Gene.primaryIdentifier)
    # head(hsa_gene_entrez, 3)
    # 
    # # get all HumanMine widgets
    # human.widgets = as.data.frame(getWidgets(im.human))
    # human.widgets
    # 
    # # get enrichment, gene-related widgets for HumanMine
    # subset(human.widgets, widgetType == 'enrichment' & targets == "SNP")
    # 
    # # Perform enrichment analysis
    # # returns: a data.frame with the results of the enrichment analysis which was performed in the defined InterMine platform. p-values are given after applying the correction algorithm. Count and populationAnnotationCount columns contain the number of genes that belong to each GO term in the list and in the background population respectively.
    # 
    # GO_enrichResult = doEnrichment(
    #   im = im.human,
    #   ids = hsa_gene_entrez,
    #   widget = "go_enrichment_for_gene"
    #   # organism = "Homo sapiens" # optional if results from more than one organism are retrieved
    # )
    # 
    # head(GO_enrichResult$data)
    # dim(GO_enrichResult$data)
    # GO_enrichResult$populationCount # size of the reference background population (populationCount)
    # GO_enrichResult$notAnalysed #number of input features that were not included in the enrichment analysis (notAnalyzed)
    # GO_enrichResult$im # name and url of the Mine (im)
    # 
    # 
    # 
    # 
    # ## get available filter values for Gene Ontology Enrichment widget
    # as.character(subset(human.widgets, name == "go_enrichment_for_gene")$filters)
    # 
    # # Perform enrichment analysis for GO molecular function terms
    # GO_MF_enrichResult = doEnrichment(
    #   im = im.human,
    #   ids = hsa_gene_entrez,
    #   widget = "go_enrichment_for_gene",
    #   filter = "molecular_function")
    # 
    # head(GO_MF_enrichResult$data)
    # 
    # ## multiple testing adjustment
    # 
    # # Perform enrichment analysis for Protein Domains in list of genes
    # PD_FDR_enrichResult = doEnrichment(
    #   im = im.human,
    #   ids = hsa_gene_entrez,
    #   widget = "prot_dom_enrichment_for_gene",
    #   correction = "Benjamini Hochberg"
    # ) 
    # 
    # head(PD_FDR_enrichResult$data)
    # 
    # # Perform enrichment analysis for Protein Domains in list of genes
    # # but without a correction algoritm this time
    # PD_None_enrichResult = doEnrichment(
    #   im = im.human,
    #   ids = hsa_gene_entrez,
    #   widget = "prot_dom_enrichment_for_gene",
    #   correction = "None"
    # )
    # 
    # head(PD_None_enrichResult$data)
    # 
    # ## visualize
    # library(GeneAnswers)
    # # convert InterMineR Gene Ontology Enrichment analysis results to GeneAnswers object
    # geneanswer_object = convertToGeneAnswers(
    #   
    #   # assign with doEnrichment result:
    #   enrichmentResult = GO_enrichResult,
    #   
    #   # assign with list of genes:
    #   geneInput = data.frame(GeneID = as.character(hsa_gene_entrez), 
    #                          stringsAsFactors = FALSE),
    #   
    #   # assign with the type of gene identifier
    #   # in our example we use Gene.primaryIdentifier values (ENTREZ IDs):
    #   geneInputType = "Gene.primaryIdentifier",
    #   
    #   # assign with Bioconductor annotation package:
    #   annLib = 'org.Hs.eg.db',
    #   
    #   # assign with annotation category type
    #   # in our example we use Gene Ontology (GO) terms:
    #   categoryType = "GO"
    #   
    #   #optional define manually if 'enrichIdentifier' is missing from getWidgets:
    #   #enrichCategoryChildName = "Gene.goAnnotation.ontologyTerm.parents.identifier"
    # )
    # 
    # class(geneanswer_object)
    # summary(geneanswer_object)
    # 
    # # convert to GeneAnswers results for GO terms associated with molecular function
    # geneanswer_MF_object = convertToGeneAnswers(
    #   enrichmentResult = GO_MF_enrichResult,
    #   geneInput = data.frame(GeneID = as.character(hsa_gene_entrez), 
    #                          stringsAsFactors = FALSE),
    #   geneInputType = "Gene.primaryIdentifier",
    #   annLib = 'org.Hs.eg.db',
    #   categoryType = "GO.MF"
    #   #enrichCategoryChildName = "Gene.goAnnotation.ontologyTerm.parents.identifier"
    # )
    # 
    # class(geneanswer_MF_object)
    # 
    # 
    # 
    # 
    # 
    # 
    # # GeneAnswers barplot
    # geneAnswersChartPlots(geneanswer_object, 
    #                       chartType='barPlot',
    #                       sortBy = 'correctedPvalue',
    #                       top = 5)
    # # generate concept-gene network
    # geneAnswersConceptNet(geneanswer_object, 
    #                       colorValueColumn=NULL,
    #                       centroidSize='correctedPvalue', 
    #                       output='interactive',
    #                       catTerm = FALSE,
    #                       catID = FALSE,
    #                       showCats = 1:5)
    
    
    graphics.off()
  }
  
  
  
  # gene enrichment for snp
  sig_symbol = c()
  if (mcf & (grepl("genotype",feat_type) | grepl("rna",feat_type)) & grepl("Xall",pval_path)) {
    if (grepl("genotype",feat_type)) {
      sig_enst = 
        str_extract_all(meta_col$`Associated Gene`[ pv_table[,grepl("none",colnames(pv_table))]<pthres ], "ENST[0-9]+")
      sig_enst = unlist(lapply(sig_enst, function(x) x[!duplicated(x)]))
      try({ sig_entr = grch38_dt$entrez[match(sig_enst,grch38_dt$enstxp)] })
      
    } else if (grepl("rna",feat_type)) {
      isensid = any(grepl("ENS", meta_col[1,]))
      
      idd_col = ifelse (isensid, grep("ENS", meta_col[1,]), id_col)
      
      sig_id = meta_col[ pv_table[,grepl("none",colnames(pv_table))]<pthres, idd_col]
      sig_enst = sig_id
      sig_entr = NULL
      if(isensid)  {
        sig_enst = sapply(strsplit(unlist(str_extract_all(sig_id, "ENSG[0-9]+")),"[.]"), function(x) x[1])
        try({ sig_entr = grch38_dt$entrez[match(sig_enst,grch38_dt$ensgene)] })
      }
      
      if (all(is.na(sig_entr))) {
        sig_enst = unlist(str_extract_all(meta_col$`Associated Gene`[ pv_table[,grepl("none",colnames(pv_table)), idd_col]<pthres ], "ENST[0-9]+"))
        sig_enst = sapply(strsplit(sig_enst,"[.]"), function(x) x[1])
        try({ sig_entr = grch38_dt$entrez[match(sig_enst,grch38_dt$enstxp)] })
      }
    }
    if (length(sig_enst)==0 | all(is.na(sig_enst))) next()
    
    sig_entr = grch38_dt$entrez[match(sig_enst,grch38_dt$enstxp)]
    sig_entr = sig_entr[!is.na(sig_entr)]
    
    engenes = enrichDO(sig_entr)
    if (nrow(engenes@result)>0) {
      engenes_table = cbind(engenes@result, 
                            t(sapply(engenes@result$ID, function(x) 
                              c(length(engenes@geneSets[[x]]), paste(grch38_dt$symbol[grch38_dt$entrez%in%engenes@geneSets[[x]]], collapse="_")) )))
      
      write.csv(engenes_table,file=paste0(pval_path,"_dose-gene.csv"))
      png(file=paste0(pval_path,"_dose-gene.png"), width=width, height=height)
      enrichMap(engenes)
      graphics.off()
    }
    
    
  }
  
  
  time_output(start1)
}) } #pi
time_output(start)









## compare goodppl and flipperDR -----------------------------------

pval_paths = sort(gsub(".Rdata","",list.files(gwas_dir, pattern=".Rdata", full.names=T)))
pval_paths = Reduce("intersect", lapply(compare_classes, function(cc) gsub(cc,"",pval_paths[grepl(cc,pval_paths)]) ))

start = Sys.time()

for (pval_path in pval_paths) { try({
  p0 = lapply(compare_classes, function(x) read.csv(paste0(gsub("-X",paste0("-",x,"X"),pval_path),".csv"), row.names=1))
  ptable = NULL
  for (pvi in 2:length(p0)) 
    ptable = merge(p0[[pvi-1]],p0[[pvi]], by=id_col, suffixes=compare_classes[c("",compare_classes[pvi])])
  ptable1 = ptable[apply(ptable[,grepl("_none",colnames(ptable))],1,function(x) all(x<pthres)),]
  
  write.csv(ptable1, file=paste0(gsub("-X",paste0("-",paste0(compare_classes,collapse="."),"X"),pval_path),".csv"))
  write.csv(ptable, file=paste0(gsub("-X",paste0("-",paste0(compare_classes,collapse="."),"X"),pval_path),"_full.csv"))
  
})}

time_output(start)
