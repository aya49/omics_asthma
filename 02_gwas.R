## input: features & meta_file
## output: p vaue features testing significance of ER DR correlation
## aya43@sfu.ca
## created 20180524



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result")



## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")#,asthma)

feat_dir = paste0(result_dir,"/feat")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir,showWarnings=F)
gwas_dir = paste0(stat_dir,"/gwas"); dir.create(gwas_dir,showWarnings=F)

## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("data.table")
libr("qqman")

libr("traseR")
data(taSNP) #subset of below
data(taSNPLD)
# data(Tcell)
taSNP0 = as.data.frame(taSNP)
taSNPLD0 = as.data.frame(taSNPLD)

libr("FunciSNP")

libr("ggplot2")
libr("entropy")
libr("limma")
libr("foreach")
libr("doMC")
libr("stringr")
libr("Matrix")



## options
options(stringsAsFactors=F)

no_cores = 5#detectCores()-3
registerDoMC(no_cores)

overwrite = F
writecsv = T #write results as csv on top of Rdata?

pthres = .01
padjust = p.adjust.methods
padjust = padjust[!padjust%in%c("none")]
val_max1t = 5e-8 #comment out if just want default, large large p value threshold guideline value

tests = c("chi2", "lmbayes")
tests_cont = c("lmbayes")
tests_cate = c("chi2")

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis
cont_col = 15 #if a column has more than this unique values, it's continuous
scale_cont = F #if a column is continuous, scale it
good_na = .75 #proportion of na more than this, then delete the column in matrix

id_col = "id"
class_cols = "response"
controls = "ER"
categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""

cid_col = "probe"

#plot size
width = 1000
height = 500



## features and indices
feat_types = list.files(feat_dir,full.names=F,pattern=".Rdata")
feat_types = feat_types[!grepl("raw",feat_types)]
feat_types = gsub(".Rdata","",feat_types)
feat_temp = str_split(feat_types,"[.]")
feat_names = sapply(feat_temp, function(x) x[1])
feat_times = sapply(feat_temp, function(x) x[2])

inds_paths = list.files(meta_dir,pattern="_id_",full.names=T)
inds_temp = str_split(gsub(".Rdata","",fileNames(inds_paths)),"_")
inds_types = sapply(inds_temp, function(x) x[3])
inds_names = sapply(inds_temp, function(x) str_split(x[1],"[-]")[[1]][2])
inds_filecol = sapply(inds_temp, function(x) str_split(x[1],"[-]")[[1]][1])

uif = unique(inds_names[inds_filecol=="col"])
col_inds0 = lapply(uif, function(x) list(all=c("")))
names(col_inds0) = uif
file_inds = list(all=c(""))

# col inds separate by feature type, file inds don't
for (i in 1:length(inds_paths)) {
  if (inds_filecol[i]=="col") {
    col_inds0[[inds_names[i]]][[inds_types[i]]] = get(load(inds_paths[i]))
  } else if (inds_filecol[i]=="file") {
    file_inds[[inds_types[i]]] = get(load(inds_paths[i]))
  }
}



## calculate p values
start = Sys.time()

meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
# meta_col0 = as.data.frame(get(load(paste0(meta_col_dir,".Rdata"))))

for (feat_type in feat_types) {
  cat("\n", feat_type, sep="")
  
  m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
  if (sum(colnames(m0)%in%meta_file0[,id_col])==ncol(m0)) m0 = t(m0)
  
  col_inds_i = sapply(names(col_inds0), function(x) grepl(x,feat_type))
  if (sum(col_inds_i)==0) {
    col_inds = list(all=c(""))
  } else {
    col_inds = col_inds0[[col_inds_i]]
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
        
        # ensure there is a class label for each sample
        file_ind = file_ind[!is.na(meta_file0[match(file_ind,meta_file0[,id_col]),class_col])]
        
        col_ind = col_inds[[col_ind_n]]
        if (file_ind_n!="all" & all(colnames(m0)%in%col_ind)) next()
        if (col_ind_n=="all") col_ind = colnames(m0)
        # col_ind_n = paste0(".",col_ind_n)
        
        m = m0[rownames(m0)%in%file_ind, colnames(m0)%in%col_ind]
        m = m[apply(m,1,function(x) any(!is.na(x))), 
              apply(m,2,function(x) sum(is.na(x))<(good_na*length(x)))]
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
        # meta_col = meta_col0[match(colnames(m),meta_col0[,cid_col]),]
        
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
          
          if ((test%in%tests_cont & sum(m_cont_col)==0) |
              (test%in%tests_cate & sum(!m_cont_col)==0) |
              !overwrite & file.exists(paste0(pname,".Rdata"))) next()
          
          ## p value calculation
          cpname = paste0(gsub(paste0("X",col_ind_n),"Xall",pname),".Rdata")
          if (col_ind_n!="all" & file.exists(cpname)) {
            pvalt_temp = get(load(cpname))
            if (!is.null(dim(pvalt_temp))) {
              pvalt_temp = pvalt_temp[rownames(pvalt_temp)%in%colnames(m),]
            } else {
              cat(" Xall a vector...")
              next()
              # # added this because... might have saved pval from an older ver that didn't make tables
              # pvalt_temp1 = data.frame(pvalt_temp)
              # rownames(pvalt_temp1) = names(pvalt_temp1)
              # colnames(pvalt_temp1) = paste0(test,"_p_none")
              # save(pvalt_temp1, file=cpname)
            }
            pval11 = pvalt_temp[,!grepl(paste(padjust,collapse="|"),colnames(pvalt_temp)),drop=F]
          } else {
            
            loop_ind = loop_ind_f(1:ncol(m),no_cores)
            
            if (test=="chi2") {
              pvalt0 = foreach(ii = loop_ind) %dopar% { 
                pvalii = apply(m[,ii],2,function(x) {
                  ind = !is.na(x)
                  return(ifelse(length(unique(x[ind]))==1, 1, chisq.test(class[ind],x[ind])$p.value)) 
                })
                for (xi in ii) {
                  x = m[,xi]
                  ind = !is.na(x)
                  a = ifelse(length(unique(x[ind]))==1, 1, chisq.test(class[ind],x[ind])$p.value)
                }
                return(pvalii)
              }
              pval11 = data.frame(chi2_p_none=unlist(pvalt0))
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
              png(file=paste0(pname,"_logfc.v.avg"), width=width, height=height*2)
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
} #feat_type

time_output(start)



## plot / enrichment
start = Sys.time()

# enrichment data
data(taSNP) #from dbGaP and NHGRI
data(taSNPLD) #linkage disequilibrium trait-associated SNP
data(CEU) #system("plink -file plink.file -freq -out chr"); to retrieve all SNPs withcorresponding MAF (minor allele frequency) from the CEU vcf files downloaded


pval_paths = gsub(".Rdata","",list.files(gwas_dir, pattern=".Rdata", full.names=T))
feat_types = sapply(str_split(fileNames(pval_paths),"[.]"),function(x) str_split(x[1],"[-]")[[1]][1])

for (pi in 1:length(pval_paths)) {
  start1 = Sys.time()
  
  pval_path = pval_paths[pi]
  feat_type = feat_types[pi]
  
  pval0 = get(load(paste0(pval_path,".Rdata")))
  
  meta_col0 = NULL
  meta_col_path = paste0(meta_col_dir,"-",feat_type,".Rdata")
  mcf = file.exists(meta_col_path)
  if (mcf) {
    meta_col0 = as.data.frame(get(load(meta_col_path)))
    meta_col = meta_col0[match(rownames(pval0),meta_col0[,id_col]),]
  } 
  
  
  ## save as csv
  pv_table = pval0[apply(pval0,1,function(x) any(x<1)),apply(pval0,2,function(x) !all(x==1))]
  if (all(dim(pv_table)>0)) {
    # if (is.null(dim(pv_table))) pv_table = matrix(pv_table,nrow=1)
    if (mcf) pv_table = cbind(meta_col[rownames(pv_table),],pv_table)
    write.csv(pv_table,file=paste0(pval_path,".csv"))
  }
  
  
  
  
  ## PLOT!
  prow = ncol(pval0)
  pcol = 1#length(pval)/prow
  
  # for (col_ind_n in names(col_inds)) {
  #   pval_path = gsub("Xall",paste0("X",col_ind_n),pval_path_)
  #   if (col_ind_n=="all") {
  #     col_ind = names(pval0[[1]])
  #   } else {
  #     col_ind = col_inds[[col_ind_n]]
  #   }
  #   
  #   pval = lapply(pval0, function(x) x[names(x)%in%col_ind])
  # m = m0[file_ind,]
  # m = m[,apply(m,2,function(x) min(table(x))>good_col & length(unique(x))>1)]
  # meta_file = meta_file0[match(rownames(m),meta_file0[,id_col]),]
  
  # qq plot
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
  
  
  if (mcf) {
    # manhattan plot
    
    png(paste0(pval_path,"_manhattan.png"), width=pcol*width, height=prow*height)
    par(mfcol=c(prow,pcol))
    
    for (pvaln in sort(colnames(pval0))) { try({
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
        pos = meta_col$pos_phys
        chr = meta_col$chromosome
      }
      if (grepl("rna",feat_type)) {
        pos = meta_col$start
        chr = meta_col$chr
      }
      
      manhattan_plot(val=-log(pvalt), pos=pos, chrom=chr, 
                     val_thres=-log(pthres), val_max2=-log(1e-4), val_max1=-log(val_max1t),
                     main=paste0(pvaln," sig=", sum(pvalt < pthres)),
                     lines=c(-log(.01),-log(.025),-log(.05),-log(.1)))
      
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
    
    
    png(paste0(pval_path,"_traseR.png"), width=pcol*width, height=prow*height)
    par(mfcol=c(prow,pcol))
    
    for (pvaln in sort(colnames(pval0))) {
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
      
      write.csv(x1$tb1, file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNP.png"))
      write.csv(x2$tb1, file=paste0(pval_path,"_",gsub("_","-",pvaln),"_traseenrichSNPLD.png"))
      
      
      
      # enrichment analysis using hypergeometric test
      sea00 <- sear(input=sigGenes, type = 'mrna')
      sea0 = filter(sea00, collection %in% c("TISSUES")) %>% 
        arrange(fdr)
      sea0$geneset <- as.character(sea0$geneset)
      sea0$geneset[sea0$subcollection == "CNS"] <- "CNS"
      sea0$geneset <- factor(sea0$geneset)
      sea <- sea0 %>% dplyr::group_by(subcollection, geneset) %>% 
        dplyr::summarise(fdr = mean(fdr))
      sea$subcollection <- factor(sea$subcollection, 
                                  levels = toupper(c("cns", "thymus", "tonsils", "spleen", "lymph", "bcells","myeloid", "nkcells","blood","rbcs","neutrophils", "tcells","th1cells","th2cells","tregs","cd8tcells","cd4tcells")))
      subCol <- sea %>% dplyr::group_by(subcollection) %>% 
        dplyr::arrange(subcollection) %>% dplyr::select(subcollection, geneset, fdr) %>% as.data.frame() 
      subCol$geneset <- factor(as.character(subCol$geneset), levels = unique(as.character(subCol$geneset)))
      
      #dev.off(); dev.off();
      cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                     brewer.pal(9, "Set1"))
      pdf(paste0(data_dir, "/results/Figures/Figure2/Figure2c.pdf"), width = 8, height = 5.5)
      ggplot(subCol, aes(x = geneset, y = -log10(fdr), fill = subcollection)) + 
        geom_bar(stat = "identity") +
        customTheme(sizeStripFont=10, xAngle=90, hjust=1, vjust=0.5, 
                    xSize=10, ySize=10, xAxisSize=10, yAxisSize=10) +
        ylab(expression("-log"[10]~"(BH-FDR)")) + xlab("") + scale_fill_manual(values=cbPalette, name = "Tissue")
      dev.off()
      
      
      
      
      
      
      library(InterMineR)
      
      # get HumanMine instance
      im.human = initInterMine(listMines()["HumanMine"])
      
      # retrieve data model for HumanMine
      hsa_model = getModel(im.human)
      
      # all targets of HumanMine enrichment widgets possess InterMine ids
      subset(hsa_model, type %in% c("Gene", "Protein", "SNP") & child_name == "id")
      
      data("PL_DiabetesGenes")
      head(PL_DiabetesGenes, 3)
      
      # get Gene.symbol
      hsa_gene_names = as.character(PL_DiabetesGenes$Gene.symbol)
      head(hsa_gene_names, 3)
      
      # get Gene.primaryIdentifiers (ENTREZ Gene identifiers)
      hsa_gene_entrez = as.character(PL_DiabetesGenes$Gene.primaryIdentifier)
      head(hsa_gene_entrez, 3)
      
      # get all HumanMine widgets
      human.widgets = as.data.frame(getWidgets(im.human))
      human.widgets
      
      # get enrichment, gene-related widgets for HumanMine
      subset(human.widgets, widgetType == 'enrichment' & targets == "SNP")
      
      # Perform enrichment analysis
      # returns: a data.frame with the results of the enrichment analysis which was performed in the defined InterMine platform. p-values are given after applying the correction algorithm. Count and populationAnnotationCount columns contain the number of genes that belong to each GO term in the list and in the background population respectively.

      GO_enrichResult = doEnrichment(
        im = im.human,
        ids = hsa_gene_entrez,
        widget = "go_enrichment_for_gene"
        # organism = "Homo sapiens" # optional if results from more than one organism are retrieved
      )
      
      head(GO_enrichResult$data)
      dim(GO_enrichResult$data)
      GO_enrichResult$populationCount # size of the reference background population (populationCount)
      GO_enrichResult$notAnalysed #number of input features that were not included in the enrichment analysis (notAnalyzed)
      GO_enrichResult$im # name and url of the Mine (im)
      
      
      
      
      ## get available filter values for Gene Ontology Enrichment widget
      as.character(subset(human.widgets, name == "go_enrichment_for_gene")$filters)
      
      # Perform enrichment analysis for GO molecular function terms
      GO_MF_enrichResult = doEnrichment(
        im = im.human,
        ids = hsa_gene_entrez,
        widget = "go_enrichment_for_gene",
        filter = "molecular_function")
      
      head(GO_MF_enrichResult$data)
      
      ## multiple testing adjustment
      
      # Perform enrichment analysis for Protein Domains in list of genes
      PD_FDR_enrichResult = doEnrichment(
        im = im.human,
        ids = hsa_gene_entrez,
        widget = "prot_dom_enrichment_for_gene",
        correction = "Benjamini Hochberg"
      ) 
      
      head(PD_FDR_enrichResult$data)
      
      # Perform enrichment analysis for Protein Domains in list of genes
      # but without a correction algoritm this time
      PD_None_enrichResult = doEnrichment(
        im = im.human,
        ids = hsa_gene_entrez,
        widget = "prot_dom_enrichment_for_gene",
        correction = "None"
      )
      
      head(PD_None_enrichResult$data)
      
      ## visualize
      library(GeneAnswers)
      # convert InterMineR Gene Ontology Enrichment analysis results to GeneAnswers object
      geneanswer_object = convertToGeneAnswers(
        
        # assign with doEnrichment result:
        enrichmentResult = GO_enrichResult,
        
        # assign with list of genes:
        geneInput = data.frame(GeneID = as.character(hsa_gene_entrez), 
                               stringsAsFactors = FALSE),
        
        # assign with the type of gene identifier
        # in our example we use Gene.primaryIdentifier values (ENTREZ IDs):
        geneInputType = "Gene.primaryIdentifier",
        
        # assign with Bioconductor annotation package:
        annLib = 'org.Hs.eg.db',
        
        # assign with annotation category type
        # in our example we use Gene Ontology (GO) terms:
        categoryType = "GO"
        
        #optional define manually if 'enrichIdentifier' is missing from getWidgets:
        #enrichCategoryChildName = "Gene.goAnnotation.ontologyTerm.parents.identifier"
      )
      
      class(geneanswer_object)
      summary(geneanswer_object)
      
      # convert to GeneAnswers results for GO terms associated with molecular function
      geneanswer_MF_object = convertToGeneAnswers(
        enrichmentResult = GO_MF_enrichResult,
        geneInput = data.frame(GeneID = as.character(hsa_gene_entrez), 
                               stringsAsFactors = FALSE),
        geneInputType = "Gene.primaryIdentifier",
        annLib = 'org.Hs.eg.db',
        categoryType = "GO.MF"
        #enrichCategoryChildName = "Gene.goAnnotation.ontologyTerm.parents.identifier"
      )
      
      class(geneanswer_MF_object)
      
      
      
      
      
      
      # GeneAnswers barplot
      geneAnswersChartPlots(geneanswer_object, 
                            chartType='barPlot',
                            sortBy = 'correctedPvalue',
                            top = 5)
      # generate concept-gene network
      geneAnswersConceptNet(geneanswer_object, 
                            colorValueColumn=NULL,
                            centroidSize='correctedPvalue', 
                            output='interactive',
                            catTerm = FALSE,
                            catID = FALSE,
                            showCats = 1:5)
    }
  }
  
  
  
  
  
  time_output(start1)
} #pi
time_output(start)












