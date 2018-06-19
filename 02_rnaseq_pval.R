## input: features & meta_file
## output: p vaue features testing significance of ER DR correlation
## aya43@sfu.ca
## created 20180524



## root directory
root = "~/projects/asthma"
setwd(root)

type = "genes" #"isoforms", "genes"
result_dir = paste0(root, "/result/RNAseq_",type); dir.create(result_dir, showWarnings=F)


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")#,asthma)

feat_dir = paste0(result_dir,"/feat")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir,showWarnings=F)
gwas_dir = paste0(result_dir,"/gwas"); dir.create(gwas_dir,showWarnings=F)

## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("data.table")
libr("MatrixEQTL")
libr("entropy")
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

pthres = .05
padjust = p.adjust.methods

tests = c("chi2")

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis

id_col = "fileName"
class_col = "response"
comparisons = list(c("time-pre","all","response"),
                   c("time-post","all","response"),
                   c("time-pre","all","response"),
                   c("time-post","all","response")) #list of (row, col, class)
control = "ER"
categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""

cid_col = "probe"

#plot size
width = 1000
height = 500



## features and indices
feat_types = list.files(feat_dir,full.names=F,pattern=".Rdata")
feat_types = gsub(".Rdata","",feat_types)

col_inds_paths = list.files(meta_dir,pattern="col_id_",full.names=T)
col_inds_names = sapply(str_split(gsub(".Rdata","",col_inds_paths),"_"), function(x) x[length(x)])
col_inds = lapply(col_inds_paths, function(x) get(load(x)))
names(col_inds) = col_inds_names
col_inds = append(list(all=""), col_inds)

file_inds_paths = list.files(meta_dir,pattern="file_id_",full.names=T)
file_inds_names = sapply(str_split(gsub(".Rdata","",file_inds_paths),"_"), function(x) x[length(x)])
file_inds = lapply(file_inds_paths, function(x) get(load(x)))
names(file_inds) = file_inds_names
file_inds = append(list(all=""), file_inds)




start = Sys.time()


meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
meta_col0 = as.data.frame(get(load(paste0(meta_col_dir,".Rdata"))))

for (feat_type in feat_types) {
  m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
  if (sum(colnames(m0)%in%meta_file0[,id_col])==ncol(m0)) m0 = t(m0)
  m0[,]
  
  for (comparison in comparisons) {
    
  }
  
  for (file_ind_n in names(file_inds)) {
    file_ind = file_inds[[file_ind_n]]
    if (file_ind_n=="all") file_ind = rownames(m0)
    # file_ind_n = paste0("-",file_ind_n)
    
    for (col_ind_n in names(col_inds)) {
      col_ind = col_inds[[col_ind_n]]
      if (col_ind_n=="all") col_ind = colnames(m0)
      # col_ind_n = paste0(".",col_ind_n)
      
      m = m0[file_ind[file_ind%in%rownames(m0)], col_ind[col_ind%in%colnames(m0)]]
      
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
        
        ## p value calculation
        if (!overwrite & file.exists(paste0(pname,".Rdata"))) next()
        
        cpname = paste0(gsub(paste0("X",col_ind_n),"Xall",pname),".Rdata")
        if (col_ind_n!="all" & file.exists(cpname)) {
          pvalt_temp = get(load(cpname))
          pvalt = pvalt_temp[[grep("_none$",names(pvalt_temp))]]
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
            pvalt = unlist(pvalt0)
          }
          pvalt = unlist(pvalt)
          names(pvalt) = colnames(m)
        }
        
        
        # pval1$chi2 = apply(m,2,function(x) { 
        #   ind = !is.na(x)
        #   return(ifelse(length(unique(x[ind]))==1, 1, chisq.test(class[ind],x[ind])$p.value)) 
        # })
        
        # pval1$mi = apply(m,2,function(x) {
        #   ind = !is.na(x)
        #   return(ifelse(length(unique(x[ind]))==1, 1, mcnemar.test(factor(class[ind]),factor(x[ind]))))
        # })
        
        ## p value adjust
        pval = foreach(padj = padjust) %dopar% {
          # pval[[paste0(pvaln,"_",padj)]] = p.adjust(pvalt,method=padj)
          return(p.adjust(pvalt,method=padj))
        }
        names(pval) = paste0(test,"_",padjust)
        
        save(pval,file=paste0(pname,".Rdata"))
      } #test
    } #col_ind
  } #row_ind
} #feat_type



## plot

pval_paths = gsub(".Rdata","",list.files(gwas_dir, pattern=".Rdata", full.names=T))

for (pval_path in pval_paths) {
  
  pval = get(load(paste0(pval_path,".Rdata")))
  meta_col = meta_col0[match(names(pval[[1]]),meta_col0[,cid_col]),]
  
  
  prow = length(pval)
  pcol = 2#length(pval)/prow
  png(paste0(pval_path,".png"), width=pcol*width, height=prow*height)
  par(mfcol=c(prow,pcol))
  
  for (pvaln in sort(names(pval))) {
    try({
      pvalt = unlist(pval[[pvaln]])
      val_max1t = 1e-4
      if (sum(pvalt>val_max1t)>5) val_max1t = 5e-6
      if (all(pvalt==1)) next()
      
      # chrom = meta_col$chromosome
      # chrom[chrom=="X"]="23"
      # chrom[chrom=="Y"]="24"
      # chrom[chrom=="XY"]="25"
      # chrom[chrom=="MT"]="26"
      # 
      # manhattan(data.frame(P=pvalt, BP=meta_col$pos_phys, CHR=as.numeric(chrom), SNP=meta_col$dbSNP))
      qq(pvalt)
      manhattan_plot(val=-log(pvalt), pos=meta_col$pos_phys, chrom=meta_col$chromosome, 
                     val_thres=-log(pthres), val_max2=-log(1e-4), val_max1=-log(val_max1t),
                     main=paste0(pvaln," sig=", sum(pvalt < pthres)),
                     lines=c(-log(.01),-log(.025),-log(.05),-log(.1)))
      
    })
  }
  # require(gridExtra)
  # grid.arrange(plots,ncol=pcol,nrow=prow)
  # plotsa = marrangeGrob(grob=plots,ncol=pcol,nrow=prow)
  # ggexport(plotsa, filename = paste0(gwas_dir,"/snp-file-genotype_norm-NA.png"))
  # ggsave(file=paste0(gwas_dir,"/snp-file-genotype_norm-NA.png"), plotsa)
  graphics.off()
  
}



time_output(start)












