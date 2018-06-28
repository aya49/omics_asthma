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
gwas_dir = paste0(stat_dir,"/gwas_source"); dir.create(gwas_dir,showWarnings=F)
qq_dir = paste0(stat_dir,"/gwas_plot-qq"); dir.create(qq_dir,showWarnings=F)
mh_dir = paste0(stat_dir,"/gwas_plot-manhattan"); dir.create(mh_dir,showWarnings=F)

## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("data.table")
libr("qqman")
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

overwrite = T
writecsv = T #write results as csv on top of Rdata?

pthres = .05
padjust = p.adjust.methods
val_max1t = 5e-8 #comment out if just want default, large large p value threshold guideline value

tests = c("chi2")

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis
cont_col = 15 #if a column has more than this unique values, it's continuous
scale_cont = T #if a column is continuous, scale it

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




start = Sys.time()


meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
# meta_col0 = as.data.frame(get(load(paste0(meta_col_dir,".Rdata"))))

for (feat_type in feat_types) {
  cat("\n", feat_type)
  
  m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
  if (sum(colnames(m0)%in%meta_file0[,id_col])==ncol(m0)) m0 = t(m0)
  
  col_inds_i = sapply(names(col_inds0), function(x) grepl(x,feat_type))
  if (sum(col_inds_i)==0) {
    col_inds = list(all=c(""))
  } else {
    col_inds = col_inds0[[col_inds_i]]
  }
  
  
  for (class_coli in 1:length(class_cols)) {
    class_col = class_cols[class_coli]
    control = controls[class_coli]
    
    for (file_ind_n in names(file_inds)) {
      file_ind = file_inds[[file_ind_n]]
      if (file_ind_n=="all") file_ind = rownames(m0)
      # file_ind_n = paste0("-",file_ind_n)
      
      # ensure there is a class label for each sample
      file_ind = file_ind[!is.na(meta_file0[match(file_ind,meta_file0[,id_col]),class_col])]
      
      for (col_ind_n in names(col_inds)) {
        col_ind = col_inds[[col_ind_n]]
        if (col_ind_n=="all") col_ind = colnames(m0)
        # col_ind_n = paste0(".",col_ind_n)
        
        m = m0[rownames(m0)%in%file_ind, colnames(m0)%in%col_ind]
        m = delna(m)
        if (any(dim(m)==0) | sum(!is.na(m))<sum(is.na(m))) next()
        m_cont_col = apply(m,2,function(x) length(unique(x))>cont_col)
        if (scale_cont) m[,m_cont_col] = scale(as.numeric(as.matrix(m[,m_cont_col])))
        m[is.na(m)] = -1
        
        cat(" (",file_ind_n, " x ",col_ind_n,"; ",nrow(m)," x ",ncol(m),") ")
        
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
          
          pv_table = Reduce("cbind",pval)
          colnames(pv_table) = names(pval)
          pv_table = pv_table[apply(pv_table,1,function(x) any(x<1)),apply(pv_table,2,function(x) any(x<1))]
          
          save(pv_table,file=paste0(pname,".Rdata"))
        } #test
      } #col_ind
    } #row_ind
  } #class_col
} #feat_type



## plot

pval_paths = gsub(".Rdata","",list.files(gwas_dir, pattern=".Rdata", full.names=T))
feat_types = sapply(str_split(fileNames(pval_paths),"[.]"),function(x) x[1])

for (pi in pval_paths) {
  pval_path = pval_paths[pi]
  feat_type = feat_types[pi]
  
  pval0 = get(load(paste0(pval_path,".Rdata")))
  
  meta_col0 = NULL
  meta_col_path = paste0(meta_col_dir,"-",feat_type,".Rdata")
  mcf = file.exists(meta_col_path)
  if (mcf) {
    meta_col0 = as.data.frame(get(load(meta_col_path)))
    meta_col = meta_col0[match(names(pval[[1]]),meta_col0[,cid_col]),]
  } 
  
  
  ## save as csv
  pv_table = pval0[apply(pval0,1,function(x) any(x<1)),apply(pval0,2,function(x) any(x<1))]
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
  png(paste0(gsub(gwas_dir,qq_dir,pval_path),".png"), width=pcol*width, height=prow*height)
  par(mfcol=c(prow,pcol))
  for (pvaln in sort(colnames(pval))) { try({
    pvalt = pval[,pvaln]
    # if (sum(pvalt>val_max1t)>5) val_max1t = 5e-6
    pv_ind = pvalt<1
    if (sum(pv_ind)==0) next()
    
    qq(pvalt)
  }) }
  graphics.off()
  
  
  # manhattan plot
  if (mcf) {
    png(paste0(gsub(gwas_dir,mh_dir,pval_path),".png"), width=pcol*width, height=prow*height)
    par(mfcol=c(prow,pcol))
    
    for (pvaln in sort(colnames(pval))) { try({
      pvalt = pval[,pvaln]
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
  
} #pi
time_output(start)












