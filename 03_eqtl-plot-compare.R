## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result")


# asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
eqtl_dir = paste0(stat_dir,"/eqtl"); dir.create(eqtl_dir, showWarnings=F)
gwas_dir = paste0(stat_dir,"/gwas"); dir.create(gwas_dir,showWarnings=F)


## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("data.table")
libr("MatrixEQTL")
libr("foreach")
libr("doMC")
libr("stringr")
libr("qdapTools") # make dummy variables
libr("Matrix")
libr("ggplot2")
libr("plyr")



## options
no_cores = 15#detectCores()-3
registerDoMC(no_cores)

overwrite = T
writecsv = T #write results as csv on top of Rdata?

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis
good_na = .75

id_col = "id"
class_col = "response"
categorical = T # is class column categorical?
interested_cols = c("response")
# interested_cols = c("age","bmi","sex","centre","batch","race","response") 
# interested_cont_cols = ""

useModels = c("modelLINEAR", "modelANOVA", "modelLINEAR_CROSS")

pthres = .01
max_plots = 150 # max number of plots to make

# plotting size
width = 500
height = 300

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




## plot eQTL ---------------------------------------------

start = Sys.time()

pval_paths = list.files(gwas_dir, full.names=T, pattern=".Rdata")
eqtl_paths = list.files(eqtl_dir, full.names=T, pattern=".Rdata")
eqtl_names = fileNames(eqtl_paths)
feat_types = lapply(str_split(eqtl_names,"_"), function(x) c(str_split(x,"[-]")[[1]][1], str_split(x,"[-]")[[1]][2]))

meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
meta_file0[,id_col] = as.character(meta_file0[,id_col])

for (ei in 1:length(eqtl_paths)) {
  try({
    # a = foreach (ei = 1:length(eqtl_paths)) %dopar% {
    eqtl_path = eqtl_paths[ei]
    eqtl_name = eqtl_names[ei]
    feat_type = feat_types[[ei]]
    bins = ""
    if (grepl("genotype",feat_type[1])) {
      if (grepl("01", feat_type[1])) bins = "01"
      if (grepl("12", feat_type[1])) bins = "12"
      feat_type[1] = gsub("[.]01|[.]12","",feat_type[1])
    }
    
    ## compare pre/diff/post
    if (!(grepl("genotype",feat_type[1]) & grepl(".post",feat_type[2]) & 
          file.exists(gsub(".post",".pre", eqtl_path)) & 
          file.exists(gsub(".post",".diff", eqtl_path)))) next()
    
    
    # load eqtl / matrix
    mecisl = NULL
    m1 = get(load(paste0(feat_dir,"/",feat_type[1],".Rdata")))
    m2l = NULL
    for (met in c("pre","post","diff")) {
      try({ mecis = get(load(gsub(".post",paste0(".",met), eqtl_path))) })
      mecis$snps = as.character(mecis$snps)
      mecis$gene = as.character(mecis$gene)
      mecisl[[met]] = mecis
      
      m2 = get(load(paste0(feat_dir,"/",gsub(".post",paste0(".",met), feat_type[2]),".Rdata")))
      m2l[[met]] = m2
    }
    row_ind = sort(intersect(rownames(m1), Reduce("intersect", lapply(m2l,rownames))))
    mecisl_names = Reduce("union", lapply(mecisl, function(x) paste(x$snps, x$gene, sep="_")))
    mecisl = lapply(mecisl, function(x) {
      x = x[match(mecisl_names,paste(x$snps, x$gene, sep="_")),]
      x[,c("snps","gene")] = Reduce("rbind",strsplit(mecisl_names,"_"))
      return(x)
    })
    mecis = mecisl[[1]]
    
    # prepare meta_file
    file_ind_n = gsub("[-]|X","",str_extract(eqtl_name,"-[a-zA-Z]+X"))
    file_ind = file_inds[[file_ind_n]]
    if (file_ind_n=="all") file_ind = meta_file0[,id_col]
    meta_file = meta_file0[meta_file0[,id_col]%in%file_ind & meta_file0[,id_col]%in%row_ind,]
    class = as.numeric(factor(meta_file[,class_col]))
    class_name = as.character(levels(factor(meta_file[,class_col])))
    class_n = as.character(meta_file[,class_col])
    
    # trim matrix rows
    m1 = m1[meta_file[,id_col], ]
    m2l = lapply(m2l, function(m2) m2[meta_file[,id_col], ])
    
    # load meta_col & trim
    meta_col10 = meta_col1 = NULL
    meta_col_path1 = paste0(meta_col_dir,"-",feat_type[1],".Rdata")
    mcf1 = file.exists(meta_col_path1)
    if (mcf1) {
      meta_col10 = as.data.frame(get(load(meta_col_path1)), stringsAsFactors=F)
      if (grepl("genotype", feat_type[1])) {
        meta_col1 = meta_col10[match(mecis$snps,meta_col10[,"dbSNP"]),]
        mecisl = lapply(mecisl, function(mecis) {
          mecis = cbind(mecis,meta_col1[,c("chromosome","pos_phys")])
          colnames(mecis)[colnames(mecis)%in%c("chromosome","pos_phys")] = paste0(c("chromosome","pos_phys"),".",feat_type[1])
          return(mecis)
        })
      } else {
        meta_col1 = meta_col10[match(mecis$snps,meta_col10[,id_col]),]
      }
      m1 = m1[,meta_col1[,id_col]]
    } else {
      m1 = m1[,mecis$snps]
    }
    meta_col20 = meta_col2 = NULL
    meta_col_path2 = paste0(meta_col_dir,"-",strsplit(feat_type[2],"[.]")[[1]][1],".Rdata")
    mcf2 = file.exists(meta_col_path2)
    if (mcf2) {
      meta_col20 = as.data.frame(get(load(meta_col_path2)))
      if (grepl("rna", feat_type[2])) {
        meta_col2 = meta_col20[match(mecis$gene,meta_col20[,"symbol"]),]
        mecisl = lapply(mecisl, function(mecis) {
          mecis = cbind(mecis,meta_col2[,c("chr","start","end")])
          colnames(mecis)[colnames(mecis)%in%c("chr","start","end")] = paste0(c("chr","start","end"),".",feat_type[1])
          return(mecis)
        })
      } else {
        meta_col2 = meta_col20[match(mecis$gene,meta_col20[,id_col]),]
      }
      m2l = lapply(m2l, function(m2) m2[,meta_col2[,id_col]])
    } else {
      m2l = lapply(m2l, function(m2) m2[,mecis$gene])
    }
    
    # load gwas
    pval1 = get(load(pval_paths[grepl(paste0(feat_type[1],"-",file_ind_n,"Xall"),pval_paths)]))
    pvalt1 = pval1[,grepl("none",colnames(pval1))]
    pvalt1 = pvalt1[match(colnames(m1),names(pvalt1))]
    
    pvalt2l = NULL
    for (met in c("pre","post","diff")) {
      pval2 = get(load(pval_paths[grepl(paste0(gsub(".post",paste0(".",met), feat_type[2]),"-",file_ind_n,"Xall"),pval_paths)]))
      pvalt2 = pval2[,grepl("none",colnames(pval2))]
      pvalt2 = pvalt2[match(colnames(m2),names(pvalt2))]
      pvalt2l[[met]] = pvalt2
    }
    
    # add p values to eqtl cis table
    
    
    # adjust m1 if needed
    if (bins=="01") m1[m1==2] = 1
    if (bins=="12") m1[m1==0] = 1
    
    dir.create(gsub(".Rdata|.post","",eqtl_path), showWarnings=F)
    dfm =  rbind.fill(lapply(names(m2l), function(x) data.frame(response=meta_file[,class_col], time=rep(x,nrow(m1)))))
    df0 = lapply(1:ncol(m1), function(i) 
      rbind.fill(lapply(names(m2l), function(x) 
        data.frame(genotype=m1[,i], cont=m2l[[x]][,i])))
    )
    fdrs = Reduce("cbind",lapply(mecisl, function(x) x$FDR))
    plot_inds = which(apply(fdrs, 1, function(y) any(y[!is.na(y)]<pthres)))
    if (any(grepl("beta",colnames(mecisl[[1]]))))
      plot_inds = plot_inds[order(mecisl$diff$beta,decreasing=T)]
    for (i in plot_inds) {
      main_ = paste0("eqtl ", paste(names(mecisl),collapse="/"),
                     ": stat=",paste(sapply(mecisl, function(mecis) round(mecis$statistic[i],4)),collapse="/"), 
                     "\nunadj.p=", paste(sapply(mecisl, function(mecis) round(mecis$pvalue[i],4)),collapse="/"), 
                     "\nfdr.p=", paste(sapply(mecisl, function(mecis) round(mecis$FDR[i],4)),collapse="/"))
      
      xlab_ = paste0(mecisl[[1]]$snps[i]," unadj.p = ", round(pvalt1[i],4))
      
      ylab_ = paste0(mecisl[[1]]$gene[i]," unadj.p = ", paste(names(pvalt2l), collapse="/")," = ", paste(sapply(pvalt2l, function(pvalt2) round(pvalt2[i],4)), collapse="/"))
      
      dfi = cbind(dfm,df0[[i]])
      
      png(file=paste0(gsub(".Rdata|.post","",eqtl_path),"/",mecis$snps[i],"-",mecis$gene[i],".png"), width=width, height=height*length(unique(meta_file[,class_col])))
      
      p = ggplot(dfi, aes(interaction(response,genotype), cont)) +
        geom_boxplot(aes(x=interaction(response,genotype), y=cont, fill=response), alpha = 0.5) +
        geom_jitter(aes(colour=response), width=0.2) +
        geom_smooth(aes(group=response), method="lm", size = 2, se = F) +
        facet_grid(time~.,scales="free_x")
      
      print(p + xlab(xlab_) + ylab(ylab_) + ggtitle(main_))
      graphics.off()
    }
  })
}
time_output(start)










