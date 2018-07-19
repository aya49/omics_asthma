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
    
    # load eqtl
    mecis = NULL
    try({ mecis = get(load(eqtl_path)) })
    if (is.null(mecis)) next()
    # mecis = me$cis$eqtls
    mecis$snps = as.character(mecis$snps)
    mecis$gene = as.character(mecis$gene)
    
    # load matrix
    m1 = get(load(paste0(feat_dir,"/",feat_type[1],".Rdata")))
    m2 = get(load(paste0(feat_dir,"/",feat_type[2],".Rdata")))
    row_ind = sort(intersect(rownames(m1), rownames(m2)))
    
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
    m2 = m2[meta_file[,id_col], ]
    
    # load meta_col & trim
    meta_col10 = meta_col1 = NULL
    meta_col_path1 = paste0(meta_col_dir,"-",feat_type[1],".Rdata")
    mcf1 = file.exists(meta_col_path1)
    if (mcf1) {
      meta_col10 = as.data.frame(get(load(meta_col_path1)), stringsAsFactors=F)
      if (grepl("genotype", feat_type[1])) {
        meta_col1 = meta_col10[match(mecis$snps,meta_col10[,"dbSNP"]),]
        mecis = cbind(mecis,meta_col1[,c("chromosome","pos_phys")])
        colnames(mecis)[colnames(mecis)%in%c("chromosome","pos_phys")] = paste0(c("chromosome","pos_phys"),".",feat_type[1])
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
        mecis = cbind(mecis,meta_col2[,c("chr","start","end")])
        colnames(mecis)[colnames(mecis)%in%c("chr","start","end")] = paste0(c("chr","start","end"),".",feat_type[2])
      } else {
        meta_col2 = meta_col20[match(mecis$gene,meta_col20[,id_col]),]
      }
      m2 = m2[,meta_col2[,id_col]]
    } else {
      m2 = m2[,mecis$gene]
    }
    
    # load gwas
    pval1 = get(load(pval_paths[grepl(paste0(feat_type[1],"-",file_ind_n,"Xall"),pval_paths)]))
    pvalt1 = pval1[,grepl("none",colnames(pval1))]
    pvalt1 = pvalt1[match(colnames(m1),names(pvalt1))]
    pval2 = get(load(pval_paths[grepl(paste0(feat_type[2],"-",file_ind_n,"Xall"),pval_paths)]))
    pvalt2 = pval2[,grepl("none",colnames(pval2))]
    pvalt2 = pvalt2[match(colnames(m2),names(pvalt2))]
    
    # add p values to eqtl cis table
    
    
    # adjust m1 if needed
    if (bins=="01") m1[m1==2] = 1
    if (bins=="12") m1[m1==0] = 1
    
    # save as table
    write.csv(mecis, file=gsub(".Rdata",".csv",eqtl_path))
    
    # plot the histogram of local and distant p-values
    # png(gsub(".Rdata","_qq.png",eqtl_path), width=width, height=height)
    # plot(me)
    # graphics.off()
    
    dir.create(gsub(".Rdata","",eqtl_path), showWarnings=F)
    for (i in 1:min(nrow(mecis),max_plots)) {
      main_ = paste0("eqtl: stat=",round(mecis$statistic[i],4), "; unadj.p=", round(mecis$pvalue[i],4), "; fdr.p=", round(mecis$FDR[i],4))
      
      xlab_more = ""
      if(grepl("genotype",feat_type[1])) {
        f1_colind = grepl(feat_type[1],colnames(mecis))
        xlab_more = paste0("\n", paste(gsub(paste0(".",feat_type[1]), "",colnames(mecis)[f1_colind]), collapse="/"), 
                           ": ", paste(mecis[i,f1_colind], collapse=" / "))
      } 
      xlab_ = paste0(feat_type[1],": ", mecis$snps[i], "; unadj.p = ", round(pvalt1[i],4), xlab_more)
      
      ylab_more = ""
      if(grepl("rna",feat_type[2])) {
        f2_colind = grepl(feat_type[2],colnames(mecis))
        ylab_more = paste0("\n", paste(gsub(paste0(".",feat_type[2]),"",colnames(mecis)[f2_colind]), collapse="/"), ": ", 
                           paste(mecis[i,f2_colind], collapse=" / "))
      } 
      ylab_ = paste0(feat_type[2],": ", mecis$gene[i], "; unadj.p = ", round(pvalt2[i],4), ylab_more)
      
      png(file=paste0(gsub(".Rdata","",eqtl_path),"/",mecis$snps[i],"-",meta_col1$ mecis$gene[i],".png"), width=width, height=height)
      if (feat_type[1]=="genotype" ) {
        p = ggplot(data.frame(a=factor(m1[,i]),b=m2[,i]), aes(interaction(class_n,a), b)) + 
          geom_boxplot(aes(colour=class_n)) + 
          geom_smooth(aes(group=class_n), method="lm", size = 2, se = F) +
          stat_summary(fun.y=mean, geom="line") + 
          geom_jitter(width=0.2,aes(colour=class_n))
        
        # ggplot(data.frame(class_n=class_n, m1i=m1[,i], m2i=m2[,i]), 
        #        aes(x=interaction(class_n, m1i), y=m2i)) +
        #   geom_boxplot(aes(fill = class_n), alpha = 0.5) +
        #   geom_line(aes(group = m1i),
        #             alpha = 0.5, colour = "darkgrey")
        # facet_grid(sp~variable,scales="free_x")
      } else {
        p = ggplot(data.frame(a=m1[,i],b=m2[,i]), aes(interaction(class_n,a), b)) + 
          geom_point(aes(colour=class_n)) +
          geom_smooth(method=lm)
      }
      print(p + xlab(xlab_) + ylab(ylab_) + ggtitle(main_))
      graphics.off()
    }
    
    
  })
}
time_output(start)










