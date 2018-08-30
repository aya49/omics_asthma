## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614


## logistics
root = "~/projects/asthma"; commandArgs <- function(...) root  # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/visualizationFunctions.R"))
libr(append(pkgs(),c("qdapTools")))

no_cores = 8#detectCores()-3
registerDoMC(no_cores)


## options
overwrite = T
writecsv = T #write results as csv on top of Rdata?

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



## plot eQTL ---------------------------------------------

start = Sys.time()

pval_paths = list.files(gwas_dir, full.names=T, pattern=".Rdata")
eqtl_paths = list.files(eqtl_dir, full.names=T, pattern=".Rdata")
eqtl_names = fileNames(eqtl_paths)
feat_type_sets = lapply(str_split(eqtl_names,"_"), function(x) c(str_split(x,"[-]")[[1]][1], str_split(x,"[-]")[[1]][2]))
feat_names = lapply(feat_type_sets, function(x) sapply(str_split(x,"[.]"), function(y) y[1]))


a = foreach (ei = 1:length(eqtl_paths)) %dopar% { try({
    # a = foreach (ei = 1:length(eqtl_paths)) %dopar% {
  # a = foreach (ei = 1:length(eqtl_paths)) %dopar% {
  eqtl_path = eqtl_paths[ei]
  eqtl_name = eqtl_names[ei]
  feat_type_set = feat_type_set_ = feat_type_sets[[ei]]
  bins = ""
  if (grepl("dna",feat_type_set[1])) {
    if (grepl("01", feat_type_set[1])) bins = "01"
    if (grepl("12", feat_type_set[1])) bins = "12"
    feat_type_set[1] = gsub("[.]01|[.]12","",feat_type_set[1])
  }
  feat_name = feat_names[[ei]]
  
  file_ind = unlist(read.csv(gsub(".Rdata","_id.csv",eqtl_path)))
  file_ind_n = gsub("[-]|X","",str_extract(eqtl_name,"-[a-zA-Z]+X"))
  # file_ind = file_inds[[file_ind_n]]
  # if (file_ind_n=="all") file_ind = meta_file0[,id_col]
  # meta_file = meta_file0[meta_file0[,id_col]%in%file_ind & meta_file0[,id_col]%in%row_ind,]
  # if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiments[class_coli]
  
    ## compare pre/diff/post
    if (!(grepl("dna",feat_type[1]) & grepl(".post",feat_type[2]) & 
          file.exists(gsub(".post",".pre", eqtl_path)) & 
          file.exists(gsub(".post",".diff", eqtl_path)))) next()
    
    
    # load eqtl / matrix
    mecisl = NULL
    m1 = get(load(paste0(feat_dir,"/",feat_type_set[1],".Rdata")))
    m2l = NULL
    for (met in c("pre","post","diff")) { try({ 
      mecis = get(load(gsub(paste0(feat_type_set[2],".post"),paste0(feat_type_set[2],".",met), eqtl_path))) 
      mecis$snps = as.character(mecis$snps)
      mecis$gene = as.character(mecis$gene)
      mecisl[[met]] = mecis
      
      m2 = get(load(paste0(feat_dir,"/",gsub(".post",paste0(".",met), feat_type_set[2]),".Rdata")))
      m2l[[met]] = m2
    }) }
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
    if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiments
    
    class = as.numeric(factor(meta_file[,class_col]))
    class_name = as.character(levels(factor(meta_file[,class_col])))
    class_n = as.character(meta_file[,class_col])
    
    # trim matrix rows
    m1 = m1[meta_file[,id_col], ]
    m2l = lapply(m2l, function(m2) m2[meta_file[,id_col], ])
    
    # load meta_col & trim
    meta_col10 = meta_col1 = NULL
    meta_col_path1 = paste0(meta_col_dir,feat_name[1],".Rdata")
    mcf1 = file.exists(meta_col_path1)
    if (mcf1) {
      meta_col10 = as.data.frame(get(load(meta_col_path1)), stringsAsFactors=F)
      if (grepl("dna", feat_type_set[1])) {
        meta_col1 = meta_col10[match(mecis$snps,meta_col10[,"dbSNP"]),]
        mecisl = lapply(mecisl, function(mecis) {
          mecis = cbind(mecis,meta_col1[,c("chromosome","pos_phys")])
          colnames(mecis)[colnames(mecis)%in%c("chromosome","pos_phys")] = paste0(c("chromosome","pos_phys"),".",feat_type_set[1])
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
    meta_col_path2 = paste0(meta_col_dir,feat_name[2],".Rdata")
    mcf2 = file.exists(meta_col_path2)
    if (mcf2) {
      meta_col20 = as.data.frame(get(load(meta_col_path2)))
      if (grepl("rna", feat_type_set[2])) {
        meta_col2 = meta_col20[match(mecis$gene,meta_col20[,"symbol"]),]
        mecisl = lapply(mecisl, function(mecis) {
          mecis = cbind(mecis,meta_col2[,c("chr","start","end")])
          colnames(mecis)[colnames(mecis)%in%c("chr","start","end")] = paste0(c("chr","start","end"),".",feat_type_set[1])
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
    pval1 = get(load(pval_paths[grepl(paste0(feat_type_set[1],"-",file_ind_n,"Xall"),pval_paths)]))
    pvalt1 = pval1[,grepl("none",colnames(pval1))]
    pvalt1 = pvalt1[match(colnames(m1),names(pvalt1))]
    
    pvalt2l = NULL
    for (met in c("pre","post","diff")) {
      pval2 = get(load(pval_paths[grepl(paste0(gsub(".post",paste0(".",met), feat_type_set[2]),"-",file_ind_n,"Xall"),pval_paths)]))
      pvalt2 = pval2[,grepl("none",colnames(pval2))]
      pvalt2 = pvalt2[match(colnames(m2),names(pvalt2))]
      pvalt2l[[met]] = pvalt2
    }
    
    # add p values to eqtl cis table
    
    
    # adjust m1 if needed
    if (bins=="01") m1[m1==2] = 1
    if (bins=="12") m1[m1==0] = 1
    
    dfm =  rbind.fill(lapply(names(m2l), function(x) data.frame(response=meta_file[,class_col], time=rep(x,nrow(m1)))))
    df0 = lapply(1:ncol(m1), function(i) 
      rbind.fill(lapply(names(m2l), function(x) 
        data.frame(dna=m1[,i], cont=m2l[[x]][,i])))
    )
    fdrs = Reduce("cbind",lapply(mecisl, function(x) x$FDR))
    plot_inds = which(apply(fdrs, 1, function(y) any(y[!is.na(y)]<pthres)))
    if (length(plot_inds)==0) plot_inds = 1:nrow(mecisl[[1]])
    if (any(grepl("beta",colnames(mecisl[[1]]))))
      plot_inds = plot_inds[order(mecisl$diff$beta[plot_inds],decreasing=T)]
    dir.create(gsub(".Rdata|.post","",eqtl_path), showWarnings=F)
    for (i in plot_inds) {
      pname = paste0(gsub(".Rdata|.post","",eqtl_path),"/",mecisl$diff$snps[i],"-",mecisl$diff$gene[i],".png")
      if(!overwrite & file.exists(pname)) next()
      png(file=pname, width=width, height=height*length(unique(meta_file[,class_col])))
      
      main_ = paste0("eqtl ", paste(names(mecisl),collapse="/"),
                     ": stat=",paste(sapply(mecisl, function(mecis) round(mecis$statistic[i],4)),collapse="/"), 
                     "\nunadj.p=", paste(sapply(mecisl, function(mecis) round(mecis$pvalue[i],4)),collapse="/"), 
                     "\nfdr.p=", paste(sapply(mecisl, function(mecis) round(mecis$FDR[i],4)),collapse="/"))
      
      xlab_ = paste0(mecisl[[1]]$snps[i]," unadj.p = ", round(pvalt1[i],4))
      
      ylab_ = paste0(mecisl[[1]]$gene[i]," unadj.p = ", paste(names(pvalt2l), collapse="/")," = ", paste(sapply(pvalt2l, function(pvalt2) round(pvalt2[i],4)), collapse="/"))
      
      dfi = cbind(dfm,df0[[i]])
      dfi$dna = factor(dfi$dna)
      
      
      p = ggplot(dfi, aes(dna, cont)) +
        geom_boxplot(aes(x=dna, y=cont, colour=response), alpha = 0.5, position=position_dodge(0)) +
        geom_jitter(aes(colour=response), width=0.2) +
        geom_smooth(aes(group=response), method="lm", size = 2, se = F) +
        facet_grid(time~.,scales="free_x")
      
      print(p + xlab(xlab_) + ylab(ylab_) + ggtitle(main_))
      graphics.off()
    }
  })
}
time_output(start)










