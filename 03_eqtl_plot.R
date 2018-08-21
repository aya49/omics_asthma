## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614


## logistics
root = "~/projects/asthma"; commandArgs <- function(...) root  # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
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
feat_types = lapply(str_split(eqtl_names,"_"), function(x) c(str_split(x,"[-]")[[1]][1], str_split(x,"[-]")[[1]][2]))


foreach (ei = 1:length(eqtl_paths)) %dopar% {
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
    if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiments
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
        if (all(mecis$gene%in%meta_col20[,id_col])) {
          mcorder = match(mecis$gene,meta_col20[,id_col])
          meta_col2 = meta_col20[mcorder,]
          mecis$gene = meta_col2[,"symbol"]
          save(mecis, file=eqtl_path)
        } else {
          mcorder = match(mecis$gene,meta_col20[,"symbol"])
          meta_col2 = meta_col20[mcorder,]
        }
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
    mecis_extra = cbind(meta_col1, mecis, meta_col2)
    write.csv(mecis, file=gsub(".Rdata",".csv",eqtl_path))
    write.csv(mecis_extra, file=gsub(".Rdata","_withmetadata.csv",eqtl_path))
    
    
    
    # plot the histogram of local and distant p-values
    # png(gsub(".Rdata","_qq.png",eqtl_path), width=width, height=height)
    # plot(me)
    # graphics.off()
    
    dir.create(gsub(".Rdata","",eqtl_path), showWarnings=F)
    for (i in 1:min(nrow(mecis),max_plots)) {
      pname = paste0(gsub(".Rdata","",eqtl_path),"/",mecis$snps[i],"-",mecis$gene[i],".png")
      if (!overwrite & file.exists(pname)) next()
      png(file=pname, width=width, height=height)
      
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
      
      if (feat_type[1]=="genotype" ) {
        p = ggplot(data.frame(a=factor(m1[,i]),b=m2[,i]), aes(a, b)) + 
          geom_boxplot(aes(colour=class_n),position=position_dodge(0)) + 
          geom_smooth(aes(group=class_n, colour=class_n), method=lm, size = 2, se = F) +
          stat_summary(fun.y=mean, geom="line") + 
          geom_jitter(width=0.2,aes(colour=class_n))
        
        # ggplot(data.frame(class_n=class_n, m1i=m1[,i], m2i=m2[,i]), 
        #        aes(x=interaction(class_n, m1i), y=m2i)) +
        #   geom_boxplot(aes(fill = class_n), alpha = 0.5) +
        #   geom_line(aes(group = m1i),
        #             alpha = 0.5, colour = "darkgrey")
        # facet_grid(sp~variable,scales="free_x")
      } else {
        p = ggplot(data.frame(a=m1[,i],b=m2[,i]), aes(a, b)) + 
          geom_point(aes(colour=class_n)) +
          geom_smooth(aes(group=class_n, colour=class_n), method=lm)
      }
      print(p + xlab(xlab_) + ylab(ylab_) + ggtitle(main_))
      graphics.off()
    }
    
    
  })
}
time_output(start)










