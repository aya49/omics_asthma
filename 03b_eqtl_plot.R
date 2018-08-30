## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
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

coln = list(c("dbSNP","chromosome","pos_phys"), c("symbol","chr","start","end"))

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

class_coli = 1 # response

foreach (ei = 1:length(eqtl_paths)) %dopar% {
  try({
    # a = foreach (ei = 1:length(eqtl_paths)) %dopar% {
    eqtl_path = eqtl_paths[ei]
    eqtl_name = eqtl_names[ei]
    feat_type_set = feat_type_sets[[ei]]
    # bins = ""
    # if (grepl("dna",feat_type_set[1])) {
    #   if (grepl("01", feat_type_set[1])) bins = "01"
    #   if (grepl("12", feat_type_set[1])) bins = "12"
    #   feat_type_set[1] = gsub("[.]01|[.]12","",feat_type_set[1])
    # }
    feat_name = feat_names[[ei]]
    
    file_ind = unlist(read.csv(gsub(".Rdata","_id.csv",eqtl_path)))
    file_ind_n = gsub("[-]|X","",str_extract(eqtl_name,"-[a-zA-Z]+X"))
    # file_ind = file_inds[[file_ind_n]]
    # if (file_ind_n=="all") file_ind = meta_file0[,id_col]
    # meta_file = meta_file0[meta_file0[,id_col]%in%file_ind & meta_file0[,id_col]%in%row_ind,]
    # if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiments[class_coli]
    
    
    # load eqtl
    mecis = NULL
    try({ mecis = get(load(eqtl_path)) })
    if (is.null(mecis)) next()
    # mecis = me$cis$eqtls
    mecis$snps = as.character(mecis$snps)
    mecis$gene = as.character(mecis$gene)
    
    # load matrix
    m_set = lapply(feat_type_set, function(x) get(load(paste0(feat_dir,"/",x,".Rdata")))[file_ind,])
    # row_ind = sort(intersect(rownames(m_set[[1]]), rownames(m_set[[2]])))
    names(m_set) = feat_type_set
    
    # # adjust m_set[[1]] if needed
    # if (bins=="01") m_set[[1]][m_set[[1]]==2] = 1
    # if (bins=="12") m_set[[1]][m_set[[1]]==0] = 1
    
    # load gwas
    pval = lapply(feat_type_set, function(x) get(load(pval_paths[grepl(paste0(x,"-",file_ind_n,"Xall"),pval_paths)])))
    pvalt = lapply(1:length(pval), function(xi) {
      x = pval[[xi]]; y = x[,grepl("_none$",colnames(x))]; names(y) = rownames(x)
      y[match(colnames(m_set[[xi]]),names(y))]
    })
    
    # add p values to eqtl cis table
    mecis = cbind(mecis, Reduce("cbind", pvalt))
    colnames(mecis)[c(ncol(mecis),ncol(mecis)-1)] = paste0("p.",feat_type_set)
    
    # prepare meta_file
    meta_file = meta_file0[match(file_ind, meta_file0[,id_col]),]
    class = as.numeric(factor(meta_file[,class_col]))
    class_name = as.character(levels(factor(meta_file[,class_col])))
    class_n = as.character(meta_file[,class_col])
    
    # load meta_col & trim
    meta_col_set = meta_col_ids = NULL; mecis_extra = mecis
    for (x in 1:length(feat_name)) {
      file_name = paste0(meta_col_dir,feat_name[x],".Rdata")
      if (!file.exists(file_name)) next()
      meta_col_setx = get(load(file_name))
      if (!all(coln[[x]]%in%colnames(meta_col_setx))) next()
      meta_col_setx1 = meta_col_setx[match(colnames(m_set[[feat_type_set[x]]]), meta_col_setx[,id_col]),]
      
      if (x==1) { mecis_ids = mecis$snps
      } else { mecis_ids = mecis$gene }
      meta_col_ids[[x]] = rep(NA,ncol(m_set[[x]]))
      if (grepl("dna",feat_name[x])) {
        meta_col_ids[[x]] = meta_col_setx1$dbSNP
      } else if (grepl("rna",feat_name[x])) {
        meta_col_ids[[x]] = meta_col_setx1$symbol
      }
      meta_col_ids[[feat_type_set[x]]][is.na(meta_col_ids[[x]])] = colnames(m_set[[x]])[is.na(meta_col_ids[[x]])]
      
      meta_col_set[[feat_type_set[x]]] = meta_col_setx1[match(mecis_ids,meta_col_ids[[x]]),,drop=F]
      colnames(meta_col_setx) = paste0(colnames(meta_col_setx),".",feat_name[x])
      mecis_extra = cbind(mecis_extra,meta_col_set)
    }
    cole = feat_type_set%in%names(meta_col_set)
    
    # trim matrices
    for (x in 1:length(feat_name)) {
      if (x==1) { mecis_ids = mecis$snps
      } else { mecis_ids = mecis$gene }
      if (!all(mecis_ids%in%colnames(m_set[[x]]))) {
        m_set[[x]] = m_set[[x]][,meta_col_set[[feat_type_set[x]]][,id_col],drop=F]
      } else { m_set[[x]] = m_set[[x]][,match(mecis_ids,colnames(m_set[[x]])),drop=F] }
    }
    
    # save eqtl as table
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
      if(grepl("dna",feat_type_set[1])) {
        f1_colind = grepl(c("chromosome","pos_phys"),colnames(mecis_extra))
        xlab_more = paste0("\n", paste(gsub(paste0(".",feat_type_set[1]), "",colnames(mecis)[f1_colind]), collapse="/"), 
                           ": ", paste(mecis_extra[i,f1_colind], collapse=" / "))
      } 
      xlab_ = paste0(feat_type_set[1],": ", mecis$snps[i], "; unadj.p = ", round(pvalt[[1]][i],4), xlab_more)
      
      ylab_more = ""
      if(grepl("rna",feat_type_set[2])) {
        f2_colind = grepl(c("chr","start","end"),colnames(mecis_extra))
        ylab_more = paste0("\n", paste(gsub(paste0(".",feat_type_set[2]),"",colnames(mecis)[f2_colind]), collapse="/"), ": ", 
                           paste(mecis_extra[i,f2_colind], collapse=" / "))
      } 
      ylab_ = paste0(feat_type_set[2],": ", mecis$gene[i], "; unadj.p = ", round(pvalt[[2]][i],4), ylab_more)
      
      if (feat_type_set[1]=="dna" ) {
        p = ggplot(data.frame(a=factor(m_set[[1]][,i]),b=m_set[[2]][,i]), aes(a, b)) + 
          geom_boxplot(aes(colour=class_n),position=position_dodge(0)) + 
          geom_smooth(aes(group=class_n, colour=class_n), method=lm, size = 2, se = F) +
          stat_summary(fun.y=mean, geom="line") + 
          geom_jitter(width=0.2,aes(colour=class_n))
        
        # ggplot(data.frame(class_n=class_n, m1i=m_set[[1]][,i], m2i=m_set[[2]][,i]), 
        #        aes(x=interaction(class_n, m1i), y=m2i)) +
        #   geom_boxplot(aes(fill = class_n), alpha = 0.5) +
        #   geom_line(aes(group = m1i),
        #             alpha = 0.5, colour = "darkgrey")
        # facet_grid(sp~variable,scales="free_x")
      } else {
        p = ggplot(data.frame(a=m_set[[1]][,i],b=m_set[[2]][,i]), aes(a, b)) + 
          geom_point(aes(colour=class_n)) +
          geom_smooth(aes(group=class_n, colour=class_n), method=lm)
      }
      print(p + xlab(xlab_) + ylab(ylab_) + ggtitle(main_))
      graphics.off()
    }
    
    
  })
}
time_output(start)










