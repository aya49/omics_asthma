## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/src/_dirs.R"))
source(paste0(root, "/src/_func.R"))
source(paste0(root, "/src/_func-asthma.R"))
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


# load features
m0s = lapply(feat_types, function(feat_type) get(load(paste0(feat_dir,"/",feat_type,".Rdata"))))
names(m0s) = feat_types
feat_names = unique(sapply(feat_types, function(x) gsub("[0-9]|[.]pre|[.]post|[.]diff","",x)))

# load meta_cols
m_col0s = lapply(feat_names, function(feat_name) {
  fn = paste0(meta_col_dir,feat_name,".Rdata")
  if (file.exists(fn)) return(get(load(fn)))
  return(NULL)
})
names(m_col0s) = feat_names


## plot eQTL ---------------------------------------------

start = Sys.time()

pval_paths = list.files(gwas_dir, full.names=T, pattern=".Rdata")
eqtl_paths = list.files(eqtl_dir, full.names=T, pattern=".Rdata")
eqtl_names = fileNames(eqtl_paths)
feat_type_sets = lapply(str_split(eqtl_names,"_"), function(x) gsub("[0-9]","", c(str_split(x,"[-]")[[1]][1], str_split(x,"[-]")[[1]][2])))
feat_names = lapply(feat_type_sets, function(x) sapply(str_split(x,"[.]"), function(y) y[1]))

class_coli = 1 # response
class_col = class_cols[class_coli]

# for (ei in 1:length(eqtl_paths)) { try({
a = foreach (ei = 1:length(eqtl_paths)) %dopar% { try({
  # a = foreach (ei = 1:length(eqtl_paths)) %dopar% {
  eqtl_path = eqtl_paths[ei]; dir.create(gsub(".Rdata","",eqtl_path), showWarnings=F)
  eqtl_name = eqtl_names[ei]
  feat_type_set = feat_type_sets[[ei]]
  # bins = ""
  # if (grepl("dna",feat_type_set[1])) {
  #   if (grepl("01", feat_type_set[1])) bins = "01"
  #   if (grepl("12", feat_type_set[1])) bins = "12"
  #   feat_type_set[1] = gsub("[.]01|[.]12","",feat_type_set[1])
  # }
  feat_name = feat_names[[ei]]
  file_ind_n = gsub("[-X]","",str_extract(eqtl_name,"-[a-z]+X"))
  
  file_ind = unlist(read.csv(gsub(".Rdata","_id.csv",eqtl_path)))
  # file_ind = file_inds[[file_ind_n]]
  # if (file_ind_n=="all") file_ind = meta_file0[,id_col]
  # meta_file = meta_file0[meta_file0[,id_col]%in%file_ind & meta_file0[,id_col]%in%row_ind,]
  # if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiments[class_coli]
  
  
  # load eqtl
  mecis = NULL
  try({ mecis = get(load(eqtl_path)) }); if (is.null(mecis)) next()
  colnames(mecis)[c(1:2)] = feat_type_set; save(mecis, file=eqtl_path) # temporary
  mecis[,c(1:2)] = apply(mecis[,c(1:2)],2,as.character)
  
  mecis_extra = mecis
  pvalt = m_set = meta_col_set = NULL; 
  for (xi in 1:length(feat_type_set)) {
    x = feat_type_set[xi]
    
    # load gwas
    pn = pval_paths[grepl(paste0(x,"-",file_ind_n,"Xall"),pval_paths)]
    if (file.exists(pn)) {
      pvali = get(load(pn))
      pvalt[[x]] = pvalti = pvali[match(mecis[,xi],rownames(pvali)), grepl("none$",colnames(pvali))]
      names(pvalti) = rownames(mecis)
      # if (length(pvalti)==0) {
      mecis_extra = cbind(mecis_extra, pvalti)
      colnames(mecis_extra)[ncol(mecis_extra)] = paste0("gwas.p_",x)
      # }
    }
    
    # load feature matrix
    m_seti = m0s[[x]]
    m_set[[x]] = m_seti[file_ind,mecis[,xi]]
    
    # load meta_col & trim
    if (!is.null(m_col0s[[feat_name[xi]]])) {
      meta_col_setx = m_col0s[[feat_name[xi]]]
      meta_col_set[[x]] = meta_col_setx1 = meta_col_setx[match(mecis[,xi], meta_col_setx[,id_col]),,drop=F]
      
      colnames(meta_col_setx1) = paste0(colnames(meta_col_setx1),"_",feat_name[xi])
      mecis_extra = cbind(mecis_extra,meta_col_setx1)
    }
  }
  cole = feat_type_set%in%names(meta_col_set)
  
  # # adjust m_set[[1]] if needed
  # if (bins=="01") m_set[[1]][m_set[[1]]==2] = 1
  # if (bins=="12") m_set[[1]][m_set[[1]]==0] = 1
  
  # prepare meta_file
  meta_file = meta_file0[match(file_ind, meta_file0[,id_col]),]
  class = as.numeric(factor(meta_file[,class_col]))
  class_name = as.character(levels(factor(meta_file[,class_col])))
  class_n = as.character(meta_file[,class_col])

  # save eqtl as table
  # write.csv(mecis, file=gsub(".Rdata",".csv",eqtl_path))
  write.csv(mecis_extra, file=gsub(".Rdata",".csv",eqtl_path))
  
  # plot the histogram of local and distant p-values
  # png(gsub(".Rdata","_qq.png",eqtl_path), width=width, height=height)
  # plot(me)
  # graphics.off()
  
  # plot out eqtl one by one; make labels symbol or rs numbers if possible
  dir.create(gsub(".Rdata","",eqtl_path), showWarnings=F)
  id1 = mecis[,1]; id2 = mecis[,2]
  if (grepl("dna", feat_type_set[1])) {
    id1 = meta_col_set[[feat_type_set[1]]][,"dbSNP"]
    id1[is.na(id1)] = mecis[is.na(id1),1]
  }
  if (grepl("rna", feat_type_set[2])) {
    id2 = meta_col_set[[feat_type_set[2]]][,"symbol"]
    id2[is.na(id2)] = mecis[is.na(id2),2]
  }
  for (i in 1:min(nrow(mecis),max_plots)) {
    pname = paste0(gsub(".Rdata","",eqtl_path),"/",
                   str_pad(i, 3, pad="0"),"_",mecis[i,1],"_",mecis[i,2],".png")
    if (!overwrite & file.exists(pname)) next()
    
    # titles and labels
    main_ = paste0("eqtl: stat=",round(mecis$statistic[i],4), "; unadj.p=", round(mecis$pvalue[i],4), "; fdr.p=", round(mecis$FDR[i],4))
    
    xlab_more = ylab_more = ""
    if(grepl("dna",feat_type_set[1])) {
      f1_colind = grepl(paste0(c("chromosome","pos_phys"),"_",feat_name[1],collapse="|"),colnames(mecis_extra))
      xlab_more = 
        paste0("\n", paste(gsub(paste0("_",feat_name[1]),"",colnames(mecis_extra)[f1_colind]), collapse="/"),
               ": ", paste(mecis_extra[i,f1_colind], collapse=" / "))
    } 
    if(grepl("rna",feat_type_set[2])) {
      f2_colind = grepl(paste0(c("chr","start","end"),"_",feat_name[2],collapse="|"),colnames(mecis_extra))
      ylab_more = 
        paste0("\n", paste(gsub(paste0("_",feat_name[2]),"",colnames(mecis_extra)[f2_colind]), collapse="/"),
               ": ", paste(mecis_extra[i,f2_colind], collapse=" / "))
    } 
    xlab_ = paste0(feat_type_set[1],": ", id1[i], "\nunadj.p = ", round(pvalt[[1]][i],4), xlab_more)
    ylab_ = paste0(feat_type_set[2],": ", id2[i], "\nunadj.p = ", round(pvalt[[2]][i],4), ylab_more)
    
    # plot
    if (grepl("dna",feat_type_set[2]) ) {
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
    png(file=pname, width=width, height=height)
    print(p + xlab(xlab_) + ylab(ylab_) + ggtitle(main_))
    graphics.off()
  }

}) }
time_output(start)










