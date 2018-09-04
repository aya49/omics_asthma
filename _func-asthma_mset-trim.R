class_col = class_cols[class_coli]
control = controls[class_coli]
experiment = experiments[class_coli]

## trim matrices
ftsri = Reduce(intersect,lapply(feat_type_set, function(x) m0_inds[[x]][[file_ind_n]]$row_ind))
if (length(ftsri)==0) next()
m_set = NULL
for(ftsi in 1:length(feat_type_set)) {
  feat_type = feat_type_set[ftsi]
  fni = feat_name[ftsi]
  m0 = m0s[[feat_type]][ftsri, m0_inds[[feat_type]][[file_ind_n]]$col_ind,drop=F]
  
  # get row/col indices
  m = m0[,apply(m0,2, function(x) length(unique(x[!is.na(x)]))>1)]
  if (grepl("dna",feat_type)) 
    m = m[,apply(m, 2,function(x) min(table(x))>good_col & length(table(x))>1)]
  if (ncol(m)==0) next()
  
  m_set[[feat_type]] = m
}


cat(" (",file_ind_n, " x ",col_ind_n,"; ",paste(sapply(m_set,nrow),collapse="/")," x ",paste(sapply(m_set,ncol),collapse="/"),") ", sep="")

# prepare meta_file
# get row files/samples & covariate
meta_file = meta_file0[match(rownames(m_set[[1]]),meta_file0[,id_col]),]
# if (file_ind_n=="flipperdr") meta_file[meta_file[,flipper_col],class_col] = experiment
if (file_ind_n=="flipperdr") meta_file[!meta_file[,id_col]%in%goodppl,class_col] = experiment
