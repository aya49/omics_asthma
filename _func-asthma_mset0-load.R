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

# index all features
m0_inds = lapply(names(m0s), function(feat_type) {
  m0 = m0s[[feat_type]]
  fin = lapply(names(file_inds), function(file_ind_n) {
    file_ind = file_inds[[file_ind_n]]
    if (file_ind_n!="all" & all(rownames(m0)%in%file_ind)) nextm = T
    if (file_ind_n=="all") file_ind = rownames(m0)
    
    # get rid of rows/cols with too many na or is homogenous
    mind = match(file_ind, rownames(m0))
    m = as.matrix(m0)[mind[!is.na(mind)],, drop=F]
    class(m) = "numeric"
    
    # get row/col indices
    mcol_ind = colnames(m)[apply(m[apply(m,1,function(x) any(!is.na(x))),,drop=F],2,
                                 function(x) sum(!is.na(x))>(good_na*length(x)))]
    mrow_ind = rownames(m)[apply(m[,mcol_ind,drop=F], 1, function(x) any(!is.na(x)))]
    if (grepl("dna",feat_type)) mcol_ind = 
      mcol_ind[apply(m[mrow_ind,mcol_ind,drop=F], 2,function(x) min(table(x))>good_col & length(table(x))>1)]
    m_ = m[mrow_ind,mcol_ind,drop=F]
    if (any(dim(m_)==0) | sum(!is.na(m_))<sum(is.na(m_))) return(NULL)
    mcol_ind = mcol_ind[apply(m_,2,function(x) length(unique(x[!is.na(x)]))>1)]
    if (length(mcol_ind)==0) return(NULL)
    m_ = m_[,mcol_ind,drop=F]
    
    # are columns continuous values?
    m_cont_col = apply(m_,2,function(x) length(unique(x))>cont_col)
    
    # scale columns?
    if (exists("scale_cont")) {if (scale_cont) m[,m_cont_col] = scale(as.numeric(as.matrix(m[,m_cont_col])))}
    
    return(list(row_ind=mrow_ind, col_ind=mcol_ind))#, cont_col=m_cont_col))
  })
  names(fin) = names(file_inds)
  fin = fin[!sapply(fin, is.null)]
  if (length(fin)==0) return(NULL)
  if (length(fin)>1) {
    keepfin = 1
    for (fini in 2:length(fin)) {
      if (all(sapply(fin[1:(fini-1)], function(x) !identical(x$row_ind,fin[[fini]]$row_ind))))
        keepfin = append(keepfin,fini)
    }
    fin = fin[keepfin]
  }
  return(fin)
})
names(m0_inds) = names(m0s)
m0_inds = m0_inds[!sapply(m0_inds, is.null)]
