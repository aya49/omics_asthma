class_col = class_cols[class_coli]
control = controls[class_coli]
experiment = experiments[class_coli]

m_set = nextmm = NULL
for (x in 1:length(m_set0)) {
  feat_type = feat_type_set[x]
  m0 = as.matrix(m_set0[[x]])[as.character(meta_file[,id_col]),, drop=F]
  col_inds_i = sapply(names(col_inds0), function(x) grepl(x,feat_type))
  if (sum(col_inds_i)==0) { col_inds = list(all="")
  } else { col_inds = col_inds0[[which(col_inds_i)]] }
  
  
  nextm = F
  
  # prep m
  m1 = m0[rownames(m0)%in%meta_file0[,id_col],, drop=F]
  # if (feat_type!="dna" & f1_bin!="") nextm = T
  
  file_ind = file_inds[[file_ind_n]]
  if (file_ind_n!="all" & all(rownames(m1)%in%file_ind)) nextm = T
  if (file_ind_n=="all") file_ind = rownames(m1)
  # file_ind_n = paste0("-",file_ind_n)
  
  ## ensure there is a class label for each sample
  # file_ind = file_ind[!is.na(meta_file0[match(file_ind,meta_file0[,id_col]),class_col])]
  
  col_ind = col_inds[[col_ind_n]]
  if (file_ind_n!="all" & all(colnames(m1)%in%col_ind)) nextm = T
  if (col_ind_n=="all") col_ind = colnames(m1)
  # col_ind_n = paste0(".",col_ind_n)
  
  m = as.matrix(m1)[rownames(m1)%in%file_ind, colnames(m1)%in%col_ind]
  class(m) = "numeric"
  m = m[apply(m,1,function(x) any(!is.na(x))),,drop=F]
  m = m[,apply(m,2,function(x) sum(!is.na(x))>(good_na*length(x))),drop=F]
  if (grepl("dna",feat_type)) {
    # if (f1_bin=="01") m[m==2] = 1
    # if (f1_bin=="12") m[m==0] = 1
    m = m[,apply(m,2,function(x) min(table(x))>good_col & length(table(x))>1)]
  } 
  if (any(dim(m)==0) | sum(!is.na(m))<sum(is.na(m))) {
    cat(" (empty or too many NA's skipped!)"); nextm = T
  }
  m = m[,apply(m,2,function(x) length(unique(x[!is.na(x)]))>1)]
  if (ncol(m)==0) nextm = T
  
  m_cont_col = apply(m,2,function(x) length(unique(x))>cont_col)
  
  m[is.na(m)] = -1
  
  cat(" (",file_ind_n, " x ",col_ind_n,"; ",nrow(m)," x ",ncol(m),") ", sep="")
  
  
  if (ifelse(exists("scale_cont"),scale_cont,F)) m[,m_cont_col] = scale(as.numeric(as.matrix(m[,m_cont_col])))
  m_set[[feat_type]] = m
  nextmm = append(nextmm,nextm)
}
if (!sum(!nextmm)>1) next()
m_set = lapply(m_set,function(x) x[Reduce("intersect",lapply(m_set, rownames)),])
