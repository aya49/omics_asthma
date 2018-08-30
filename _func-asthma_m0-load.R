# load matrix and adjust column indices to be used
m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))

col_inds_i = sapply(names(col_inds0), function(x) grepl(x,feat_type))
if (sum(col_inds_i)==0) { col_inds = list(all="")
} else { col_inds = col_inds0[[which(col_inds_i)]] }
  