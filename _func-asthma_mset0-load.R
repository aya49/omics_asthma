m_set0 = lapply(feat_type_set, function(x) get(load(paste0(feat_dir,"/",x,".Rdata"))) )
feat_name = sapply(strsplit(feat_type_set,"[.]"), function(x) x[1])

col_ind_n = "all"
