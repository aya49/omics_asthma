## input: load pancancer, elements, metab data
## output: reformatted data
## aya43@sfu.ca
## created 20170807
## last modified 20180614

## root directory
root = "~/projects/asthma"
# setwd(root)
result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)
meta_dir = paste0(result_dir, "/meta")
meta_file_dir = paste0(meta_dir, "/file")
meta_col_dir = paste0(meta_dir, "/col")

feat_dir = paste0(result_dir, "/feat");

library("annotables")
library("biomaRt")
library("stringr")

id_col = "id"

## input directory
data_dir = paste0(root, "/data")

meta_file = get(load(paste0(meta_file_dir,".Rdata")))
chr_map = read.table(paste0(data_dir,"/GRCh37_UCSC2ensembl.txt"))
iddates = paste0(meta_file$id, "_", meta_file$date)





## pan cancer
load(paste0(data_dir, "/RNApancancer/data/preprocessedData/panCancerDatasets.RDATA"))

meta_file_pc = demoPre
ids_pc = as.character(meta_file_pc$NAME)
ids_pc[ids_pc=="WF (WRF)"] = "WRF"
iddates_pc = paste0(ids_pc,"_",as.character(meta_file_pc$AIC_YMD))
iddates_pc[!iddates_pc%in%iddates]
meta_file_pc_d = meta_file_pc[!iddates_pc%in%iddates,]
write.csv(meta_file_pc_d, file=paste0(meta_file_dir,".rnapancancer.setdiff.Rdata"))

feat_pc = t(genEset)
rownames(feat_pc) = ids_pc
save(feat_pc, file=paste0(feat_dir, "/rnapcgenes.pre.Rdata"))
write.csv(feat_pc, file=paste0(feat_dir, "/rnapcgenes.pre.csv"))

meta_col = read.csv(paste0(data_dir, "/Cells_nCounter_Human_PanCancer_Immune_Profiling_Panel_Gene_List.csv"), row.names=1)
colnames(meta_col) = c("name_cd", id_col, "probe_ns", "gene_class", "cell", "ir_type", "trait_class", "trait", "synonyms")
grch38_dt <- merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")
meta_col_pc = meta_col[match(colnames(feat_pc), as.character(meta_col$id)),]
meta_col_pc$id = colnames(feat_pc)
meta_col_pc = cbind(meta_col_pc, 
                    grch38_dt[match(meta_col_pc[,id_col], grch38_dt$symbol),c("symbol","chr","start","end","biotype","description")])
meta_col_pc$symbol = meta_col_pc$id

# convert chromosomes with mappings here: https://github.com/dpryan79/ChromosomeMappings/blob/58f75d5a6b2d4d7274c7246d10db666033926ee5/GRCh37_UCSC2ensembl.txt
meta_col_pc_ = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position', 'gene_biotype', "description"),
      filters = 'hgnc_symbol',
      values = meta_col_pc$symbol,
      mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                     host="grch37.ensembl.org", 
                     path="/biomart/martservice", 
                     dataset="hsapiens_gene_ensembl"))
matchorder = match(meta_col_pc$symbol,meta_col_pc_$hgnc_symbol)
meta_col_pc[!is.na(matchorder),c("id","symbol","chr","start","end","biotype","description")] = meta_col_pc_[matchorder[!is.na(matchorder)],]

save(meta_col_pc, file=paste0(meta_col_dir,"-rnapcgenes.Rdata"))
write.csv(meta_col_pc, file=paste0(meta_col_dir,"-rnapcgenes.csv"))




## gene/isoform meta_cols (add columns from pan cancer and annotatables)

col_paths = list.files(meta_dir, full.names=T)
col_g_paths = col_paths[grepl("genes|iso", col_paths) & !grepl(".raw", col_paths) & 
                          grepl("col-", col_paths) &  grepl(".Rdata", col_paths)]

grch38_dt <- merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")
for (col_g_path in col_g_paths) {
  meta_col_temp = get(load(col_g_path))
  if (any(grepl("trait",colnames(meta_col_temp)))) next()
  
  ensg_col = apply(meta_col_temp, 2, function(x) any(grepl("ENSG", x)))
  enst_col = apply(meta_col_temp, 2, function(x) any(grepl("ENST", x)))
  if (any(enst_col)) {
    enst_colind = which(enst_col)[1]
    meta_col_temp = cbind(meta_col_temp, 
                          grch38_dt[match(str_extract(meta_col_temp[,enst_colind], "ENST[0-9]+"), grch38_dt$enstxp),
                                    c("symbol","chr","start","end","biotype","description")])
  } else if (any(ensg_col)) {
    ensg_colind = which(ensg_col)[1]
    meta_col_temp = cbind(meta_col_temp, 
                          grch38_dt[match(str_extract(meta_col_temp[,enst_colind], "ENSG[0-9]+"), grch38_dt$ensgene),
                                    c("symbol","chr","start","end","biotype","description")])
  }
  if (any(grepl("symbol", colnames(meta_col_temp))) &
      !any(grepl("trait", colnames(meta_col_temp))))
    meta_col_temp = cbind(meta_col_temp, meta_col_pc[match(meta_col_temp$symbol, meta_col_pc$id),])
  
  save(meta_col_temp, file=col_g_path)
  write.csv(meta_col_temp, file=gsub("Rdata","csv",col_g_path))
}



