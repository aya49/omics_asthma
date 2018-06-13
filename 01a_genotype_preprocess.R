## input: RSEM output .isoforms/genes.results files
## aya43@sfu.ca
## created 20180509
## last modified 20180509

## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result/genotype"); dir.create(result_dir, showWarnings=F)

## input directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")
# meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_genotyperaw_dir = paste0(feat_dir,"_snp-file-genotyperaw")

## output directory
feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")
# stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)

## libraries
source("code/_func.R")
# libr("data.table")
# libr("limma")
# libr("edgeR")
# libr("foreach")
# libr("doMC")
libr("stringr")
libr("Matrix")

## options
writecsv = F #write results as csv on top of Rdata?
good_col = 3 #each gene must have > good_col samples of each genotype; else delete
good_col_nap = .5 #each gene must have > good_col_na x 100% samples with data
# good_count = 10 #each gene must have >10 abundence in more than half the samples; else delete
id_col = "fileName"

## load data
meta_file = get(load(paste0(meta_file_dir,".Rdata")))
# meta_col = get(load(paste0(meta_col_dir,".Rdata")))
m0 = get(load(paste0(feat_genotyperaw_dir,".Rdata")))
if (sum(colnames(m0)%in%meta_file[,id_col])>=ncol(m0)) m0 = t(m0)

good_col_na = good_col_nap * nrow(m0)
col_ind = sapply(1:ncol(m0), function(xi) {
  x = m0[,xi]
  ind_x = !is.na(x)
  if (sum(ind_x) >= good_col_na) {
    x = x[ind_x]
    a = min(table(x))>good_col & length(unique(x))>1
  } else {
    a = F
  }
  return(a)
})
# for (xi in 1:ncol(m0)) {
#   x = m0[,xi]
#   ind_x = !is.na(x)
#   if (sum(ind_x) >= good_col_na) {
#     x = x[ind_x]
#     a = min(table(x))>good_col & length(unique(x))>1
#   } else {
#     a = F
#   }
# }

m = m0[,col_ind]

save(m, file=paste0(feat_genotype_dir,".Rdata"))
if (writecsv) write.csv(m, file=paste0(feat_genotype_dir,".csv"))