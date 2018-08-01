## input: meta data and cells
## output: well formatted data and meta files
## aya43@sfu.ca
## created 20180509
## last modified 20180509



## root directory
root = "~/projects/asthma"
setwd(root)
result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)


## input directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")

rnaseq_dir = paste0(root,"/data/RNAseq")
meta_file_temp1_dir = paste0(rnaseq_dir,"/RNASeq.asthma.clinical_sequencing_merged.csv")


## output directory
feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_cell_dir = paste0(feat_dir,"/cell")

## libraries
source("code/_func.R")
libr(c("data.table", "Matrix", "stringr"))


## options
writecsv = T #write results as csv on top of Rdata?
# save cells



## start

meta_file = get(load(paste0(meta_file_dir,".Rdata")))

meta_file1_temp0 = fread(meta_file_temp1_dir, data.table=F)
meta_file1_temp0[grepl("WRF",meta_file1_temp0[,"Name"]),"Name"] = "WRF"
meta_file1_temp0[,"Filename"] = gsub(".bam","",meta_file1_temp0[,"Filename"])
meta_file1_temp0$Sample.ID[grepl("L_ST_025",meta_file1_temp0$Sample.ID)] = "L_ST_025"

cells = meta_file1_temp0[,c(35:48,73:97)]

cells_pre = cells[meta_file1_temp0$Time=="pre",]
rnpre = meta_file$id[match(meta_file1_temp0$id_namedate[meta_file1_temp0$Time=="pre"],meta_file$id_mstlst)]
cells_pre = cells_pre[!is.na(rnpre),]
rownames(cells_pre) = rnpre[!is.na(rnpre)]

cells_post = cells[meta_file1_temp0$Time=="post",]
rnpost = meta_file$id[match(meta_file1_temp0$Sample.ID[meta_file1_temp0$Time=="post"],meta_file$id_mstlst)]
cells_post = cells_post[!is.na(rnpost),]
rownames(cells_post) = rnpost[!is.na(rnpost)]

save(cells_pre, file=paste0(feat_cell_dir,".pre.Rdata"))
if (writecsv) write.csv(cells_pre, file=paste0(feat_cell_dir,".pre.csv"))
save(cells_post, file=paste0(feat_cell_dir,".post.Rdata"))
if (writecsv) write.csv(cells_post, file=paste0(feat_cell_dir,".post.csv"))
