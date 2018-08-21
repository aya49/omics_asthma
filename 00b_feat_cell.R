## input: meta data and cells
## output: well formatted cell data and meta files
## aya43@sfu.ca
## created 20180509
## last modified 20180509


## logistics
root = "~/projects/asthma"; commandArgs <- function(...) root  # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
libr(pkgs())


## options
writecsv = T #write results as csv on top of Rdata?


# load meta file with  cell counts
meta_file1_temp0 = fread(meta_file_temp1_dir, data.table=F)
meta_file1_temp0[grepl("WRF",meta_file1_temp0[,"Name"]),"Name"] = "WRF"
meta_file1_temp0[,"Filename"] = gsub(".bam","",meta_file1_temp0[,"Filename"])
meta_file1_temp0$Sample.ID[grepl("L_ST_025",meta_file1_temp0$Sample.ID)] = "L_ST_025"

cells = meta_file1_temp0[,c(35:48,73:97)]

# split data into pre/post allergic attack samples for each subject
cells_pre = cells[meta_file1_temp0$Time=="pre",]
rnpre = meta_file0$id[match(meta_file1_temp0$id_namedate[meta_file1_temp0$Time=="pre"],meta_file0$id_mstlst)]
cells_pre = cells_pre[!is.na(rnpre),]
rownames(cells_pre) = rnpre[!is.na(rnpre)]

cells_post = cells[meta_file1_temp0$Time=="post",]
rnpost = meta_file0$id[match(meta_file1_temp0$Sample.ID[meta_file1_temp0$Time=="post"],meta_file0$id_mstlst)]
cells_post = cells_post[!is.na(rnpost),]
rownames(cells_post) = rnpost[!is.na(rnpost)]

# save
save(cells_pre, file=paste0(feat_cell_dir,".pre.Rdata"))
if (writecsv) write.csv(cells_pre, file=paste0(feat_cell_dir,".pre.csv"))
save(cells_post, file=paste0(feat_cell_dir,".post.Rdata"))
if (writecsv) write.csv(cells_post, file=paste0(feat_cell_dir,".post.csv"))
