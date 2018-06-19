## input: RSEM output .isoforms/genes.results files
## output: raw cell counts matrix
## aya43@sfu.ca
## created 20180509
## last modified 20180522

## root directory
root = "~/projects/asthma"
setwd(root)

dir.create(paste0(root, "/result"), showWarnings=F)

type = "genes" #"isoforms", "genes"
result_dir = paste0(root, "/result/RNAseq_",type); dir.create(result_dir, showWarnings=F)

## input directory
data_dir = paste0(root,"/data/RNAseq")
rsem_dir = paste0(data_dir,"/rsem")
meta_file_temp1_dir = paste0(data_dir,"/RNASeq.asthma.clinical_sequencing_merged.csv")
# meta_file_temp2_dir = paste0(data_dir,"/asthmaDemo_allsite.csv")
meta_col_temp_dir = paste0(data_dir,"/HuGene-2_1-st-v1.na36.hg19.transcript.csv")
meta_col_tr_temp_dir = paste0(data_dir,"/HuGene-2_1-st-v1.na36.hg19.probeset.csv")

## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")
meta_cell_dir = paste0(meta_dir,"/cell")

feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_feature_dir = paste0(feat_dir,"_",type,"-file-featraw")


## options




## load & save matrix -----------------------------------

start = Sys.time()

data_paths = list.files(rsem_dir, pattern=paste0(type,".results"), full.names=T)
data_filenames = fileNames(data_paths, ext=paste0(type,".results"))

if (type=="genes") {
  # counts0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f5 ",x))))
  counts0 = lapply(data_paths, function(x) fread(x, select=5))
  counts = Reduce("cbind",counts0)
  
  # counts_rownames = read.table(pipe(paste0("cut -f1 ",data_paths[1])))[-1,c(2,1)]
  counts_rownames = fread(data_paths[1], select = c(1,2), data.table=F)
  colnames(counts_rownames) = c("id","transcript") #g3n3=id
  save(counts_rownames, file=paste0(meta_col_dir,".Rdata"))
  
  rownames(counts) = counts_rownames$id
  colnames(counts) = data_filenames
  
  save(counts, file=paste0(feat_feature_dir,".Rdata"))
  
} else if (type=="isoforms")  {
  start1 = Sys.time()
  # countsis0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f5 ",x))))
  countsis0 = lapply(data_paths, function(x) fread(x, select=5))
  countsis = Reduce("cbind",countsis0)
  # countsis = foreach(x=data_paths, .combine=cbind) %dopar% { return(fread(x, select=5)) }
  
  # isopct0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f8 ",x))))
  isopct0 = lapply(data_paths, function(x) fread(x, select=8))
  isopct = Reduce("cbind",isopct0)
  # isopct = foreach(x=data_paths, .combine=cbind) %dopar% { return(fread(x, select=8)) }
  time_output(start1)
  
  # counts_rownames = read.table(pipe(paste0("cut -f1,2 ",data_paths[1])))
  counts_rownames = fread(data_paths[1], select = c(1,2), data.table=F)
  colnames(counts_rownames) = c("id","gene") #transcript=id
  save(counts_rownames, file=paste0(meta_col_dir,".Rdata"))
  
  rownames(countsis) = rownames(isopct) = counts_rownames$id
  colnames(countsis) = colnames(isopct) = data_filenames
  
  save(countsis, file=paste0(feat_feature_dir,".Rdata"))
  save(isopct, file=paste0(gsub("featraw","isopctraw",feat_feature_dir),".Rdata"))
}

time_output(start)








## save meta_file; meta_col ----------------------------------------

meta_file1_temp0 = fread(meta_file_temp1_dir, data.table=F)
ucol = col_probe(meta_file1_temp0)
meta_file1_temp1 = meta_file1_temp0[,-ucol$u1]
ucol = col_probe(meta_file1_temp1)

# meta_file2_temp0 = fread(meta_file_temp2_dir, data.table=F)
# ucol = col_probe(meta_file2_temp0)
# meta_file2_temp1 = meta_file2_temp0[,-ucol$u1]
# ucol = col_probe(meta_file2_temp1)



meta_file1_temp = meta_file1_temp1[,c("Filename", "Subject", "Phenotype", "Time", "Allergen", "SITE", 
                                    "RACE", "SEX", "HT.cm.", "AGE", "BLFEV", 
                                    "Volume", "Concentration", "Quantity", "Extracted")]
colnames(meta_file1_temp) = c("fileName", "sample", "response", "time", "allergen", "centre", 
                              "race", "sex", "height", "age", "blfev", 
                              "vol", "conc", "qty", "extracted")
meta_file1_temp$fileName = gsub(".bam","",meta_file1_temp[,"fileName"])

meta_cell = cbind(meta_file1_temp1[meta_file1_temp$fileName,c(31:ncol(meta_file1_temp1))])

meta_file = meta_file1_temp[match(data_filenames,meta_file1_temp$fileName),]
meta_cell = meta_cell[match(data_filenames,meta_cell$fileName),]

save(meta_file, file=paste0(meta_file_dir,".Rdata"))
write.csv(meta_file, file=paste0(meta_file_dir,".csv"))
save(meta_cell, file=paste0(meta_cell_dir,".Rdata"))














