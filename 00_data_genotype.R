## input: Asthma CEL files
## aya43@sfu.ca
## created 20180509
## last modified 20180509



## root directory
root = "H:/asthma"
setwd(root)
result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)



## input directory
gt1_call_dir = paste0("d:/asthma/genotype/output/AxiomGT1.calls.txt")

probe_to_snp_annotation_dir = paste0("d:/asthma/genotype/Axiom_PMRA.na35.annot.csv")

gt1_meta_dir = paste0("d:/asthma/genotype/output/additional_sample_data.txt")
meta_file_temp_dir = paste0("d:/asthma/genotype/meta_file_temp.csv")



## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")



## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
library(data.table)
library(entropy)
library(foreach)
library(doMC)
library(stringr)
library(Matrix)
source(paste0(root,"/codes/_func.R"))



## options
no_cores = 15#detectCores()-3
registerDoMC(no_cores)

writecsv = T #write results as csv on top of Rdata?



start = Sys.time()

## load data --------------------------------------
gt1_calls = read.table(gt1_call_dir, sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1)

probe_to_snp_annotation = fread(probe_to_snp_annotation_dir)

gt1_meta = as.matrix(fread(gt1_meta_dir))
meta_file_temp = as.matrix(fread(meta_file_temp_dir))




## make meta_file ----------------------------------
meta_file = cbind(meta_file_temp[, !colnames(meta_file_temp)%in%c("Inferred Gender","Random order", "Number", 
                                                                  "1000 Genome Concordance", "Reproducibility", "Replicate Info", 
                                                                  "Cluster CR (Axiom inlier cluster call rate)", "dishQC")], 
                  gt1_meta[match(meta_file_temp[,"Sample ID"],gt1_meta[,"Name"]), !colnames(gt1_meta)%in%c("Name")])
meta_file[meta_file[,"Kit"]=="Mini kit","Batch"] = 0
meta_file = gsub(" ","",meta_file)

# display column names and how many unique elements in each, of meta_file
ucol = col_probe(meta_file)
meta_file = meta_file[,-ucol$u1]

# trim control samples
meta_file = meta_file[meta_file[, "Sample Type"]=="Normal",]
meta_file = meta_file[, !colnames(meta_file)%in%c("Sample Type",
                                                  "dishQC",
                                                  "Cluster CR (Axiom inlier cluster call rate)",
                                                  "Array",
                                                  "Replicate Info",
                                                  "Number",
                                                  "Random order") &
                        !grepl("[/]|Plate|Het|Rate|Concordance|Reproducibility|Filename",colnames(meta_file))]
colnames(meta_file) = c("well","gender","centre","sample","response","kit","batch","race")

# save
save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"))



time_output(start)



## annotate Affymetric probes to SNP ----------------------------------
annot0 = probe_to_snp_annotation[match(rownames(gt1_calls),unlist(probe_to_snp_annotation[,"Probe Set ID"])),]

# display column names and how many unique elements in each, of meta_file
ucol = col_probe(annot0)
annot = annot0[,!colnames(annot0)%in%colnames(annot0)[ucol$u1]]

annot = annot0[,]




## TO BE CONTINUED

foreach









