## input: features & meta_file
## output: p vaue features testing significance of ER DR correlation
## aya43@sfu.ca
## created 20180524
## last modified 20180524



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result/genotype")


asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col",asthma)

feat_dir = paste0(result_dir,"/feat")
feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype",asthma)


## output directory



## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
library(data.table)
library(MatrixEQTL)
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
id_col = "well"
class_col = "response"
categorical = T # is class column categorical?
confound_col = c("sex","centre","batch","race") #"centre"



start = Sys.time()







meta_file = get(load(paste0(meta_file_dir,".Rdata")))
meta_col = as.data.frame(get(load(paste0(meta_col_dir,".Rdata"))))
m = get(load(paste0(feat_genotype_dir,".Rdata")))

df = data.frame(class=meta_file[,id_col])
loop_ind = loop_ind_f(1:nrow(m),no_cores)
result = foreach (i=loop_ind) {
  df$val
}









time_output(start)












