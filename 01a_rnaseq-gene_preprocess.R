## input: RSEM output .isoforms/genes.results files
## aya43@sfu.ca
## created 20180509
## last modified 20180509

## root directory
root = "~/projects/asthma"
setwd(root)

dir.create(paste0(root, "/result"), showWarnings=F)
result_dir = paste0(root, "/result/RNAseq/gene"); dir.create(result_dir, showWarnings=F)

## input directory
data_dir = paste0(root,"/data/RNAseq")

## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_gene_dir = paste0(meta_dir,"/gene")
feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_countraw_dir = paste0(feat_dir,"/gene-file-countraw")

## libraries
library(data.table)
library(entropy)
library(foreach)
library(doMC)
library(stringr)
library(Matrix)
source(paste0(root,"/codes/_func.R"))













start = Sys.time()

