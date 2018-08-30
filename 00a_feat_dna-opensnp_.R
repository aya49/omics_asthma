## input: raw dna and meta files from lab and affymetrix axiom
## output: well formatted data and meta files
## aya43@sfu.ca
## created 20180615
## last modified 20180615



## root directory
root = "~/projects/asthma"
setwd(root)

dir.create(paste0(root, "/result"), showWarnings=F)
result_dir = paste0(root, "/result/dna"); dir.create(result_dir, showWarnings=F)



## input directory
data_dir = "data"
dna_dir = paste0(data_dir, "/dna")
opensnp_meta_dir = paste0(dna_dir, "/opensnp/phenotypes_201805150935.csv")
opensnp_zip_dir = paste0(dna_dir, "/opensnp/opensnp_datadump.current.zip")

hapmap_call_dir = paste0(dna_dir, "/hapmap_dna/Nsp_HapMap270.brlmm.call")
hapmap_conf_dir = paste0(dna_dir, "/hapmap_dna/Nsp_HapMap270.brlmm.conf")
hapmap_call0_dir = paste0(dna_dir, "/hapmap_dna/Sty_HapMap270.brlmm.call")
hapmap_conf0_dir = paste0(dna_dir, "/hapmap_dna/Sty_HapMap270.brlmm.conf")
hapmap_all_dir = paste0(dna_dir, "/hapmap_dna/all_CEU.txt")


## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")
# meta_colasthma_dir = paste0(meta_dir,"/colasthma")

feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_dnaraw_dir = paste0(feat_dir,"_snp-file-dnaraw")
feat_dnaasthma_dir = paste0(feat_dir,"/snp-file-dnaasthma")
feat_dnagoodppl_dir = paste0(feat_dir,"/snp-file-dnagoodppl")

## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("TxDb.Hsapiens.UCSC.hg19.knownGene")
libr("org.Hs.eg.db") #org.* annotation packages; can forge own and interact with using library("AnnotationDbi")
libr("gwascat") # interface to the [NHGRI's](http://www.genome.gov/) database of gwas

libr("data.table")
libr("entropy")
libr("foreach")
libr("doMC")
libr("stringr")
libr("Matrix")



## options
options(stringsAsFactors=F)

no_cores = 15#detectCores()-3
registerDoMC(no_cores)

writecsv = T #write results as csv on top of Rdata?

good_col = 3
key_word = "asthma"



start = Sys.time()

## load data --------------------------------------

meta_file_opensnp = read.csv(opensnp_meta_dir,header=T)

acol = grep(key_word,colnames(meta_file_opensnp),ignore.case=T)
acoli = grep(key_word,str_split(colnames(meta_file_opensnp)[acol],"[.]")[[1]],ignore.case=T)
acolval = sapply(str_split(meta_file_opensnp[,acol],"[;]"), function(x) x[acoli])

read.table(unz("temp.zip", ""))




hapmap_call = fread(hapmap_call_dir)
hapmap_conf = fread(hapmap_conf_dir)
hapmap_call0 = fread(hapmap_call0_dir)
hapmap_conf0 = fread(hapmap_conf0_dir)
hapmap_all = fread(hapmap_all_dir)






