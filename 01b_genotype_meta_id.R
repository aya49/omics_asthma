## input: SNP references & file information
## output: indices to test the feature matrix on
## aya43@sfu.ca
## created 20180614



## root directory
root = "~/projects/asthma"
setwd(root)

dir.create(paste0(root, "/result"), showWarnings=F)
result_dir = paste0(root, "/result/genotype"); dir.create(result_dir, showWarnings=F)


## input directory
data_dir = "data"
genotype_dir = paste0(data_dir, "/genotype")
meta_snp_idrod_dir = paste0(genotype_dir, "/rod/rod.csv")

meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

## output directory



## libraries
source("code/_func.R")
libr("TxDb.Hsapiens.UCSC.hg19.knownGene")
libr("org.Hs.eg.db") #org.* meta_colation packages; can forge own and interact with using library("meta_colationDbi")
libr("gwascat") # interface to the [NHGRI's](http://www.genome.gov/) database of gwas
libr("foreach")
libr("doMC")
libr("stringr")


## options
options(stringsAsFactors=F)

no_cores = 15#detectCores()-3
registerDoMC(no_cores)

cid_col = "id"
id_col = "fileName"




start = Sys.time()


meta_col = get(load(paste0(meta_col_dir,".Rdata")))
meta_file = get(load(paste0(meta_file_dir,".Rdata")))


# get probes/SNP with asthma affiliated gene nearby?
start1 = Sys.time()
asthmarows = apply(meta_col, 1, function(x) any(grepl("asthma",paste(x,collapse=""),ignore.case=T))) & !duplicated(meta_col$dbSNP)
sum(asthmarows)
time_output(start1)

meta_col_asthma_id = meta_col[asthmarows, cid_col] #save indices
save(meta_col_asthma_id, file=paste0(meta_col_dir,"_id_asthma.Rdata"))

asthmarows_gwrns = grepl("asthma",meta_col$Disease.Trait,ignore.case=T) & grepl("European", meta_col$Initial.Sample.Size)
meta_col_asthma_id_gwrns = meta_col[asthmarows_gwrns, cid_col] #save indices
save(meta_col_asthma_id_gwrns, file=paste0(meta_col_dir,"_id_asthma-gwrns.Rdata"))

asthmarows_rod0 = read.csv(meta_snp_idrod_dir)
meta_col_asthma_id_rod = meta_col[meta_col$dbSNP%in%asthmarows_rod0[,"SNP"], cid_col]
save(meta_col_asthma_id_rod, file=paste0(meta_col_dir,"_id_asthma-rod.Rdata"))

meta_col_asthma_id_st = meta_col[meta_col$dbSNP%in%c("rs993076","rs1800777") | 
                                apply(meta_col, 1, function(x) 
                                  any(grepl("PLPP3|CETP|farp1|tlr1|bcl10|mal1|card9|card11|pcsk9|serping1",paste(x,collapse=""),ignore.case=T)))
                                , cid_col]
save(meta_col_asthma_id_st, file=paste0(meta_col_dir,"_id_asthma-st.Rdata"))



# get good files
goodpplcols_ind = meta_file[!is.na(meta_file[,"bmi"]),id_col]
save(goodpplcols_ind, file=paste0(meta_file_dir,"_id_goodppl.Rdata"))
