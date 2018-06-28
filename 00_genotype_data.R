## input: raw genotype and meta files from lab and affymetrix axiom
## output: well formatted data and meta files
## aya43@sfu.ca
## created 20180509
## last modified 20180509



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)
type = "genotype"


## input directory
data_dir = "data"
genotype_dir = paste0(data_dir, "/",type)
gt1_call_dir = paste0(genotype_dir, "/AxiomGT1.calls.txt")

meta_snp_temp_dir = paste0(genotype_dir, "/Axiom_PMRA.na35.annot.csv")

gt1_meta_dir = paste0(genotype_dir, "/additional_sample_data.txt")
meta_file_temp_dir = paste0(genotype_dir, "/meta_file_temp.csv")
meta_file_temp2_dir = paste0(genotype_dir, "/asthmaDemo_allsite.csv")
meta_file_extra_dir = paste0(data_dir, "/RNAseq/asthmaDemo_allsite.xlsx")

meta_snp_idrod_dir = paste0(genotype_dir, "/rod/rod.csv")


## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col-",type)

feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_genotype_dir = paste0(feat_dir,"/",type)

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
libr("gdata") #read xls
libr("Matrix")



## options
options(stringsAsFactors=F)

no_cores = 15#detectCores()-3
registerDoMC(no_cores)

writecsv = T #write results as csv on top of Rdata?
# for trimming meta and matrix
good_col = 3 #each gene must have > good_col samples of each genotype; else delete
good_col_nap = .5 #each gene must have > good_col_na x 100% samples with data

id_col = "id"



start = Sys.time()


## load matrix ---------------------------------------

gt1_calls = read.table(gt1_call_dir, sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1)





## save meta_file ----------------------------------

# load files
gt1_meta = fread(gt1_meta_dir, data.table=F)
meta_file_temp = fread(meta_file_temp_dir, data.table=F)
# meta_file_temp2 = fread(meta_file_temp2_dir, data.table=F)
meta_file_extra = read.xls(meta_file_extra_dir, check.names=T)

# extract well = unique identifiers for each sample
wells = gsub(".CEL","",sapply(str_split(colnames(gt1_calls),"_"), function(x) x[length(x)]))

# keep some useful? columns
meta_file0 = cbind(meta_file_temp[, !colnames(meta_file_temp)%in%c("Inferred Gender","Random order", "Number", 
                                                                  "1000 Genome Concordance", "Reproducibility", "Replicate Info", 
                                                                  "Cluster CR (Axiom inlier cluster call rate)", "dishQC")], 
                  gt1_meta[match(meta_file_temp[,"Sample ID"],gt1_meta[,"Name"]), !colnames(gt1_meta)%in%c("Name")]
                  )

# adjust values
meta_file0[meta_file0[,"Kit"]=="Mini kit","Batch"] = 0
meta_file0 = gsub(" ","",as.matrix(meta_file0))

# display column names and how many unique elements in each, of meta_file; delete cols with only one unique element
ucol = col_probe(meta_file0)
meta_file0 = meta_file0[,-ucol$u1]

# remove control samples
meta_file0 = meta_file0[meta_file0[, "Sample Type"]=="Normal",]

# delete more columns
meta_file0 = meta_file0[, !colnames(meta_file0)%in%c("Sample Type",
                                                  "dishQC",
                                                  "Cluster CR (Axiom inlier cluster call rate)",
                                                  "Array",
                                                  "Replicate Info",
                                                  "Number",
                                                  "Random order") &
                        !grepl("[/]|Plate|Het|Rate|Concordance|Reproducibility|Filename",colnames(meta_file0))]

# rename columns
colnames(meta_file0) = c(paste0("filename_",type),"sex",paste0("centre_",type),id_col,"response",paste0("kit_",type),paste0("batch_",type),"race")
# meta_file0 = meta_file0[,-c("kit")] #batch = 0 is a minikit; else it's a DNeasy

# order rows according to genotype matrix
meta_file0 = as.data.frame(meta_file0)
roworder = match(wells,meta_file0[,paste0("filename_",type)])
meta_file1 = meta_file0[roworder[!is.na(roworder)],]

# merge meta_files
meta_file_extra$NAME[grepl("WRF",meta_file_extra$NAME)] = "WRF"
meta_file_extra1 = meta_file_extra[!duplicated(meta_file_extra$NAME),]

meta_file_extra2 = merge(meta_file1, meta_file_extra1, by.x="id", by.y="NAME", all=T)

# reconcile duplicate values
meta_file_extra2$SITE[meta_file_extra2$SITE=="MAC"] = "McMaster"
meta_file_extra2$SITE = tolower(meta_file_extra2$SITE)
meta_file_extra2$centre_genotype = tolower(meta_file_extra2$centre_genotype)
meta_file_extra2$SITE[is.na(meta_file_extra2$SITE)] = 
  meta_file_extra2$centre_genotype[is.na(meta_file_extra2$SITE)]

meta_file_extra2$response[is.na(meta_file_extra2$response)] = meta_file_extra2$CorrectResponse[is.na(meta_file_extra2$response)]
meta_file_extra2$response[is.na(meta_file_extra2$response)] = meta_file_extra2$Mac_Response[is.na(meta_file_extra2$response)]

meta_file_extra2$race[meta_file_extra2$race==""] = "Unknown"

meta_file_extra2$sex[is.na(meta_file_extra2$sex)] = meta_file_extra2$SEX[is.na(meta_file_extra2$sex)]
meta_file_extra2$sex[meta_file_extra2$sex=="Female"] = "F"
meta_file_extra2$sex[meta_file_extra2$sex=="Male"] = "M"
meta_file_extra2$sex[meta_file_extra2$sex==""] = NA


# delete more columns
meta_file2 = meta_file_extra2[,c(1:3,5,7:11,13,17:21,23:35,42,48,
                                 53:65)] # data availability

colnames(meta_file2)[c(7:16,30)] = c("sponsor","centre","drug","date","weight","height","age","blfev","prfev","allergen","cohort")
meta_file2[meta_file2=="" | meta_file2=="NA"] = NA
meta_file2$bmi = meta_file2$weight/(meta_file2$height^2)


# save
save(meta_file2, file=paste0(meta_file_dir,".raw.Rdata"))
if (writecsv) write.csv(meta_file2, file=paste0(meta_file_dir,".raw.csv"))

# meta_file = meta_file2







## save meta_col ----------------------------------
meta_snp_temp = fread(meta_snp_temp_dir)

annot0 = meta_snp_temp[match(rownames(gt1_calls),unlist(meta_snp_temp[,"Probe Set ID"])),]

# display column names and how many unique elements in each, delete those with only 1 unique element
ucol = col_probe(annot0)
annot = as.data.frame(annot0[,colnames(annot0)[ucol$u1]:=NULL])
ucol = col_probe(annot)

if (sum(annot[,"Physical Position"]!=annot[,"Position End"])==0) 
  annot[,!colnames(annot)%in%c("Position End")]
# levels(annot$chromosome) = paste("chr", c(1:22, "X", "Y", "M"), sep="") #convert to bioconductor format
gwrngs.emd = as.data.frame(get(data(gwrngs38)))
colnames(gwrngs.emd) = paste0("gwrngs_",colnames(gwrngs.emd))
# risk.alleles = gsub("[^\\-]*-([ATCG?])", "\\1", dm$Strongest.SNP.Risk.Allele)

#find asthma genes
# dm$Link[grepl("asthma",dm$Disease.Trait,ignore.case=T) & grepl("European", dm$Initial.Sample.Size) & dm$Initial.Sample.Size=="6,685 European ancestry cases, 14,091 European ancestry controls"]
# dm$Link[grepl("asthma",dm$Disease.Trait,ignore.case=T) & grepl("European", dm$Initial.Sample.Size) & dm$Initial.Sample.Size=="12,475 European ancestry cases, 19,967 European ancestry controls"]
# annot = merge(annot, gwrngs.emd, by.x="dbSNP", by.y="SNPs")
gw = gwrngs.emd[match(annot$`dbSNP RS ID`, gwrngs.emd$gwrngs_SNPs),]
# gw_ind = c(1:nrow(gw))
# gw_ind[!grepl("European",gw$gwrngs_Initial.Sample.Size)] = NA
annot = cbind(annot, gw)

# compare with NHS gwas studies
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene #transcriptDb; behind the scenes, everythign is SQLite
tx.by.gene = transcriptsBy(txdb, "gene") #list names are Entrez gene ID's
# columns(org.Hs.eg.db)
# # keys: APOE gene
# select(org.Hs.eg.db, keys="APOE", columns=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="SYMBOL") #keytypes()
# # look up Gene ID
# tx.by.gene["348"]
# apoe.i <- findOverlaps(tx.by.gene["348"], my.snps) #RangesMatching class; if don't give chr name, warning sequence names don't match

colnames(annot)[1:6] = c(id_col,"affySNP","dbSNP","dbSNPloctype","chromosome","pos_phys")
save(annot, file=paste0(meta_col_dir,".raw.Rdata"))
if (writecsv) write.csv(annot, file=paste0(meta_col_dir,".raw.csv"))







## save matrix ---------------------------------------

colnames(gt1_calls) = wells
gt1_calls1 = t(gt1_calls[,wells%in%meta_file2[,paste0("filename_",type)]])
gt1_calls1[gt1_calls1==-1] = NA
save(gt1_calls1, file=paste0(feat_genotype_dir,".raw.Rdata"))
# if (writecsv) write.csv(gt1_calls1, file=paste0(feat_genotype_dir,".raw.csv"))

# gt1_calls_asthma = gt1_calls1[asthmarows,]
# save(gt1_calls_asthma, file=paste0(feat_genotypeasthma_dir,".Rdata"))
# 
# gt1_calls_goodppl = gt1_calls1[,goodpplcols]
# save(gt1_calls_goodppl, file=paste0(feat_genotypegoodppl_dir,".Rdata"))



time_output(start)
start = Sys.time()






## trim ---------------------------------------------------

# get rid of column with less than a certain number of NA
good_col_na = good_col_nap * nrow(gt1_calls1)
col_ind = apply(gt1_calls1, 2, function(x) {
  a = table(x) 
  sum(!is.na(x)) >= good_col_na & 
    min(a)>good_col & length(a)>1
})
row_ind = meta_file2$batch_genotype!="0" # make patient name the unique id

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

meta_file = meta_file2[row_ind,!colnames(meta_file2)%in%c("kit_genotype")]
meta_col = annot[col_ind,]
m = gt1_calls1[rownames(gt1_calls1)%in%meta_file[,paste0("filename_",type)],col_ind]
rownames(m) = meta_file[match(rownames(m),meta_file[,paste0("filename_",type)]),id_col]

save(meta_col, file=paste0(meta_col_dir,".Rdata"))
if (writecsv) write.csv(meta_col, file=paste0(meta_col_dir,".csv"))

save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"))

save(m, file=paste0(feat_genotype_dir,".Rdata"))
if (writecsv) write.csv(m, file=paste0(feat_genotype_dir,".csv"))





## save some indices ----------------------------------------------------------------

# get probes/SNP with asthma affiliated gene nearby?
start1 = Sys.time()
asthmarows = apply(meta_col, 1, function(x) any(grepl("asthma",paste(x,collapse=""),ignore.case=T))) & !duplicated(meta_col$dbSNP)
sum(asthmarows)
time_output(start1)

meta_col_asthma_id = meta_col[asthmarows, id_col] #save indices
save(meta_col_asthma_id, file=paste0(meta_col_dir,"_id_asthma-ebiclinvaromim.Rdata"))

asthmarows_gwrns = 
  grepl("asthma",meta_col$gwrngs_Disease.Trait,ignore.case=T) & 
  grepl("European", meta_col$gwrngs_Initial.Sample.Size)
meta_col_asthma_id_gwrns = meta_col[asthmarows_gwrns, id_col] #save indices
save(meta_col_asthma_id_gwrns, file=paste0(meta_col_dir,"_id_asthma-gwrns.Rdata"))

asthmarows_rod0 = read.csv(meta_snp_idrod_dir)
meta_col_asthma_id_rod = meta_col[meta_col$dbSNP%in%asthmarows_rod0[,"SNP"], id_col]
save(meta_col_asthma_id_rod, file=paste0(meta_col_dir,"_id_asthma-rod.Rdata"))

gwrngs.emd$gwrngs_Disease.Trait[gwrngs.emd$gwrngs_SNPs%in%asthmarows_rod0[,"SNP"]]
meta_col$gwrngs_Disease.Trait[meta_col$gwrngs_SNPs%in%asthmarows_rod0[,"SNP"]]

meta_col_asthma_id_st = meta_col[meta_col$dbSNP%in%c("rs993076","rs1800777") | 
                                   apply(meta_col, 1, function(x) 
                                     any(grepl("PLPP3|CETP|farp1|tlr1|bcl10|mal1|card9|card11|pcsk9|serping1",paste(x,collapse=""),ignore.case=T)))
                                 , id_col]
save(meta_col_asthma_id_st, file=paste0(meta_col_dir,"_id_asthma-st.Rdata"))



# get good files
goodpplcols_ind = meta_file[!is.na(meta_file[,"bmi"]),id_col]
save(goodpplcols_ind, file=paste0(meta_file_dir,"_id_goodppl.Rdata"))








