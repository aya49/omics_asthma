## input: raw genotype and meta files
## output: well formatted meta files
## aya43@sfu.ca
## created 20180509
## last modified 20180509



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)


## input directory
data_dir = paste0(root,"/data")
genotype_dir = paste0(data_dir, "/genotype")
gt1_call_dir = paste0(genotype_dir, "/AxiomGT1.calls.txt")

meta_snp_temp_dir = paste0(genotype_dir, "/Axiom_PMRA.na35.annot.csv")
rnaseq_dir = paste0(root,"/data/RNAseq")

meta_file_temp2_dir = paste0(rnaseq_dir,"/rnaseq_demo.Rdata")
meta_file_rnaseq_dir = paste0(root,"/data/asthmaDemo_allsite_withSampleInfo_DH_v5.csv")
meta_file_rnapc_dir = paste0(data_dir,"/RNAelements/data/allsitesDemo_clean.rds")

gt1_meta_dir = paste0(genotype_dir, "/additional_sample_data.txt")
meta_file_temp_dir = paste0(genotype_dir, "/meta_file_temp.csv")
meta_file_temp3_dir = paste0(data_dir, "/allsitesDemo.csv")
meta_file_extra_dir = paste0(data_dir, "/RNAseq/asthmaDemo_allsite.xlsx")
meta_file_data_dir = paste0(root,"/data/asthmaDemo_allsite_withSampleInfo_DH_v5.csv")


## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")

## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr(c("TxDb.Hsapiens.UCSC.hg19.knownGene", 
       "org.Hs.eg.db", #org.* annotation packages; can forge own and interact with using library("AnnotationDbi")
       "annotables",
       "gwascat", # interface to the [NHGRI's](http://www.genome.gov/) database of gwas
       #"entropy",
       "data.table", "Matrix", "gdata", #read xls
       #"foreach", "doMC",
       "stringr", "pracma"))


## options
options(stringsAsFactors=F)

no_cores = 15#detectCores()-3
registerDoMC(no_cores)

writecsv = T #write results as csv on top of Rdata?

id_col = "id"



start = Sys.time()







## meta 1: genotype ----------------------------------

# load genotype filenames
gt1_calls_colnames = names(read.table(gt1_call_dir, sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1, nrows=1))
wells = gsub(".CEL","",sapply(str_split(gt1_calls_colnames,"_"), function(x) x[length(x)]))

# load and bind genotype affymetrix & lab meta files (same rows)
gt1_meta = fread(gt1_meta_dir, data.table=F)
gt1lab_meta = fread(meta_file_temp_dir, data.table=F)
meta_file0 = 
  cbind(gt1lab_meta, 
        gt1_meta[match(gt1lab_meta[,"Sample ID"],gt1_meta[,"Name"]), !colnames(gt1_meta)%in%c("Name")])

# display column names and how many unique elements in each, of meta_file; delete cols with only one unique element
ucol = col_probe(meta_file0)
meta_file0 = meta_file0[,-ucol$u1]

# remove control samples
meta_file0 = meta_file0[meta_file0[, "Sample Type"]=="Normal",]

# adjust values
meta_file0[meta_file0[,"Kit"]=="Mini kit","Batch"] = 0
meta_file0 = gsub(" ","",as.matrix(meta_file0))
# # get rid of batch 0; there's a copy of it in other batches
# row_ind = meta_file0[,"Batch"]!="0" & !grepl("Normal",meta_file0[,"Sample ID"],ignore.case=T) # make patient name the unique id
# meta_file = meta_file0[row_ind,]


# delete columns
meta_file0 = 
  meta_file0[, !colnames(meta_file0) %in%
               c("Sample Type", "dishQC",
                 "Cluster CR (Axiom inlier cluster call rate)",
                 "Array", "Replicate Info", "Number", "Random order",
                 "Inferred Gender","Random order", "Number", 
                 "1000 Genome Concordance", "Reproducibility", "Replicate Info", 
                 "Cluster CR (Axiom inlier cluster call rate)", "dishQC") &
               !grepl("[/]|Plate|Het|Rate|Concordance|Reproducibility|Filename",colnames(meta_file0))]

# rename columns
colnames(meta_file0) = 
  c("filename_genotype","sex","centre",id_col,
    "response","kit_genotype", "batch_genotype","race")
# meta_file0 = meta_file0[,-c("kit")] #batch = 0 is a minikit; else it's a DNeasy

# order rows according to genotype matrix
meta_file0 = as.data.frame(meta_file0)
roworder = match(wells,meta_file0[,paste0("filename_genotype")])
meta_file1 = meta_file0[roworder[!is.na(roworder)],]

# check for flippers
response_table = rowSums(table(meta_file1$id, meta_file1$response) != 0)
flippers <- names(response_table[response_table >1])
flippers_non <- names(response_table[response_table ==1])




## meta 2: general ------------------------------
meta_file_extra = read.xls(meta_file_extra_dir, check.names=T)
meta_file_extra = meta_file_extra[!grepl("Normal",meta_file_extra$NAME),]

#check if any rows missing for the pancancer data set
meta_file_extra_id = paste(meta_file_extra$NAME, meta_file_extra$AIC_YMD, meta_file_extra$Time)
meta_file_pc = readRDS(meta_file_rnapc_dir)
meta_file_pc = meta_file_pc[!grepl("Normal",meta_file_pc$NAME),] #delete controls
meta_file_pc$NAME[grepl("WRF",meta_file_pc$NAME)] = "WRF"
meta_file_pc_id = paste(meta_file_pc$NAME, meta_file_pc$AIC_YMD, meta_file_pc$Time)
sum(!meta_file_pc_id%in%meta_file_extra_id)

#bmi
meta_file_extra$bmi = as.numeric(meta_file_extra[,"Wt..Kg."])/(meta_file_extra[,"HT.cm."]^2)

# id_namedate = subject names _ date 
meta_file_extra$id_namedate = paste("uniqueID",
                                    as.numeric(factor(as.character(meta_file_extra$NAME))),
                                    as.numeric(factor(as.character(meta_file_extra$AIC_YMD))), sep="_")

# id_name = subject names 
meta_file_extra$id_name = factor(paste("ID", as.numeric(factor(as.character(meta_file_extra$NAME))), sep="_"))

# scale breathing test at diff times columns based on blfev
dropFEV1 = as.data.frame(t(100*scale(t(
  meta_file_extra[, c("BLFEV","F10L","F20L","F30L","F45L","F60L","F90L","F120L",
                      "F180L","F240L","F300L","F360L","F420L")]), 
  center=meta_file_extra$BLFEV, scale=meta_file_extra$BLFEV)))

# calculate response
ear = suppressWarnings( apply(dropFEV1, 1, function(x) {min(x[1:8], na.rm=T)}) )
ear[ear == "Inf"] = NA
lar = suppressWarnings( apply(dropFEV1, 1, function(x) {min(x[9:13], na.rm=T)}) )
lar[lar == "Inf"] = NA
meta_file_extra$EAR = round(ear, 1)
meta_file_extra$LAR = round(lar, 1)

time = gsub("L", "", substring(colnames(dropFEV1), 2))
# time = vapply(str_split(colnames(dropFEV1), "_"),'[',3, FUN.VALUE=character(1)) %>% as.numeric
time[1] = 0;
time = as.numeric(time)/60

# AUR trapz: area of a function for each patient, all and LAS time vs drops
meta_file_extra$AUC = apply(dropFEV1, 1, function(y) trapz(x=time, y=y) )
meta_file_extra$AUC_LAR = apply(dropFEV1, 1, function(y) trapz(x=time[9:13], y=y[9:13]) )

#AIS: PC20 of pre/post
ais = meta_file_extra[meta_file_extra$Time == "Pre", "PC20"] / meta_file_extra[meta_file_extra$Time == "Post", "PC20"]
names(ais) = meta_file_extra$id_namedate[meta_file_extra$Time == "Pre"]
meta_file_extra$AIS = ais

# calculate a reponse
meta_file_extra$response_calc = rep(NA, nrow(meta_file_extra))
meta_file_extra$response_calc[meta_file_extra$LAR > -15] = "ER"
meta_file_extra$response_calc[meta_file_extra$LAR <= -15] = "DR"
meta_file_extra$response_calc[as.character(meta_file_extra$Mac_Response) == "Control"] = "Control"

# check for flippers
response_table = rowSums(table(meta_file_extra$NAME, meta_file_extra$response_calc) != 0)
flippers <- names(response_table[response_table >1])
flippers_non <- names(response_table[response_table ==1])

# delete people without a response
meta_file_extra1 = meta_file_extra1[meta_file_extra1$response!="",]

# adjust id
meta_file_extra$NAME[grepl("WRF",meta_file_extra$NAME)] = "WRF"

# isolate one sample only for each subject's pre/post
u_ind = as.vector(sapply(unique(meta_file_extra$NAME), function(x) {
  n_indpre = which(meta_file_extra$NAME==x & meta_file_extra$Time=="Pre")
  n_indpost = which(meta_file_extra$NAME==x & meta_file_extra$Time=="Post")
  return(c(n_indpost[ which.max(apply(meta_file_extra[n_indpost,],1, function(y) sum(y!="" & !is.na(y))))],
           n_indpre[ which.max(apply(meta_file_extra[n_indpre,],1, function(y) sum(y!="" & !is.na(y))))]))
}))
meta_file_extra1 = meta_file_extra[u_ind,]

# reconcile duplicate values
meta_file_extra1$SITE[meta_file_extra1$SITE=="MAC"] = "McMaster"
meta_file_extra1$SITE = tolower(meta_file_extra1$SITE)
meta_file_extra1$centre = tolower(meta_file_extra1$centre)
meta_file_extra1$SITE[is.na(meta_file_extra1$SITE)] = 
  meta_file_extra1$centre[is.na(meta_file_extra1$SITE)]

meta_file_extra1$response[is.na(meta_file_extra1$response)] = meta_file_extra1$CorrectResponse[is.na(meta_file_extra1$response)]
meta_file_extra1$response[is.na(meta_file_extra1$response)] = meta_file_extra1$Mac_Response[is.na(meta_file_extra1$response)]

# meta_file_extra1$race[meta_file_extra1$race==""] = "Unknown"
# 
# meta_file_extra1$sex[is.na(meta_file_extra1$sex)] = meta_file_extra1$SEX[is.na(meta_file_extra1$sex)]
# meta_file_extra1$sex[meta_file_extra1$sex=="Female"] = "F"
# meta_file_extra1$sex[meta_file_extra1$sex=="Male"] = "M"
# meta_file_extra1$sex[meta_file_extra1$sex==""] = NA


## merge meta genotype and meta general with democlean ---------------------------------
meta_file_extra2 = cbind(meta_file1[match(meta_file_extra1$NAME,meta_file1$id),],meta_file_extra1)

match_col = match(colnames(meta_file_extra2),colnames(meta_file_pc))
match_col[c(2)] = c(9)

colnames(meta_file_extra2)[colnames(meta_file_extra2)%in%"response_calc"] = "response"
meta_file_extra2$flipper_calc = F
meta_file_extra2$flipper_calc[meta_file_extra2$id%in%flippers] = T
meta_file_extra2$race[grepl("Unknown",meta_file_extra2$race) | meta_file_extra2$race=="" | is.na(meta_file_extra2$race)] = meta_file_extra2$RACE[grepl("Unknown",meta_file_extra2$race) | meta_file_extra2$race=="" | is.na(meta_file_extra2$race)]
meta_file_extra2$race[meta_file_extra2$race=="White"] = "Caucasian"
meta_file_extra2$race = tolower(meta_file_extra2$race)
meta_file_extra2$SEX[meta_file_extra2$SEX==""] = meta_file_extra2$sex[meta_file_extra2$SEX==""]
meta_file_extra2$SEX[meta_file_extra2$SEX=="Female"] = "F"
meta_file_extra2$SEX[meta_file_extra2$SEX=="Male"] = "M"
meta_file_extra2$SEX[meta_file_extra2$SEX==""] = NA

meta_file_extra2 = meta_file_extra2[!is.na(meta_file_extra2$response),]


## merge meta genotype & meta general with 3: meta rnaseq filnames ---------------------------------
meta_file1_temp2 = get(load(meta_file_temp2_dir)) # contains rnaseq ids
meta_file_data = read.csv(meta_file_rnaseq_dir)
meta_file_data$UniqueID[meta_file_data$rnaseq=="Y"]

meta_file1_temp2 = data.frame(lapply(meta_file1_temp2, as.character), stringsAsFactors=FALSE)
meta_file1_temp2[grepl("WRF",meta_file1_temp2[,"NAME"]),"NAME"] = "WRF"
meta_file_extra2$filename_rnaseq = meta_file1_temp2$Read.Set.Id[
  match(paste(meta_file_extra2$id,meta_file_extra2$Time),
        paste(meta_file1_temp2$NAME,meta_file1_temp2$Time))]


## delete columns from merged meta_file ---------------------------------------
# (save raw version with both response=pre/post, another one with just pre)
meta_file2 = meta_file_extra2[,c(1,3,13,6:9,11,14,17:22,24:36,43,44,49,54:66,70,73:78,80,51,81)] # data availability
colnames(meta_file2)[c(3,7:13,16,30,31,51:53)] = 
  c("id","sponsor","drug","date","sex","weight","height",
    "age","allergen","time","cohort","response","flipper","id_mstlst")

# isolate meta table with id as subject name
meta_file = meta_file2[meta_file2$time=="Pre", c(1:13,16,30:45,51:54)]
colnames(meta_file)[colnames(meta_file)%in%"filename_rnaseq"] = "filename_rnaseq.pre"
meta_file$filename_rnaseq.post = meta_file2$filename_rnaseq[grep("Post",meta_file2$time)[
  match(meta_file$id,meta_file2$id[meta_file2$time=="Post"])]]

## save --------------------------------------------
save(meta_file2, file=paste0(meta_file_dir,".raw.Rdata"))
if (writecsv) write.csv(meta_file2, file=paste0(meta_file_dir,".raw.csv"))
save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"))

goodppl = meta_file$id[!meta_file$flipper]
save(goodppl, file=paste0(meta_file_dir,"_id_goodppl.Rdata"))

