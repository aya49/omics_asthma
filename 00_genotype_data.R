## input: raw genotype and meta files from lab and affymetrix axiom
## output: well formatted data and meta files
## aya43@sfu.ca
## created 20180509
## last modified 20180509



## root directory
root = "~/projects/asthma"
setwd(root)

dir.create(paste0(root, "/result"), showWarnings=F)
result_dir = paste0(root, "/result/genotype"); dir.create(result_dir, showWarnings=F)



## input directory
data_dir = "data"
genotype_dir = paste0(data_dir, "/genotype")
gt1_call_dir = paste0(genotype_dir, "/AxiomGT1.calls.txt")

meta_snp_temp_dir = paste0(genotype_dir, "/Axiom_PMRA.na35.annot.csv")

gt1_meta_dir = paste0(genotype_dir, "/additional_sample_data.txt")
meta_file_temp_dir = paste0(genotype_dir, "/meta_file_temp.csv")
meta_file_temp2_dir = paste0(genotype_dir, "/asthmaDemo_allsite.csv")
meta_file_extra_dir = paste0(data_dir, "/RNAseq/asthmaDemo_allsite.csv")



## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")
# meta_colasthma_dir = paste0(meta_dir,"/colasthma")

feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_genotyperaw_dir = paste0(feat_dir,"_snp-file-featraw")

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




start = Sys.time()

## load data --------------------------------------
gt1_calls = read.table(gt1_call_dir, sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1)
# gt1_calls = as.matrix(fread(gt1_call_dir))

meta_snp_temp = fread(meta_snp_temp_dir)

gt1_meta = fread(gt1_meta_dir, data.table=F)
meta_file_temp = fread(meta_file_temp_dir, data.table=F)
# meta_file_temp2 = fread(meta_file_temp2_dir, data.table=F)
meta_file_extra = read.csv(meta_file_extra_dir)









## save meta_col ----------------------------------
annot0 = meta_snp_temp[match(rownames(gt1_calls),unlist(meta_snp_temp[,"Probe Set ID"])),]

# display column names and how many unique elements in each, delete those with only 1 unique element
ucol = col_probe(annot0)
annot = as.data.frame(annot0[,colnames(annot0)[ucol$u1]:=NULL])
ucol = col_probe(annot)


if (sum(annot[,"Physical Position"]!=annot[,"Position End"])==0) 
  annot[,!colnames(annot)%in%c("Position End")]
# levels(annot$chromosome) = paste("chr", c(1:22, "X", "Y", "M"), sep="") #convert to bioconductor format
gwrngs.emd = as.data.frame(get(data(gwrngs38)))
# risk.alleles = gsub("[^\\-]*-([ATCG?])", "\\1", dm$Strongest.SNP.Risk.Allele)

#find asthma genes
# dm$Link[grepl("asthma",dm$Disease.Trait,ignore.case=T) & grepl("European", dm$Initial.Sample.Size) & dm$Initial.Sample.Size=="6,685 European ancestry cases, 14,091 European ancestry controls"]
# dm$Link[grepl("asthma",dm$Disease.Trait,ignore.case=T) & grepl("European", dm$Initial.Sample.Size) & dm$Initial.Sample.Size=="12,475 European ancestry cases, 19,967 European ancestry controls"]
# annot = merge(annot, gwrngs.emd, by.x="dbSNP", by.y="SNPs")
annot = cbind(annot, gwrngs.emd[match(annot[,"dbSNP"],gwrngs.emd[,"SNPs"]),])

# compare with NHS gwas studies
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene #transcriptDb; behind the scenes, everythign is SQLite
tx.by.gene = transcriptsBy(txdb, "gene") #list names are Entrez gene ID's
# columns(org.Hs.eg.db)
# # keys: APOE gene
# select(org.Hs.eg.db, keys="APOE", columns=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="SYMBOL") #keytypes()
# # look up Gene ID
# tx.by.gene["348"]
# apoe.i <- findOverlaps(tx.by.gene["348"], my.snps) #RangesMatching class; if don't give chr name, warning sequence names don't match




colnames(annot)[1:6] = c("id","affySNP","dbSNP","dbSNPloctype","chromosome","pos_phys")
save(annot, file=paste0(meta_col_dir,".Rdata"))






## save meta_file ----------------------------------


wells = gsub(".CEL","",sapply(str_split(colnames(gt1_calls),"_"), function(x) x[length(x)]))

meta_file = cbind(meta_file_temp[, !colnames(meta_file_temp)%in%c("Inferred Gender","Random order", "Number", 
                                                                  "1000 Genome Concordance", "Reproducibility", "Replicate Info", 
                                                                  "Cluster CR (Axiom inlier cluster call rate)", "dishQC")], 
                  gt1_meta[match(meta_file_temp[,"Sample ID"],gt1_meta[,"Name"]), !colnames(gt1_meta)%in%c("Name")]
                  )
meta_file[meta_file[,"Kit"]=="Mini kit","Batch"] = 0
meta_file = gsub(" ","",as.matrix(meta_file))

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
colnames(meta_file) = c("fileName","sex","centre","sample","response","kit","batch","race")
# meta_file = meta_file[,-c("kit")] #batch = 0 is a minikit; else it's a DNeasy

meta_file0 = as.data.frame(meta_file)
roworder = match(wells,meta_file0[,"fileName"])
meta_file1 = meta_file0[roworder[!is.na(roworder)],]

meta_file_extra$NAME[grepl("WRF",meta_file_extra$NAME)] = "WRF"
meta_file_extra1 = meta_file_extra[!duplicated(meta_file_extra$NAME),]

meta_file2 = cbind(meta_file1, meta_file_extra1[match(meta_file1$sample, meta_file_extra1$NAME),c(10:12,15)])
colnames(meta_file2)[9:12] = c("weight","height","age","allergen")
meta_file2[meta_file2=="" | meta_file2=="NA"] = NA
meta_file2$bmi = meta_file2$weight/(meta_file2$height^2)



# save
save(meta_file2, file=paste0(meta_file_dir,".Rdata"))
write.csv(meta_file2, file=paste0(meta_file_dir,".csv"))














## save matrix ---------------------------------------
colnames(gt1_calls) = wells
gt1_calls1 = gt1_calls[,!wells%in%setdiff(wells,meta_file[,"fileName"])]
gt1_calls1[gt1_calls1==-1] = NA
save(gt1_calls1, file=paste0(feat_genotyperaw_dir,".Rdata"))

# gt1_calls_asthma = gt1_calls1[asthmarows,]
# save(gt1_calls_asthma, file=paste0(feat_genotypeasthma_dir,".Rdata"))
# 
# gt1_calls_goodppl = gt1_calls1[,goodpplcols]
# save(gt1_calls_goodppl, file=paste0(feat_genotypegoodppl_dir,".Rdata"))








time_output(start)






