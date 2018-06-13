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
feat_genotyperaw_dir = paste0(feat_dir,"_snp-file-genotyperaw")
feat_genotypeasthma_dir = paste0(feat_dir,"/snp-file-genotypeasthma")
feat_genotypegoodppl_dir = paste0(feat_dir,"/snp-file-genotypegoodppl")

## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
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
annot = annot0[,colnames(annot0)[ucol$u1]:=NULL]
ucol = col_probe(annot)

# column meanings
# Probe Set ID :  902560                          e.g. (matches with calls file rowname) "AFFX-SP-000001"
# Affy SNP ID :  888799                           e.g. "Affx-2643493"
# dbSNP RS ID :  852861                           e.g. "rs10466213" or "---"
# dbSNP Loctype :  6;  2 1 3 5 6 4                e.g. 1: replace + possibley insertion before the location on the subject sequence. 
#                                                      2: replace (True SNP)
#                                                      3: delete + possible delete before the location
#                                                      4: range insertion i.e. part of flank surrounding SNP is replaced with longer seq
#                                                      5: range replace i.e.   part of flank surrounding SNP is replaced with same len seq
#                                                      6: range delete i.e.    part of flank surrounding SNP is replaced with shorter seq
# Chromosome :  25                                e.g. which chromosome SNP is on i.e. 1-22 somatic, X, Y, MT (List of variations that map to the mitochondria )
# Physical Position :  886175                     e.g. "123096468"
# Position End :  886175                          e.g. "123096468"
# ChrX pseudo-autosomal region 1 :  2;  0 1       e.g.  homologous sequences of nucleotides on the X and Y chromosomes
#                                                       Genetic recombination (occurring during sexual reproduction) 
#                                                       is known to be limited only to the pseudoautosomal regions (PAR1 and PAR2) 
#                                                       of the X and Y chromosomes
#                                                       note: all SNPs here on the X,Y are pseudo-autosomal
# Cytoband :  256                                 e.g. Cytoband name to draw chromosome ideograms
# Flank :  888799                                 e.g. sequence around SNP
# Allele A :  179                                 e.g. ...
# Allele B :  4753                                e.g. ...
# Ref Allele :  4184                              e.g. ...
# Alt Allele :  1185                              e.g. ...
# Associated Gene :  571387                       e.g. "ENST00000429809 // downstream // 150988 // Hs.385516 // LINC01153 // 101927889 // long intergenic non-protein coding RNA 1153 ///
#                                                      separated by ///; ensemble gene id // downstream // ? // // transcript: LINC long introgenic non-protein coding rna // ? // trasncript: description
# Genetic Map :  866990                           e.g. 
# Microsatellite :  888169                        e.g. Microsatellites are di-, tri-, or tetra nucleotide tandem repeats in DNA sequences
#                                                      D10S294 (locus) // downstream // 145255 /// D10S1679 // upstream // 37026
# Allele Frequencies :  356858                    e.g. 0.3882 // 0.6118 // CEU /// 0.5464 // 0.4536 // CHB /// 0.5899 // 0.4101 // JPT /// 0.5227 // 0.4773 // YRI
#                                                      CEU (Utah residents with Northern and Western European ancestry from the CEPH collection)
#                                                      CHB (Han Chinese in Beijing, China)
#                                                      JPT (Japanese in Tokyo, Japan)
#                                                      YRI (Yoruba in Ibadan, Nigeria)
# Heterozygous Allele Frequencies :  318462       e.g. similar to above
# Number of individuals :  134                    e.g. hapmap is the organization who wants to develop a haplotype map of the human genome
# In Hapmap :  2;  YES ---                        e.g. 
# Strand Versus dbSNP :  3;  same --- reverse     e.g. 
# Probe Count :  10;  2 1 4 3 10 8 16 5 12 6      e.g. 
# ChrX pseudo-autosomal region 2 :  2;  0 1       e.g. 
# Minor Allele :  2065                            e.g. 
# Minor Allele Frequency :  341762                e.g. 
# OMIM :  26701                                   e.g. 176804 // {Asthma, aspirin-induced, susceptibility to} // 208550 // frameshift /// 176804 // {Asthma, aspirin-induced, susceptibility to} // 208550 // intron
# Ordered Alleles :  5015                         e.g. 
# ClinVar VariantID :  27858                      e.g. 
# ClinVar RSID :  27732                           e.g. 
# ClinVar ClinicalSignificance :  282             e.g. "---"        "PATHOGENIC"
# ClinVar GeneSymbol :  3121                      e.g. 
# ClinVar Traits :  5776                          e.g.  [1] "---"; [2] "PLATELET-ACTIVATING FACTOR ACETYLHYDROLASE DEFICIENCY"; [3] "MENTAL RETARDATION, AUTOSOMAL RECESSIVE 51" 
# ClinVar OMIM Gene :  11567                      e.g. 
# ClinVar OMIM Phenotype :  2289                  e.g. "---"    "614278" "600807"
# ClinVar OMIM Description :  2438                e.g. [1] "---"; [2] "Platelet-activating factor acetylhydrolase deficiency"; [3] "{Asthma, susceptibility to}"
# ClinVar MIM :  2932                             e.g. "---"    "601690" "605238"
# EBI PUBMEDID :  3266                            e.g. 
# EBI DISEASE/TRAIT :  2574                       e.g. "Psoriasis // Psoriasis"
# EBI MAPPED GENE(S) :  11274                     e.g. "LOC105373896 - EPHA4" "PLA2G7 // PLA2G7"  "LOC105379121 - TSLP // LOC105379121 - TSLP // LOC105379121 - TSLP" "SCGB1A1, LOC102723765"   
# EBI SNPS :  17087                               e.g. 
# EBI SNP_ID_CURRENT :  17026                     e.g. 
# EBI MAPPED_TRAIT :  2376                        e.g. "pulmonary function measurement, forced expiratory volume, asthma // body height // prostate carcinoma // body height // body height"
# EBI MAPPED_TRAIT_URI :  2376                    e.g. 
## TO BE CONTINUED

if (sum(annot[,"Physical Position"]!=annot[,"Position End"])==0) annot[,c("Position End"):=NULL]
colnames(annot)[1:6] = c("probe","affySNP","dbSNP","dbSNPloctype","chromosome","pos_phys")
save(annot, file=paste0(meta_col_dir,".Rdata"))

# get probes/SNP with asthma affiliated gene nearby?
start1 = Sys.time()
asthmarows = apply(annot, 1, function(x) any(grepl("asthma",paste(x,collapse=""),ignore.case=T)))
sum(asthmarows)
time_output(start1)

annot_asthma_id = annot$probe[asthmarows] #save indices
save(annot_asthma_id, file=paste0(meta_col_dir,"_id_asthma.Rdata"))





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

goodpplcols = unique(match(meta_file_extra$NAME, meta_file2$sample))
goodpplcols = sort(goodpplcols[!is.na(goodpplcols)])
goodpplcols_ind = as.character(meta_file2[goodpplcols,"fileName"])
save(goodpplcols_ind, file=paste0(meta_file_dir,"_id_goodppl.Rdata"))


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






