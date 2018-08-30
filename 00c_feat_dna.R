## input: raw dna and meta files from affymetrix axiom product
## output: well formatted data and meta col files
## aya43@sfu.ca
## created 20180509
## last modified 20180509


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
libr(append(pkgs(),"traseR"))

data(taSNP) #subset of below
taSNP0 = as.data.frame(taSNP)
data(taSNPLD)
taSNPLD0 = as.data.frame(taSNPLD)



## options
writecsv = T #write results as csv on top of Rdata?
# for trimming meta and matrix
type = "dna"


## start
start = Sys.time()


## load matrix & meta ---------------------------------------
gt1_calls = read.table(gt1_call_dir, sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1)

gt1meta = read.table(gt1_meta_dir, sep="\t", header=T)
gt1meta = gt1meta[match(colnames(gt1_calls),gt1meta[,1]),]
meta_file2 = get(load(paste0(meta_file_dir,".raw.Rdata")))
meta_file = get(load(paste0(meta_file_dir,".Rdata")))


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
annot = apply(annot,2,function(x) { x[x=="---"] = NA; x })

save(annot, file=paste0(meta_col_gt_dir,".raw.Rdata"))
if (writecsv) write.csv(annot, file=paste0(meta_col_gt_dir,".raw.csv"))


## save matrix ---------------------------------------
wells = str_extract(colnames(gt1_calls), "[A-Z][0-9][0-9]")
colnames(gt1_calls) = wells
gt1_calls1 = t(gt1_calls[,wells%in%meta_file2[,paste0("filename_",type)]])
gt1_calls1[gt1_calls1==-1] = NA
save(gt1_calls1, file=paste0(feat_dna_dir,".raw.Rdata"))
# if (writecsv) write.csv(gt1_calls1, file=paste0(feat_dna_dir,".raw.csv"))

# gt1_calls_asthma = gt1_calls1[asthmarows,]
# save(gt1_calls_asthma, file=paste0(feat_dnaasthma_dir,".Rdata"))
# 
# gt1_calls_goodppl = gt1_calls1[,goodpplcols]
# save(gt1_calls_goodppl, file=paste0(feat_dnagoodppl_dir,".Rdata"))



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

meta_col_f = meta_col = as.data.frame(annot[col_ind,])

start = Sys.time()
grch38_dt <- merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")
meta_col$symbols = sapply(str_extract_all( meta_col$`Associated Gene`, "ENST[0-9]+" ), function(x) {
  symb = grch38_dt$symbol[match(x,grch38_dt$enstxp)]
  symb = symb[!duplicated(symb) & !is.na(symb)]
  if (length(symb)==0) return("")
  return(paste(symb,collapse="_"))
})

# enst_split = str_extract_all(meta_col[,"Associated Gene"], "ENST[0-9]+")
# enst = unique(unlist(enst_split))
# for (enst_i in enst) {
#   symb = grch38_dt[grch38_dt$enstxp==enst_i, "symbol"]
#   if (length(symb)>0) {
#     colind_w_symb = sapply(enst_split, function(x) any(x==enst_i) | any(is.na(x)))
#     meta_col[colind_w_symb,"symbols"] = paste0(meta_col[colind_w_symb,"symbols"],"_",symb)
#   }
# }
time_output(start)

meta_col_g = meta_col

meta_col_pc = read.csv(paste0(data_dir, "/Cells_nCounter_Human_PanCancer_Immune_Profiling_Panel_Gene_List.csv"), row.names=1)
if (!"pc_traits"%in%colnames(meta_col)) {
  start = Sys.time()
  meta_col$pc_traits =  meta_col$pc_traits_class = ""
  for (i in 1:nrow(meta_col_pc)) {
    gene = meta_col_pc$id[i]
    g_ind = grepl(gene, meta_col[,"symbol"])
    meta_col$pc_traits[g_ind] = 
      paste(meta_col$pc_traits[g_ind], "_", meta_col_pc[i,"trait"])
    meta_col$pc_traits_class[g_ind] = 
      paste(meta_col$pc_traits_class[g_ind], "_", meta_col_pc[i,"trait_class"])
  }
  time_output(start)
}
meta_col_tmp = meta_col
meta_col[meta_col=="" | meta_col=="---"] = NA

meta_col = cbind(meta_col, taSNPLD0[match(meta_col$dbSNP, taSNPLD0$SNP_ID), c("Trait", "Trait_Class")])
colnames(meta_col)[colnames(meta_col)%in%c("Trait", "Trait_Class")] = c("trait","trait_class")

m = gt1_calls1[rownames(gt1_calls1)%in%meta_file[,paste0("filename_",type)],col_ind]
rownames(m) = meta_file[match(rownames(m),meta_file[,paste0("filename_",type)]),id_col]

save(meta_col, file=paste0(meta_col_gt_dir,".Rdata"))
if (writecsv) write.csv(meta_col, file=paste0(meta_col_gt_dir,".csv"))

save(m, file=paste0(feat_dna_dir,".Rdata"))
if (writecsv) write.csv(m, file=paste0(feat_dna_dir,".csv"))

## save dominant/recessive versions of dna -----------
m01 = m12 = m
m01[m01==2] = 1
m12[m12==0] = 1
save(m01, file=paste0(feat_dna_dir,"01.Rdata"))
if (writecsv) write.csv(m01, file=paste0(feat_dna_dir,"01.csv"))
save(m12, file=paste0(feat_dna_dir,"12.Rdata"))
if (writecsv) write.csv(m12, file=paste0(feat_dna_dir,"12.csv"))







## save some indices ----------------------------------------------------------------

# get probes/SNP with asthma affiliated gene nearby?
start1 = Sys.time()
asthmarows = apply(meta_col, 1, function(x) any(grepl("asthma",paste(x,collapse=""),ignore.case=T))) & !duplicated(meta_col$dbSNP)
sum(asthmarows)
time_output(start1)

meta_col_asthma_id = meta_col[asthmarows, id_col] #save indices
save(meta_col_asthma_id, file=paste0(meta_col_gt_dir,"_id_asthma-ebiclinvaromim.Rdata"))

asthmarows_gwrns = 
  grepl("asthma",meta_col$gwrngs_Disease.Trait,ignore.case=T) & 
  grepl("European", meta_col$gwrngs_Initial.Sample.Size)
meta_col_asthma_id_gwrns = meta_col[asthmarows_gwrns, id_col] #save indices
save(meta_col_asthma_id_gwrns, file=paste0(meta_col_gt_dir,"_id_asthma-gwrns.Rdata"))

asthmarows_rod0 = read.csv(meta_snp_idrod_dir)
meta_col_asthma_id_rod = meta_col[meta_col$dbSNP%in%asthmarows_rod0[,"SNP"], id_col]
save(meta_col_asthma_id_rod, file=paste0(meta_col_gt_dir,"_id_asthma-rod.Rdata"))

gwrngs.emd$gwrngs_Disease.Trait[gwrngs.emd$gwrngs_SNPs%in%asthmarows_rod0[,"SNP"]]
meta_col$gwrngs_Disease.Trait[meta_col$gwrngs_SNPs%in%asthmarows_rod0[,"SNP"]]

meta_col_asthma_id_st = meta_col[meta_col$dbSNP%in%c("rs993076","rs1800777") | 
                                   apply(meta_col, 1, function(x) 
                                     any(grepl("PLPP3|CETP|farp1|tlr1|bcl10|mal1|card9|card11|pcsk9|serping1",paste(x,collapse=""),ignore.case=T)))
                                 , id_col]
save(meta_col_asthma_id_st, file=paste0(meta_col_gt_dir,"_id_asthma-st.Rdata"))



# get good files
goodpplcols_ind = meta_file[!is.na(meta_file[,"bmi"]),id_col]
save(goodpplcols_ind, file=paste0(meta_file_dir,"_id_goodppl.Rdata"))





## check if all subjects are recorded --------------------

# :( not recorded for dna
# # load data availability, ensure all subjects with rnaseq data are recorded in this script
# meta_file_data = read.csv(meta_file_data_dir)
# length(meta_file_data$UniqueID[meta_file_data$rnaseq=="Y"]) = 
#   sum(meta_file_data$UniqueID[meta_file_data$rnaseq=="Y"] %in% 
#         meta_file$id[!is.na(meta_file$filename_dna)])



