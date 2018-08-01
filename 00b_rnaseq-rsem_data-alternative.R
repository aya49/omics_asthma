## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result")


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = passte0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
eqtl_dir = paste0(stat_dir,"/eqtl"); dir.create(eqtl_dir, showWarnings=F)
gwas_dir = paste0(stat_dir,"/gwas"); dir.create(gwas_dir,showWarnings=F)



source("code/_func.R")
libr(c("amritr","limma","edgeR",
       "annotables","biomaRt",
       "Matrix"))










grch38_dt = merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")




## import clinical datasets
demo = readRDS("~/projects/asthma/data/RNAelements/data/demo/allsitesDemo.rds")
rownames(demo) = paste(demo$concealedID, demo$Time, sep=".")

## load RNA-Seq datasets
load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_rawData.RDATA")
# load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_normalized.RDATA")

datalist = list(rnaseqstarensgenes=starEnsemblExp, rnaseqtrinisoforms=trinityGeneIsoCounts, rnasequcscgenes=ucscGeneCounts, rnasequcscisoforms=ucscGeneIsoCounts)
for (expdata_n in names(datalist)) {
  expdata = datalist[[expdata_n]]
  
  dim(expdata); dim(demoRNASeq)
  # dim(genDats$expdata); dim(demo)
  
  ## use ensembl datasets
  # ttm normalization
  expdata = 
    as.matrix(calcNormFactors(DGEList(counts=expdata[-(60156:nrow(expdata)),], 
                                      group=factor(paste0(demoRNASeq$NAME,demoRNASeq$Time,demoRNASeq$CorrectResponse)))))
  
  # remove ERCC controls
  # ensemblDat = normalizelibSum(expdata[-(60156:nrow(expdata)),])
  ensemblDat = t(cpm(expdata[-(60156:nrow(expdata)),], log=T, na.omit=T))
  ensemblDat = ensemblDat[!duplicated(paste(demoRNASeq$NAME,demoRNASeq$Time)),]
  demoRNASeq = demoRNASeq[!duplicated(paste(demoRNASeq$NAME,demoRNASeq$Time)),]
  m0 = ensemblDat[,apply(ensemblDat,2,function(x) median(x)>3)]
  m0[m0<3] = NA
  meta_file0 = demoRNASeq[rownames(ensemblDat), ]
  
  raw_pre = ensemblDat[demoRNASeq$Time=="Pre",]
  raw_post = ensemblDat[demoRNASeq$Time=="Post",]
  
  m_pre = m0[demoRNASeq$Time=="Pre",]
  m_post = m0[demoRNASeq$Time=="Post",]
  
  rownames(raw_pre) = rownames(m_pre) = as.character(demoRNASeq$NAME[demoRNASeq$Time=="Pre"])
  rownames(raw_post) = rownames(m_post) = as.character(demoRNASeq$NAME[demoRNASeq$Time=="Pre"])
  
  save(raw_pre, file=paste0(feat_dir,"/",expdata_n,".pre.raw.Rdata"))
  save(raw_post, file=paste0(feat_dir,"/",expdata_n,".post.raw.Rdata"))
  raw_diff = raw_post-raw_pre
  save(raw_diff, file=paste0(feat_dir,"/",expdata_n,".diff.raw.Rdata"))
  save(m_pre, file=paste0(feat_dir,"/",expdata_n,".pre.Rdata"))
  save(m_post, file=paste0(feat_dir,"/",expdata_n,".post.Rdata"))
  m_diff = m_post-m_pre
  save(m_diff, file=paste0(feat_dir,"/",expdata_n,".diff.Rdata"))
  
  
  
  ensemIDs = 
    getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                         'start_position', 'end_position', 'description'),
          filters = 'ensembl_gene_id',
          values = str_extract(colnames(raw_pre), "ENSG[0-9]+"),
          mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                         host="grch37.ensembl.org", 
                         path="/biomart/martservice", 
                         dataset="hsapiens_gene_ensembl"))
  # meta_col_ensembl.raw = grch38_dt[match(str_extract(colnames(ensemblDat), "ENSG[0-9]+"), grch38_dt$ensgene),
  #                              c("symbol","chr","start","end","biotype","description")]
  # meta_col_ensembl = grch38_dt[match(str_extract(colnames(m0), "ENSG[0-9]+"), grch38_dt$ensgene),
  #                            c("symbol","chr","start","end","biotype","description")]
  colnames(ensemIDs) = c("id","symbol","chr","start","end","description")
  
  grch38_dt <- merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")
  grch38_dts = grch38_dt[match(str_extract(colnames(raw_pre), "ENSG[0-9]+"),grch38_dt$ensgene),]
  ensemIDs = ensemIDs[match(str_extract(colnames(raw_pre), "ENSG[0-9]+"),ensemIDs$id),]
  ensemIDs[is.na(ensemIDs$id),] = grch38_dts[match(str_extract(colnames(raw_pre), "ENSG[0-9]+")[is.na(ensemIDs$id)],grch38_dts$ensgene),c("ensgene","symbol","chr","start","end","description")]
  ensemIDs$id = colnames(raw_pre)
  ensemIDs$symbol[is.na(ensemIDs$symbol)] = ensemIDs$id[is.na(ensemIDs$symbol)]
  save(ensemIDs, file=paste0(meta_col_dir,"-",expdata_n,".Rdata"))
  
}