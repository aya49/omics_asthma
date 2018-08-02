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

trinity_map_dir = paste0(root, "/data/RNAseq/asthma.trinity.blastx.outfmt6.txt")


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
demoRNASeq$NAME = as.character(demoRNASeq$NAME)
demoRNASeq$NAME[grepl("WRF",demoRNASeq$NAME)] = "WRF"

expdata_group = factor(paste0(demoRNASeq$NAME,demoRNASeq$Time,demoRNASeq$CorrectResponse))

meta_file0 = demoRNASeq
preind = as.character(meta_file0$Time)=="Pre"
postind = as.character(meta_file0$Time)=="Post"

## use ensembl datasets
# ttm normalization
# meta_col
mart = useMart(biomart="ensembl",  
               # path="biomart/martservice",
               dataset="hsapiens_gene_ensembl") 
# host = "www.ensembl.org", 
# ensemblRedirect = FALSE)


## load RNA-Seq datasets
load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_rawData.RDATA")
# load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_normalized.RDATA")

datalist = list(rnaseqstarensgenes=starEnsemblExp, rnaseqtrinisoforms=trinityGeneIsoCounts, rnasequcscgenes=ucscGeneCounts, rnasequcscisoforms=ucscGeneIsoCounts)

for (expdata_n in names(datalist)) {
  expdata = datalist[[expdata_n]]
  
  dim(expdata); dim(demoRNASeq)
  # dim(genDats$expdata); dim(demo)
  
  features=rownames(expdata)
  gene <- rep(NA, length(features))
  
  if (expdata_n == "rnaseqstarensgenes") {
    # remove ERCC controls
    expdata_temp = expdata[!grepl("ERCC-",rownames(expdata)),]
    
    # features = sapply(strsplit(features,"[.]"), function(x) x[1])
    ensemIDs = 
      getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                           'start_position', 'end_position', 'description'),
            filters = 'ensembl_gene_id',
            values = str_extract(features, "ENSG[0-9]+"),
            mart = mart)
    # meta_col_ensembl.raw = grch38_dt[match(str_extract(colnames(expdata), "ENSG[0-9]+"), grch38_dt$ensgene),
    #                              c("symbol","chr","start","end","biotype","description")]
    # meta_col_ensembl = grch38_dt[match(str_extract(colnames(m0), "ENSG[0-9]+"), grch38_dt$ensgene),
    #                            c("symbol","chr","start","end","biotype","description")]
    colnames(ensemIDs) = c("id","symbol","chr","start","end","description")
    
    grch38_dt <- merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")
    grch38_dts = grch38_dt[match(str_extract(features, "ENSG[0-9]+"),grch38_dt$ensgene),]
    ensemIDs = ensemIDs[match(str_extract(features, "ENSG[0-9]+"),ensemIDs$id),]
    ensemIDs[is.na(ensemIDs$id),] = grch38_dts[match(str_extract(features, "ENSG[0-9]+")[is.na(ensemIDs$id)],grch38_dts$ensgene),c("ensgene","symbol","chr","start","end","description")]
    
  } else if (grepl("ucsc",expdata_n)) {
    ensemIDs <- getBM(attributes = c('refseq_mrna','ucsc','hgnc_symbol', 'chromosome_name',
                                     'start_position', 'end_position', 'description'),
                      filters = "ucsc", 
                      values = features, mart = mart)
    ensemIDs = ensemIDs[match(rownames(expdata),ensemIDs$ucsc),]
    ensemIDs$id = rownames(expdata)
    colnames(ensemIDs) = c("refseq_mrna","ucsc","symbol","chr","start","end")
  } else if (expdata_n == "rnaseqtrinisoforms") {
    expdata_temp = expdata
    
    trinityMapFile <- read.delim(trinity_map_dir)
    trinityMapFile$Contig <- 
      unlist(lapply(strsplit(as.character(trinityMapFile$query_id), 
                             "_"), function(i) paste(i[1], i[2], sep = "_")))
    trinityMapFile$UniProt <- 
      unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(trinityMapFile$subject_id), 
                                                    "\\|"), function(i) i[[2]])), split = "_"), function(x) x[1]))
    trinityMapFile$GenSym <- 
      unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(trinityMapFile$subject_id), 
                                                    "\\|"), function(i) i[[3]])), split = "_"), function(x) x[1]))
    gene <- trinityMapFile$GenSym[match(features,trinityMapFile$query_id)]
    
    ensemIDs_ = 
      getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                           'start_position', 'end_position', 'description'),
            filters = 'hgnc_symbol', values = gene, mart = mart)
    ensemIDs = data.frame(id=features, symbol=gene)
    ensemIDs = cbind(ensemIDs, ensemIDs_[match(ensemIDs$symbol,ensemIDs_$hgnc_symbol),c(3:6)])
    colnames(ensemIDs) = c("id","symbol","chr","start","end","description")
  }
  
  # save col_meta
  ensemIDs$id = rownames(expdata)
  ensemIDs$symbol[is.na(ensemIDs$symbol)] = ensemIDs$id[is.na(ensemIDs$symbol)]
  save(ensemIDs, file=paste0(meta_col_dir,"-",expdata_n,".Rdata"))
  
  # filter data matrix
  expdata = as.matrix(calcNormFactors(DGEList(counts=expdata_temp, group=expdata_group)))
  expdata = t(expdata)
  
  m0 = expdata[,apply(expdata,2,function(x) median(x)>3)]
  m0[m0<3] = NA
  
  # prep and save data matrix
  raw_pre = expdata[preind,]
  raw_post = expdata[postind,]
  
  m_pre = m0[preind,]
  m_post = m0[postind,]
  
  rownames(raw_pre) = rownames(m_pre) = as.character(meta_file0$NAME[preind])
  rownames(raw_post) = rownames(m_post) = as.character(meta_file0$NAME[postind])
  
  save(raw_pre, file=paste0(feat_dir,"/",expdata_n,".pre.raw.Rdata"))
  save(raw_post, file=paste0(feat_dir,"/",expdata_n,".post.raw.Rdata"))
  raw_diff = raw_post-raw_pre
  save(raw_diff, file=paste0(feat_dir,"/",expdata_n,".diff.raw.Rdata"))
  save(m_pre, file=paste0(feat_dir,"/",expdata_n,".pre.Rdata"))
  save(m_post, file=paste0(feat_dir,"/",expdata_n,".post.Rdata"))
  m_diff = m_post-m_pre
  save(m_diff, file=paste0(feat_dir,"/",expdata_n,".diff.Rdata"))
}