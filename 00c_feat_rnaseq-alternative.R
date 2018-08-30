## input: rnaseq data from rnaseq_analysis-meta
## output: well formatted data and meta col files
## aya43@sfu.ca
## created 20180720


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
libr(append(pkgs(), c("TxDb.Hsapiens.UCSC.hg19.knownGene", "annotables", "biomaRt", "org.Hs.eg.db"))) #org.* annotation packages; can forge own and interact with using library("AnnotationDbi"))))


## options
writecsv = T
expr_cutoff = 3


## load RNA-Seq datasets
load(rnaseqa_data_dir)
# load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_normalized.RDATA")

datalist = list(rnaseqstarensgenes=starEnsemblExp, rnaseqtrinisoforms=trinityGeneIsoCounts, rnasequcscgenes=ucscGeneCounts, rnasequcscisoforms=ucscGeneIsoCounts)


## import clinical datasets
demo = readRDS(meta_file_rnaseqa_dir)
rownames(demo) = paste(demo$concealedID, demo$Time, sep=".")
demo = demo[match(colnames(starEnsemblExp),rownames(demo)),]
demo$NAME = as.character(demo$NAME)
demo$NAME[grepl("WRF",demo$NAME)] = "WRF"

expdata_group = factor(paste0(demo$NAME,demo$Time,demo$CorrectResponse))

preind = as.character(demo$Time)=="Pre"
postind = as.character(demo$Time)=="Post"

# ensembl datasets
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
mart = useDataset("hsapiens_gene_ensembl",mart=ensembl)
reqcol = c('ensembl_gene_id','entrezgene','refseq_mrna','ucsc', 'hgnc_symbol',
           'chromosome_name','start_position','end_position', 'description')
reqcoln = c("ensgene","entrez","refseq_mrna","ucsc","symbol","chr","start","end","description")

## annotatable data set
grch38_dt = merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")
save(grch38_dt, file=paste0(grch38_dir,".Rdata"))
if (writecsv) write.csv(grch38_dt, file=paste0(grch38_dir,".csv"))

for (expdata_n in names(datalist)) {
  cat("\n",expdata_n); Sys.time()
  
  expdata = datalist[[expdata_n]]
  dim(expdata); dim(demo)
  
  if (expdata_n == "rnaseqstarensgenes")
    # remove ERCC controls
    expdata = expdata[!grepl("ERCC-",rownames(expdata)),]
  
  
  # filter data matrix
  expdata = t(as.matrix(calcNormFactors(DGEList(counts=expdata, group=expdata_group))))
  
  m0 = expdata[,apply(expdata,2,function(x) median(x)>expr_cutoff)]
  m0[m0<expr_cutoff] = 0
  
  # prep and save data matrix
  raw_pre = expdata[preind,]
  raw_post = expdata[postind,]
  
  m_pre = m0[preind,]
  m_post = m0[postind,]
  
  rownames(raw_pre) = rownames(m_pre) = as.character(demo$NAME[preind])
  rownames(raw_post) = rownames(m_post) = as.character(demo$NAME[postind])
  
  save(raw_pre, file=paste0(feat_dir,"/",expdata_n,".pre.raw.Rdata"))
  save(raw_post, file=paste0(feat_dir,"/",expdata_n,".post.raw.Rdata"))
  raw_diff = as.matrix(raw_post) - as.matrix(raw_pre)
  save(raw_diff, file=paste0(feat_dir,"/",expdata_n,".diff.raw.Rdata"))
  save(m_pre, file=paste0(feat_dir,"/",expdata_n,".pre.Rdata"))
  write.csv(m_pre, file=paste0(feat_dir,"/",expdata_n,".pre.csv"))
  save(m_post, file=paste0(feat_dir,"/",expdata_n,".post.Rdata"))
  write.csv(m_post, file=paste0(feat_dir,"/",expdata_n,".post.csv"))
  m_diff = m_post-m_pre
  save(m_diff, file=paste0(feat_dir,"/",expdata_n,".diff.Rdata"))
  write.csv(m_diff, file=paste0(feat_dir,"/",expdata_n,".diff.csv"))
}
Sys.time()

## meta_col -------------------------------

for (expdata_n in names(datalist)) {
  m_pre = get(load(paste0(feat_dir,"/",expdata_n,".pre.Rdata")))
  
  cat("\n",expdata_n); Sys.time()
  
  ## make meta_col
  features = colnames(m_pre)
  gene = rep(NA, length(features))
  
  if (grepl("ucsc",expdata_n)) {
    # features = sapply(strsplit(features,"[.]"), function(x) x[1])
    features_ = features
    colidn = 'ucsc'
    colid = grep(colidn,reqcol)
    meta_col = meta_col_ = Reduce('merge', lapply(reqcol[-colid], function(reqcol_) 
      getBM(attributes=c(colidn, reqcol_),
            filters = colidn,
            values = features_,
            mart = mart)) )
    colnames(meta_col) = append(reqcoln[colid],reqcoln[-colid])
    meta_col = meta_col[match(features_,meta_col[,reqcoln[colid]]),]
    grch38_dts = grch38_dt[match(meta_col[,"ensgene"],grch38_dt["ensgene"]),]
    
    intercol = intersect(colnames(meta_col),colnames(grch38_dts))
    meta_col[,intercol] = sapply(intercol, function(x) {
      colx = meta_col[,x]
      colxna = is.na(colx)
      colx[colxna] = grch38_dts[colxna,x]
      colx
    })
    colnames(meta_col)[1] = "id"
    
  } else if (expdata_n == "rnaseqstarensgenes") {
    # features = sapply(strsplit(features,"[.]"), function(x) x[1])
    features_ = str_extract(features, "ENSG[0-9]+")
    colidn = 'ensembl_gene_id'
    colid = grep(colidn,reqcol)
    meta_col = meta_col_ = Reduce('merge', lapply(reqcol[-colid], function(reqcol_) 
      getBM(attributes=c(colidn, reqcol_),
            filters = colidn,
            values = features_,
            mart = mart)) )
    colnames(meta_col) = append(reqcoln[colid],reqcoln[-colid])
    meta_col = meta_col[match(features_,meta_col[,reqcoln[colid]]),]
    grch38_dts = grch38_dt[match(meta_col[,"ensgene"],grch38_dt["ensgene"]),]
    
    intercol = intersect(colnames(meta_col),colnames(grch38_dts))
    meta_col[,intercol] = sapply(intercol, function(x) {
      colx = meta_col[,x]
      colxna = is.na(colx)
      colx[colxna] = grch38_dts[colxna,x]
      colx
    })
    colnames(meta_col)[1] = "id"
    
  } else if (expdata_n == "rnaseqtrinisoforms") {
    trinitymapfile <- read.delim(trinity_map_dir)
    trinitymapfile$Contig <- unlist(lapply(strsplit(
      as.character(trinitymapfile$query_id), "_"), 
      function(i) paste(i[1], i[2], sep="_")))
    trinitymapfile$UniProt <- unlist(lapply(strsplit(unlist(lapply(strsplit(
      as.character(trinitymapfile$subject_id), "\\|"), 
      function(i) i[[2]])), split="_"), function(x) x[1]))
    trinitymapfile$GenSym <- unlist(lapply(strsplit(unlist(lapply(strsplit(
      as.character(trinitymapfile$subject_id), "\\|"), 
      function(i) i[[3]])), split="_"), function(x) x[1]))
    trinitymapfile_ = trinitymapfile
    trinitymapfile = trinitymapfile[match(features,trinitymapfile$query_id),]
    gene <- trinitymapfile$GenSym
    names(gene) <- features 
    trinitymapfile$query_id = features
    colnames(trinitymapfile)[c(1)] = "id"
    
    colidn = 'hgnc_symbol'
    colid = grep(colidn,reqcol)
    meta_col = meta_col_ = Reduce('merge', lapply(reqcol[-colid], function(reqcol_) 
      getBM(attributes=c(colidn, reqcol_),
            filters = colidn,
            values = unique(gene[!is.na(gene)]),
            mart = mart)) )
    colnames(meta_col) = append(reqcol[colid],reqcol[-colid])
    
    meta_col = cbind(trinitymapfile, meta_col[match(trinitymapfile$GenSym,meta_col[,reqcol[colid]]),])
    meta_col$GenSym[is.na(meta_col$GenSym)] = meta_col$hgnc_symbol[is.na(meta_col$GenSym)]
    meta_col = meta_col[,!colnames(meta_col)%in%"hgnc_symbol"]
    colnames(meta_col)[match(c("GenSym","start_position","end_position","chromosome_name"),colnames(meta_col))] = c("symbol","start","end","chr")
  }
  
  # save col_meta
  meta_col$id = features
  meta_col$symbol[is.na(meta_col$symbol)] = meta_col$id[is.na(meta_col$symbol)]
  save(meta_col, file=paste0(meta_col_dir,expdata_n,".Rdata"))
  if (writecsv) write.csv(meta_col, file=paste0(meta_col_dir,expdata_n,".csv"))
}

Sys.time()
