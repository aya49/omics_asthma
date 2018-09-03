## input: rnaseq RSEM output .isoforms/genes.results files
## output: well formatted data and meta col files
## aya43@sfu.ca
## created 20180509
## last modified 20180522


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/src/_dirs.R"))
source(paste0(root, "/src/_func.R"))
source(paste0(root, "/src/_func-asthma.R"))
libr(append(pkgs(), c("TxDb.Hsapiens.UCSC.hg19.knownGene", "annotables", "biomaRt", "org.Hs.eg.db")))


## options
writecsv = T

#plot
wdth=600
height=500

pthres = .025

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

#match cutoffs with amrit's analysis
expr_cutoff = list()
expr_cutoff$genes = 3 #log2 expression sum across samples must be above expr_cutoff
expr_cutoff$isoforms = 3 
expr_cutoff$isopct = 3 


## start
start = Sys.time()


## make meta_file; meta_cell ----------------------------------------

meta_file = get(load(paste0(meta_file_dir,".Rdata")))
meta_fileraw = get(load(paste0(meta_fileall_dir,".Rdata")))
meta_file_pre = meta_fileraw[meta_fileraw[,"time"]=="Pre",]
meta_file_post = meta_fileraw[meta_fileraw[,"time"]=="Post",]


## load & save matrix -----------------------------------
feats = list()
meta_cols = list()

# genes
data_paths = sort(list.files(rsem_dir, pattern=paste0("genes.results"), full.names=T))
data_filenames = fileNames(data_paths, ext=paste0("genes.results"))

# counts0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f5 ",x))))
counts = Reduce("cbind", llply(data_paths, function(x) fread(x, select=5)))

# counts_rownames = read.table(pipe(paste0("cut -f1 ",data_paths[1])))[-1,c(2,1)]
counts_rownames = fread(data_paths[1], select = c(1,2), data.table=F)
colnames(counts_rownames) = c(id_col,"transcript") #g3n3=id

counts = t(counts)

colnames(counts) = counts_rownames[,id_col]
rownames(counts) = data_filenames

save(counts,file=paste0(feat_rnaseq_dir,"genes.raw.Rdata"))
feats$genes = counts
meta_cols$genes = counts_rownames


# isoforms
data_paths = sort(list.files(rsem_dir, pattern=paste0("isoforms.results"), full.names=T))
data_filenames = fileNames(data_paths, ext=paste0("isoforms.results"))

start1 = Sys.time()
# countsis0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f5 ",x))))
countsis = Reduce("cbind",llply(data_paths, function(x) fread(x, select=5)))
# countsis = foreach(x=data_paths, .combine=cbind) %dopar% { return(fread(x, select=5)) }

# isopct0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f8 ",x))))
isopct0 = lapply(data_paths, function(x) fread(x, select=8))
isopct = Reduce("cbind",isopct0)
# isopct = foreach(x=data_paths, .combine=cbind) %dopar% { return(fread(x, select=8)) }
time_output(start1)

# counts_rownames = read.table(pipe(paste0("cut -f1,2 ",data_paths[1])))
counts_rownames = fread(data_paths[1], select = c(1,2), data.table=F)
colnames(counts_rownames) = c(id_col,"gene") #transcript=id
# save(counts_rownames, file=paste0(meta_col_rnaseq_dir,"isoforms.raw.Rdata"))
# save(counts_rownames, file=paste0(meta_col_rnaseq_dir,"isopct.raw.Rdata"))

countsis = t(countsis)
isopct = t(isopct)

rownames(countsis) = rownames(isopct) = data_filenames
colnames(countsis) = colnames(isopct) = counts_rownames[,id_col]

save(countsis,file=paste0(feat_rnaseq_dir,"isoforms.raw.Rdata"))
save(isopct,file=paste0(feat_rnaseq_dir,"isopct.raw.Rdata"))
feats$isoforms = countsis
feats$isopct = isopct
meta_cols$isoforms = meta_cols$isopct = counts_rownames


time_output(start)


# mart = useMart(biomart="ensembl",  
#                # path="biomart/martservice",
#                dataset="hsapiens_gene_ensembl") 

for (feat_type in names(meta_cols)) {
  counts = feats[[feat_type]]
  counts_rownames = meta_cols[[feat_type]]
  by.x = ifelse(feat_type=="genes",id_col,"gene")
  
  geneid = sapply(strsplit(counts_rownames[,by.x], ".", fixed=T), function(x) x[1])
  counts_rownames1 = cbind(counts_rownames, 
                           as.data.frame(grch38)[match(geneid,unlist(grch38[,"ensgene"])),])
  
  # # features = sapply(strsplit(features,"[.]"), function(x) x[1])
  # ensemIDs = 
  #   getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
  #                        'start_position', 'end_position', 'description',
  #                        'entrezgene_trans_name', 'entrezgene'),
  #         filters = 'ensembl_gene_id',
  #         values = geneid,
  #         mart = mart)
  # # meta_col_ensembl.raw = grch38_dt[match(str_extract(colnames(expdata), "ENSG[0-9]+"), grch38_dt$ensgene),
  # #                              c("symbol","chr","start","end","biotype","description")]
  # # meta_col_ensembl = grch38_dt[match(str_extract(colnames(m0), "ENSG[0-9]+"), grch38_dt$ensgene),
  # #                            c("symbol","chr","start","end","biotype","description")]
  # colnames(ensemIDs) = c("id","symbol","chr","start","end","description","entrez_trans","entrez")
  # 
  # grch38_dt <- merge(as.data.frame(grch38_tx2gene), as.data.frame(grch38), by="ensgene")
  # grch38_dts = grch38_dt[match(str_extract(features, "ENSG[0-9]+"),grch38_dt$ensgene),]
  # ensemIDs = ensemIDs[match(str_extract(features, "ENSG[0-9]+"),ensemIDs$id),]
  # ensemIDs[is.na(ensemIDs$id),] = grch38_dts[match(str_extract(features, "ENSG[0-9]+")[is.na(ensemIDs$id)],grch38_dts$ensgene),c("ensgene","symbol","chr","start","end","description")]
  # 
  # ensemIDs$id = counts_rownames
  # ensemIDs$symbol[is.na(ensemIDs$symbol)] = ensemIDs$id[is.na(ensemIDs$symbol)]
  
  save(counts_rownames1, file=paste0(meta_col_rnaseq_dir,feat_type,".raw.Rdata"))
  if (writecsv) write.csv(counts_rownames1, file=paste0(meta_col_rnaseq_dir,feat_type,".raw.csv"))
}


## preprocess --------------------------------------------------------------------
design = model.matrix(~response * centre + sex + time, data=meta_fileraw)

# 1. keep samples where medium cpm > expr_cutoff[[m0n]]
# 3. ttm normalization: mdge = DGEList(counts=m1, group=group); m2 = calcNormFactors(mdge)

feats0 = feats

for (m0n in names(feats)) {
  m0 = t(feats[[m0n]])
  
  # keep  probes  that  are expressed  above  background  on  at  least
  # n arrays,  where n is  the  smallest  number  of  replicates
  # assigned  to  any  of  the  treatment  combinations. 
  # 
  # filtering methods involving variances should not be used. 
  # The limma algorithm analyses the spread of the genewise variances
  
  ## filter: get rid of gene with too many 0
  # no_0 = apply(m0, 1, function(x) sum(x>0)>good_col)
  # m1 = m0[,no_0]
  
  cpm_log = cpm(m0, log=T) #counts per million
  # cpm = apply(m0, 2, function(x) (x/sum(x))*1000000)
  # # the 1 added to log function is to avoid log 0 values
  # log.cpm = log(cpm + 1, 2)
  median_log2_cpm = apply(cpm_log, 1, median)
  png(file=paste0(preprocess_dir,"/", fileNames(feat_rnaseq_dir),m0n, "_rawhist.png"), width=wdth, height=height)
  hist(median_log2_cpm)
  abline(v=expr_cutoff[[m0n]], col = "red")
  graphics.off()
  
  large_count_ind = median_log2_cpm > expr_cutoff[[m0n]]
  m1 = m0[large_count_ind,]
  rownames(m1) = rownames(m0)[large_count_ind]
  meta_cols[[m0n]] = meta_cols[[m0n]][large_count_ind,]
  
  # after removing all genes with a median log2 cpm below r expr_cutoff, 
  # we have r sum(median_log2_cpm > expr_cutoff) genes remaining. 
  # a good rule of thumb when analyzing RNA-seq data from a single cell type 
  # is to expect 9-12 thousand expressed genes.
  
  # recalculate cutoff after filtering
  cpm_log = cpm(m1, log=T, na.omit=T)
  # m1[cpm_log<expr_cutoff[[m0n]]] = NA
  # m1 = na.omit(m1)

  png(file=paste0(preprocess_dir,"/", fileNames(feat_rnaseq_dir),m0n, "_cpmfiltered-heat.png"), width=wdth, height=height)
  heatmap(cor(cpm_log[,order(meta_fileraw[match(colnames(cpm_log),meta_fileraw$filename_rnaseq),class_col])]))
  graphics.off()
  
  # pca <- prcomp(t(cpm_log), scale. = TRUE)
  # x11()
  # plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
  # text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log), col=as.numeric(factor(group)))
  # summary(pca)
  
  # 2 group comparison
  group = meta_fileraw[match(colnames(m1),meta_fileraw[,"filename_rnaseq"]),class_col]
  mdge = DGEList(counts=m1, group=group)
  
  
  # ttm normalization
  m2 = calcNormFactors(mdge)
  # m2$samples
  
  # shrinks variance of the read counts per gene
  # poisson :( assumes the mean and variance are identical, 
  # but it has been found empirically that the variance 
  # in RNA-seq measurements of gene expression are 
  # larger than the mean (termed "overdispersion")
  # so negative binomial distribution
  # 
  # any technical biases are also included in this estimate
  # calculates a dispersion estimate per gene and 
  # shrinks it towards the trended dispersion
  # 
  # shares information across genes to determine a common dispersion
  
  # m3 = estimateDisp(m2,design)
  # sqrt(m4$common.dispersion) # biological coefficient of variation
  # png(file=paste0(preprocess_dir,"/", filenames(feat_rnaseq_dir),ifelse(m0n=="isopct","pct",""), "_cpmfiltered-disper.png"), width=width, height=height)
  # plotBCV(m4)
  # graphics.off()
  
  # confounders: age sex bmi race
  # png(file=paste0(preprocess_dir,"/", fileNames(feat_rnaseq_dir),m0n, "_cpmfiltered-voom.png"), width=wdth, height=height)
  # m3 = voom(m2,design,plot=T)
  # graphics.off()
  
  rownames(m2) = rownames(m1)
  
  # m3log = log2(as.matrix(m3))
  # 
  # png(file=paste0(preprocess_dir,"/", fileNames(feat_rnaseq_dir),ifelse(m0n=="isopct","pct",""), "_hist_2.png"), width=width, height=height)
  # hist(m3log)
  # abline(v=expr_cutoff_2, col = "red")
  # graphics.off()
  
  m4 = as.matrix(cpm(m2, log=T, na.omit=T))
  m4[m4<expr_cutoff[[m0n]]] = NA
  m4 = m4[apply(m4,1,function(x) any(!is.na(x))),
          apply(m4,2,function(x) sum(!is.na(x))>(good_na*nrow(m4)))]
  feats[[m0n]] = t(m4)
  
  # save(m4, file=paste0(feat_rnaseq_dir,".Rdata"))
  # if (writecsv) write.csv(m4, file=paste0(feat_rnaseq_dir,".csv"))
  
  
  ### onwards is tests... skip this part
  try ({
    # test DE; similar to fisher exact test
    m4 = DGEList(counts=m4, group=group)
    et = exactTest(m4)
    results_edgeR = topTags(et, n=nrow(m4), sort.by="none")
    head(results_edgeR$table)
    
    # how many genes are differentially expressed at an FDR of 10%?
    sum(results_edgeR$table$FDR < .1)
    #  MA plot above plots the log2 fold change on the y-axis 
    #  versus the average log2 counts-per-million on the x-axis. 
    #  The red dots are genes with an FDR less than 10%. 
    #  The blue lines represent a four-fold change in expression.
    png(file=paste0(preprocess_dir,"/", fileNames(feat_rnaseq_dir), "_cpmfiltered-smear.png"), width=wdth, height=height)
    plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .1])
    abline(h = c(-2, 2), col = "blue")
    graphics.off()
    
    # add covariates with glm
    y = DGEList(m1)
    y = calcNormFactors(y)
    
    # coef = 2 corresponds to testing the second column of the design matrix, 
    # which in this case is whether the sample is from group
    y = estimateDisp(y, design)
    fit = glmFit(y, design)
    lrt = glmLRT(fit, coef=3)
    lrtt = topTags(lrt)
    
    # example gene
    boxplot(as.numeric(m1[10,]) ~ group)
    
    ## limma filtering normalization
    
    ## filter: get rid of lowly expressed genes in more than half the samples
    no_2 = colSums(cpm(m1)>10) >2 #???
    m1 = m1[,no_2]
    
    #cluster libraries
    x11()
    plotMDS(m2, xlim=c(-2.5,-2.5))
    
    #lin model on de
    fit = eBayes(lmFit(m2,design))
    tt = topTable(fit,coef=2)
    
  })
}


## save -----------------------------
for (m0n in names(feats)) {
  m0 = feats[[m0n]]
  counts_pre = m0[rownames(m0)%in%meta_file$filename_rnaseq.pre,]
  counts_post = m0[rownames(m0)%in%meta_file$filename_rnaseq.post,]
  rownames(counts_pre) = meta_file[match(rownames(counts_pre), meta_file$filename_rnaseq.pre),id_col]
  rownames(counts_post) = meta_file[match(rownames(counts_post), meta_file$filename_rnaseq.post),id_col]
  
  counts_diff = counts_post - counts_pre
  
  save(counts_pre, file=paste0(feat_rnaseq_dir,m0n,".pre.Rdata"))
  if (writecsv) write.csv(counts_pre, file=paste0(feat_rnaseq_dir,m0n,".pre.csv"))
  save(counts_post, file=paste0(feat_rnaseq_dir,m0n,".post.Rdata"))
  if (writecsv) write.csv(counts_post, file=paste0(feat_rnaseq_dir,m0n,".post.csv"))
  
  save(counts_diff, file=paste0(feat_rnaseq_dir,m0n,".diff.Rdata"))
  if (writecsv) write.csv(counts_diff, file=paste0(feat_rnaseq_dir,m0n,".diff.csv"))

  mcol = meta_cols[[m0n]]
  
  features = mcol[,apply(mcol,2,function(x) any(grepl("ENSG",x)))]
  
  # get more data on features
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
  meta_col = cbind(mcol,meta_col)

  save(meta_col, file=paste0(meta_col_rnaseq_dir,m0n,".Rdata"))
  if (writecsv) write.csv(meta_col, file=paste0(meta_col_rnaseq_dir,m0n,".csv"))
  
}


## check if all subjects are recorded --------------------

# load data availability, ensure all subjects with rnaseq data are recorded in this script
meta_file_data = read.csv(meta_file_data_dir)
length(meta_file_data$UniqueID[meta_file_data$rnaseq=="Y"]) == 
  sum(meta_file_data$UniqueID[meta_file_data$rnaseq=="Y"] %in% 
        meta_file$id[!is.na(meta_file$filename_rnaseq.pre)])
