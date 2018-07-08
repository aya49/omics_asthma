## input: RSEM output .isoforms/genes.results files; must do genotype data first if we want to merge demographics meta_file
## output: raw cell counts matrix
## aya43@sfu.ca
## created 20180509
## last modified 20180522

## root directory
root = "~/projects/asthma"
setwd(root)

# dir.create(paste0(root, "/result"), showWarnings=F)

# type_ = "isoforms" #"isoforms", "genes"
# type = paste0("rnaseq",type_)
result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)

## input directory
data_dir = paste0(root,"/data/RNAseq")
rsem_dir = paste0(data_dir,"/rsem")
meta_file_temp1_dir = paste0(data_dir,"/RNASeq.asthma.clinical_sequencing_merged.csv")
meta_file_temp2_dir = paste0(data_dir,"/rnaseq_demo.Rdata")
meta_file_data_dir = paste0(root,"/data/asthmaDemo_allsite_withSampleInfo_DH_v5.csv")
meta_col_temp_dir = paste0(data_dir,"/HuGene-2_1-st-v1.na36.hg19.transcript.csv")
meta_col_tr_temp_dir = paste0(data_dir,"/HuGene-2_1-st-v1.na36.hg19.probeset.csv")

## output directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col-rnaseq")

feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_feature_dir = paste0(feat_dir,"/rnaseq")
feat_cell_dir = paste0(feat_dir,"/cell")

stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
preprocess_dir = paste0(stat_dir,"/stat"); dir.create(preprocess_dir, showWarnings=F)

## libraries
source("code/_func.R")
libr("data.table")
libr("annotables") #grch38; https://github.com/stephenturner/annotables
libr("limma")
libr("edgeR")
libr("stringr")
libr("gdata") #read xls
libr("Matrix")


writecsv = T
#what to name the id column
cid_col = "id"
id_col = "id"
class_col = "response"

#plot
wdth=600
ht=500

pthres = .025

#match with amrit's analysis
expr_cutoff = list()
expr_cutoff$genes = 3 #log2 expression sum across samples must be above expr_cutoff
expr_cutoff$isoforms = 3 
expr_cutoff$isopct = 3 


good_col = 3 #each gene must have > good_col samples with >0 abundence; else delete
good_count = 10 #each gene must have >10 abundence in more than half the samples; else delete




start = Sys.time()


## load & save matrix -----------------------------------

feats = list()
meta_cols = list()

# genes

data_paths = sort(list.files(rsem_dir, pattern=paste0("genes.results"), full.names=T))
data_filenames = fileNames(data_paths, ext=paste0("genes.results"))

# counts0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f5 ",x))))
counts0 = lapply(data_paths, function(x) fread(x, select=5))
counts = Reduce("cbind",counts0)

# counts_rownames = read.table(pipe(paste0("cut -f1 ",data_paths[1])))[-1,c(2,1)]
counts_rownames = fread(data_paths[1], select = c(1,2), data.table=F)
colnames(counts_rownames) = c(id_col,"transcript") #g3n3=id

counts = t(counts)

colnames(counts) = counts_rownames[,id_col]
rownames(counts) = data_filenames

save(counts,file=paste0(feat_feature_dir,"genes.raw.Rdata"))
feats$genes = counts
meta_cols$genes = counts_rownames


# isoforms

data_paths = sort(list.files(rsem_dir, pattern=paste0("isoforms.results"), full.names=T))
data_filenames = fileNames(data_paths, ext=paste0("isoforms.results"))

start1 = Sys.time()
# countsis0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f5 ",x))))
countsis0 = lapply(data_paths, function(x) fread(x, select=5))
countsis = Reduce("cbind",countsis0)
# countsis = foreach(x=data_paths, .combine=cbind) %dopar% { return(fread(x, select=5)) }

# isopct0 = lapply(data_paths, function(x) read.table(pipe(paste0("cut -f8 ",x))))
isopct0 = lapply(data_paths, function(x) fread(x, select=8))
isopct = Reduce("cbind",isopct0)
# isopct = foreach(x=data_paths, .combine=cbind) %dopar% { return(fread(x, select=8)) }
time_output(start1)

# counts_rownames = read.table(pipe(paste0("cut -f1,2 ",data_paths[1])))
counts_rownames = fread(data_paths[1], select = c(1,2), data.table=F)
colnames(counts_rownames) = c(id_col,"gene") #transcript=id
# save(counts_rownames, file=paste0(meta_col_dir,"isoforms.raw.Rdata"))
# save(counts_rownames, file=paste0(meta_col_dir,"isopct.raw.Rdata"))

countsis = t(countsis)
isopct = t(isopct)

rownames(countsis) = rownames(isopct) = data_filenames
colnames(countsis) = colnames(isopct) = counts_rownames[,id_col]

save(countsis,file=paste0(feat_feature_dir,"isoforms.raw.Rdata"))
save(isopct,file=paste0(feat_feature_dir,"isopct.raw.Rdata"))
feats$isoforms = countsis
feats$isopct = isopct
meta_cols$isoforms = meta_cols$isopct = counts_rownames


time_output(start)



for (feat_type in names(meta_cols)) {
  counts = feats[[feat_type]]
  counts_rownames = meta_cols[[feat_type]]
  by.x = ifelse(feat_type=="genes",id_col,"gene")
  
  geneid = sapply(strsplit(counts_rownames[,by.x], ".", fixed=T), function(x) x[1])
  counts_rownames1 = cbind(counts_rownames, 
                           as.data.frame(grch38)[match(geneid,unlist(grch38[,"ensgene"])),])
  save(counts_rownames1, file=paste0(meta_col_dir,feat_type,".raw.Rdata"))
  if (writecsv) write.csv(counts_rownames1, file=paste0(meta_col_dir,feat_type,"raw.csv"))
}








## make meta_file; meta_cell ----------------------------------------

meta_file1_temp0 = fread(meta_file_temp1_dir, data.table=F)

ucol = col_probe(meta_file1_temp0)
meta_file1_temp1 = meta_file1_temp0[,-ucol$u1]
meta_file1_temp1[,"Filename"] = gsub(".bam","",meta_file1_temp1[,"Filename"])
ucol = col_probe(meta_file1_temp1)

meta_file1_temp2 = get(load(meta_file_temp2_dir))
meta_file1_temp2 = data.frame(lapply(meta_file1_temp2, as.character), stringsAsFactors=FALSE)

# meta_file2_temp0 = fread(meta_file_temp2_dir, data.table=F)
# ucol = col_probe(meta_file2_temp0)
# meta_file2_temp1 = meta_file2_temp0[,-ucol$u1]
# ucol = col_probe(meta_file2_temp1)


# meta_file1_temp = meta_file1_temp1[,c("Filename", "Subject", "Phenotype", "Time", "Allergen", "SITE", 
#                                       "RACE", "SEX", "HT.cm.", "AGE", "BLFEV", 
#                                       "Volume", "Concentration", "Quantity", "Extracted")]
# colnames(meta_file1_temp) = c("filename", "sample", "response", "time", "allergen", "centre", 
#                               "race", "sex", "height", "weight", "age", "blfev", 
#                               "vol", "conc", "qty", "extracted")

meta_file1_temp = meta_file1_temp2[,c("Read.Set.Id", "UniqueID", "CorrectResponse", "Time", "Allergen_cleanLabel", "SITE",
                                      "RACE", "SEX", "HT.cm.", "Wt..Kg.", "AGE", "PRFEV", "BLFEV", "Cohort")]
colnames(meta_file1_temp) = c("filename", id_col, "response", "time", "allergen", "centre", 
                              "race", "sex", "height", "weight", "age", "prfev", "blfev","cohort")
meta_file1_temp[grepl("WRF",meta_file1_temp[,id_col]),id_col] = "WRF"
# meta_file1_temp[,id_col] = gsub(".bam","",meta_file1_temp[,"filename"])

cell = meta_file1_temp1[match(meta_file1_temp[,"filename"],meta_file1_temp1[,"Filename"]),c(31:ncol(meta_file1_temp1))]
rownames(cell) = meta_file1_temp[,"filename"]

mfc_order = match(data_filenames,meta_file1_temp[,"filename"])
meta_fileraw = meta_file1_temp[mfc_order,]


cell = cell[mfc_order,]
cell = delna(cell)
save(cell, file=paste0(feat_cell_dir,".raw.Rdata"))

cell_pre = cell[meta_fileraw[match(rownames(cell), meta_fileraw[,"filename"]),"time"]=="Pre",]
rownames(cell_pre) = meta_fileraw[match(rownames(cell_pre), meta_fileraw[,"filename"]),id_col]
save(cell_pre, file=paste0(feat_cell_dir,".pre.Rdata"))
if (writecsv) write.csv(cell_pre, file=paste0(feat_cell_dir,".pre.csv"))

cell_post = cell[meta_fileraw[match(rownames(cell), meta_fileraw[,"filename"]),"time"]=="Post",]
rownames(cell_post) = meta_fileraw[match(rownames(cell_post), meta_fileraw[,"filename"]),id_col]
save(cell_post, file=paste0(feat_cell_dir,".post.Rdata"))
if (writecsv) write.csv(cell_post, file=paste0(feat_cell_dir,".post.csv"))









## preprocess ------------------------------------------------------------------------------------------------------

design = model.matrix(~response * centre + sex + time, data=meta_fileraw)

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
  png(file=paste0(preprocess_dir,"/", fileNames(feat_feature_dir),m0n, "_rawhist.png"), width=wdth, height=ht)
  hist(median_log2_cpm)
  abline(v=expr_cutoff[[m0n]], col = "red")
  graphics.off()
  
  large_count_ind = median_log2_cpm > expr_cutoff[[m0n]]
  m1 = m0[large_count_ind,]
  m1[m1<3] = NA
  meta_cols[[m0n]] = meta_cols[[m0n]][large_count_ind,]
  
  rownames(m1) = rownames(m0)[large_count_ind]
  
  # after removing all genes with a median log2 cpm below r expr_cutoff, 
  # we have r sum(median_log2_cpm > expr_cutoff) genes remaining. 
  # a good rule of thumb when analyzing RNA-seq data from a single cell type 
  # is to expect 9-12 thousand expressed genes.
  
  # recalculate cutoff after filtering
  cpm_log = cpm(m1, log=T)
  # m1[cpm_log<expr_cutoff[[m0n]]] = NA
  # m1 = na.omit(m1)
  
  png(file=paste0(preprocess_dir,"/", fileNames(feat_feature_dir),m0n, "_cpmfiltered-heat.png"), width=wdth, height=ht)
  heatmap(cor(cpm_log[,order(meta_fileraw[class_col,])]))
  graphics.off()
  
  # pca <- prcomp(t(cpm_log), scale. = TRUE)
  # x11()
  # plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
  # text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log), col=as.numeric(factor(group)))
  # summary(pca)
  
  # 2 group comparison
  group = meta_fileraw[match(meta_fileraw[,"filename"],colnames(m1)),class_col]
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
  # png(file=paste0(preprocess_dir,"/", filenames(feat_feature_dir),ifelse(m0n=="isopct","pct",""), "_cpmfiltered-disper.png"), width=width, height=height)
  # plotBCV(m4)
  # graphics.off()
  
  # confounders: age sex bmi race
  # png(file=paste0(preprocess_dir,"/", fileNames(feat_feature_dir),m0n, "_cpmfiltered-voom.png"), width=wdth, height=ht)
  # m3 = voom(m2,design,plot=T)
  # graphics.off()
  
  rownames(m2) = rownames(m1)
  
  
  
  # m3log = log2(as.matrix(m3))
  # 
  # png(file=paste0(preprocess_dir,"/", fileNames(feat_feature_dir),ifelse(m0n=="isopct","pct",""), "_hist_2.png"), width=width, height=height)
  # hist(m3log)
  # abline(v=expr_cutoff_2, col = "red")
  # graphics.off()
  
  m4 = m2
  m4[m4<expr_cutoff[[m0n]]] = NA
  feats[[m0n]] = t(as.matrix(m4))
  
  # save(m4, file=paste0(feat_feature_dir,".Rdata"))
  # if (writecsv) write.csv(m4, file=paste0(feat_feature_dir,".csv"))
  
  
  ### onwards is tests... incomplete
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
    png(file=paste0(preprocess_dir,"/", fileNames(feat_feature_dir), "_cpmfiltered-smear.png"), width=wdth, height=hy)
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
    no_2 = colSums(cpm(m1)>10) >2
    m1 = m1[,no_2]
    
    #cluster libraries
    x11()
    plotMDS(m2, xlim=c(-2.5,-2.5))
    
    #lin model on de
    fit = eBayes(lmFit(m2,design))
    tt = topTable(fit,coef=2)
    
  })
}







## split and save meta_file meta_cell -------------------------------------------
meta_file_pre = meta_fileraw[meta_fileraw[,"time"]=="Pre",]
meta_file_post = meta_fileraw[meta_fileraw[,"time"]=="Post",]

if (file.exists(paste0(meta_file_dir,".Rdata"))) {
  meta_file_old = get(load(paste0(meta_file_dir,".Rdata")))
  meta_file_old$id[meta_file_old$rnaseq_coreSet_biomarkerAnalysis=="Y"]
  meta_file_pre = merge.data.frame(meta_file_pre[,append(id_col,setdiff(colnames(meta_file_pre),colnames(meta_file_old)))],meta_file_old,all=T,by=id_col)
  meta_file_post = merge.data.frame(meta_file_post[,append(id_col,setdiff(colnames(meta_file_post),colnames(meta_file_old)))],meta_file_old,all=T,by=id_col)
}
meta_file_pre = meta_file_pre[order(meta_file_pre[,id_col]),]
meta_file_post = meta_file_post[order(meta_file_post[,id_col]),]
meta_file = meta_file_post
meta_file$filename_rnaseq.pre = meta_file_pre[,"filename"]
meta_file$filename_rnaseq.post = meta_file_post[,"filename"]

meta_file = meta_file[,!colnames(meta_file)%in%c("filename","time")]


save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"))



## save -----------------------------
for (m0n in names(feats)) {
  m0 = feats[[m0n]]
  counts_pre = m0[meta_file$filename_rnaseq.pre[!is.na(meta_file$filename_rnaseq.pre)],]
  counts_post = m0[meta_file$filename_rnaseq.post[!is.na(meta_file$filename_rnaseq.post)],]
  rownames(counts_pre) = rownames(counts_post) = meta_file[!is.na(meta_file$filename_rnaseq.post),id_col]
  
  counts_diff = counts_post - counts_pre
  
  save(counts_pre, file=paste0(feat_feature_dir,m0n,".pre.Rdata"))
  if (writecsv) write.csv(counts_pre, file=paste0(feat_feature_dir,m0n,".pre.csv"))
  save(counts_post, file=paste0(feat_feature_dir,m0n,".post.Rdata"))
  if (writecsv) write.csv(counts_post, file=paste0(feat_feature_dir,m0n,".post.csv"))
  
  save(counts_diff, file=paste0(feat_feature_dir,m0n,".diff.Rdata"))
  if (writecsv) write.csv(counts_diff, file=paste0(feat_feature_dir,m0n,".diff.csv"))

  mcol = meta_cols[[m0n]]
  save(mcol, file=paste0(meta_col_dir,m0n,".Rdata"))
  if (writecsv) write.csv(mcol, file=paste0(meta_col_dir,m0n,".csv"))
  
}




## check if all subjects are recorded --------------------

# load data availability, ensure all subjects with rnaseq data are recorded in this script
meta_file_data = read.csv(meta_file_data_dir)
length(meta_file_data$UniqueID[meta_file_data$rnaseq=="Y"]) = 
  sum(meta_file_data$UniqueID[meta_file_data$rnaseq=="Y"] %in% 
        meta_file$id[!is.na(meta_file$filename_rnaseq.pre)])
