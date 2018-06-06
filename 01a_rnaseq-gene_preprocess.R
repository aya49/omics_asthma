## input: RSEM output .isoforms/genes.results files
## aya43@sfu.ca
## created 20180509
## last modified 20180509

## for convenience, we will analyze the data transposed i.e. col x sample

## based on: https://gist.github.com/jdblischak/11384914

## root directory
root = "~/projects/asthma"
setwd(root)

type = "genes" #"isoforms", "genes"
feature = "count" #"count" or "isopct" for isoforms only
result_dir = paste0(root, "/result/RNAseq/",type); dir.create(result_dir, showWarnings=F)

## input directory
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)
feat_featureraw_dir = paste0(feat_dir,"/",type,"-file-", feature, "raw")

## output directory
feat_feature_dir = paste0(feat_dir,"/",type,"-file-", feature)



## libraries
libr("data.table")
libr("limma")
libr("edgeR")
libr("foreach")
libr("doMC")
libr("stringr")
libr("Matrix")
source("code/_func.R")



## options
no_cores = 15#detectCores()-3
registerDoMC(no_cores)

writecsv = T #write results as csv on top of Rdata?

pthres = .025
padjust = p.adjust.methods

good_col = 3 #each gene must have > good_col samples with >0 abundence; else delete
good_count = 10 #each gene must have >10 abundence in more than half the samples; else delete

id_col = "fileName"
time_col = "time"
class_col = "response"
control = "ER"
categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""










start = Sys.time()

## load data
meta_file = get(load(paste0(meta_file_dir,".Rdata")))
meta_col = get(load(paste0(meta_col_dir,".Rdata")))
m0 = get(load(paste0(feat_featureraw_dir,".Rdata")))
if (!sum(colnames(m0)%in%meta_file[,id_col])>=ncol(m0)) m0 = t(m0)


# keep  probes  that  are expressed  above  background  on  at  least
# n arrays,  where n is  the  smallest  number  of  replicates
# assigned  to  any  of  the  treatment  combinations. 
# 
# filtering methods involving variances should not be used. 
# The limma algorithm analyses the spread of the genewise variances

## filter: get rid of gene with too many 0
no_0 = apply(m0, 2, function(x) sum(x>0)>good_col)
m1 = m0[,no_0]

# alternative
cpm_log <- cpm(m0, log=T)
median_log2_cpm <- apply(cpm_log, 1, median)

x11()
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)

sum(median_log2_cpm > expr_cutoff)
m1 <- m0[median_log2_cpm > expr_cutoff,]

# after removing all genes with a median log2 cpm below r expr_cutoff, 
# we have r sum(median_log2_cpm > expr_cutoff) genes remaining. 
# a good rule of thumb when analyzing RNA-seq data from a single cell type 
# is to expect 9-12 thousand expressed genes.

# recalculate cutoff after filtering
cpm_log <- cpm(m1, log = TRUE)

# 2 group comparison
group = meta_file[,class_col]
y <- DGEList(counts=m1, group=group)

# ttm normalization
y <- calcNormFactors(y)
y$samples

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
y <- estimateDisp(y)
sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y)

# test DE; similar to fisher exact test
et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(m1), sort.by = "none")
head(results_edgeR$table)

# how many genes are differentially expressed at an FDR of 10%?
sum(results_edgeR$table$FDR < .1)
#  MA plot above plots the log2 fold change on the y-axis 
#  versus the average log2 counts-per-million on the x-axis. 
#  The red dots are genes with an FDR less than 10%. 
#  The blue lines represent a four-fold change in expression.
plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .1])
abline(h = c(-2, 2), col = "blue")




# add covariates with glm
y <- DGEList(m1)
y <- calcNormFactors(y)
design <- model.matrix(~group + batch + rin)
design

# coef = 2 corresponds to testing the second column of the design matrix, 
# which in this case is whether the sample is from group
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# example gene
boxplot(as.numeric(m1["HBB", ]) ~ group)









## limma filtering normalization





## filter: get rid of lowly expressed genes in more than half the samples
no_2 = colSums(cpm(m1)>10) >2
m1 = m1[,no_2]






## confounders: age sex bmi race
m2 = voom(m1,design,plot=T)

#cluster libraries
plotMDS(m2, xlim=c(-2.5,-2.5))

#lin model on de
fit = eBayes(lmFit(m2,design))
tt = topTable(fit,coef=2)



