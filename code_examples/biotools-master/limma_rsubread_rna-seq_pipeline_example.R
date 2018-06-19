# Case study: using a Bioconductor R pipeline to analyze RNA-seq data
#
# Wei Shi (shi at wehi dot edu dot au), Yang Liao and Gordon K Smyth
# Bioinformatics Division, Walter and Eliza Hall Institute, Melbourne, Australia
#
# This file is adapted from: http://bioinf.wehi.edu.au/RNAseqCaseStudy/
#     Code: http://bioinf.wehi.edu.au/RNAseqCaseStudy/code.R
#     Data: http://bioinf.wehi.edu.au/RNAseqCaseStudy/data.tar.gz

# The version of Rsubread package should be 1.12.1 or later. 
# You should run R version 3.0.2 or later. 

# Record the time
print(Sys.time())

# load libraries
library(Rsubread)
library(limma)
library(edgeR)

# read in target file
options(digits=2)
targets <- readTargets()
targets

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index for reference sequence (Chr1 in hg19)
buildindex(basename="chr1",reference="hg19_chr1.fa")

# align reads
align(index="chr1", readfile1=targets$InputFile, 
      input_format="gzFASTQ", output_format="BAM", 
      output_file=targets$OutputFile, 
      tieBreakHamming=TRUE, unique=TRUE, indels=5)

# count numbers of reads mapped to NCBI Refseq genes
fc <- featureCounts(files=targets$OutputFile,annot.inbuilt="hg19")
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

# generate RPKM values if you need them
x_rpkm <- rpkm(x,x$genes$Length)

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# perform voom normalization
y <- voom(x,design,plot=TRUE)

# cluster libraries
plotMDS(y,xlim=c(-2.5,2.5))

# fit linear model and assess differential expression
fit <- eBayes(lmFit(y,design))
topTable(fit,coef=2)

# Record the time
print(Sys.time())
