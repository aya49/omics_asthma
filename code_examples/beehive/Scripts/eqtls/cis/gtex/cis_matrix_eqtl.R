##################################################
#  cis_matrix_eqtl.R
#
#  $proj/Scripts/eqtls/cis/gtex/cis_matrix_eqtl.R
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# Example
# args = c(1:5)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/muscleskeletal_nonverlapping_certain_autosomes_normalized.txt'
# args[2] = '1'
# args[3] = '1000000'
# args[4] = 'muscleskeletal'
# args[5] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
# args[6] = '/tigress/BEE/RNAseq/Output/cis-mapping/gtex/MatrixEQTL/gct-normalized-1M/muscleskeletal/muscleskeletal_nonverlapping_certain_autosomes_normalized_MatrixEQTL'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
geno_option = 'continuous'
chr_number = args[2]
cis_dist = as.numeric(args[3])
tissue_name = args[4]
cov_dir = args[5]
out_file = args[6]
# Take all SNPs for the trans runs
# genotype_file_name = paste(args[8], 'Genotype/GTEx/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_', disc_or_cont , '_Chr', chr_number, '_Final.txt', sep="")
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_Chr', chr_number, '_Final.txt', sep="")
SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', chr_number, '.txt', sep="")
gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chr',chr_number,'.txt')

# Need to check that the subjects read in from expression matrices match the subjects in genotype data
library(dplyr)

expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]

# Make sure genotype and expression matrix columns line up
genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

# Also obtain genes and snps that are in the expression and genotype files, respectively.
gene_positions = read.table(file = gene_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gene_positions = gene_positions[,c('gene_id', 'chr', 'start', 'end')]
gene_positions = mutate(gene_positions, chr = paste("chr", chr, sep = ""))
gene_positions = select(gene_positions, c(gene_id, chr, start, end))

snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(snp_positions) = snp_positions$rsID
snp_positions = snp_positions[rownames(snp_positions) %in% rownames(genotype_matrix),]
genotype_matrix = genotype_matrix[rownames(genotype_matrix) %in% rownames(snp_positions),]

# Make sure the list of genes and their ordering are the same
# We are also filtering out non-autosomal genes
expression_matrix = expression_matrix[rownames(expression_matrix) %in% gene_positions$gene_id,]
rownames(gene_positions) = gene_positions$gene_id
gene_positions = gene_positions[rownames(expression_matrix),]
# gene_positions = gene_positions[sapply(rownames(expression_matrix), function(x) {match(x, gene_positions$gene_id)}),]
# rownames(gene_positions) = gene_positions$gene_id

# Incorporate a tissue-specific MAF filter so that we are not testing variants that actually have a low MAF in our test set
gene_set_size = dim(expression_matrix)[1]
snp_set_size_init = dim(genotype_matrix)[1]

MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
genotype_matrix = genotype_matrix[MAF>0.05 & MAF<0.95,]
snp_positions = snp_positions[MAF>0.05 & MAF<0.95,]

snp_set_size_final = dim(genotype_matrix)[1]

# if(method == 'gtex_gct_v6p') {

  suffix = '_Analysis.covariates.txt'
  # covar_mats = list.files(path = cov_dir, pattern = suffix)
  # ind = which(as.character(sapply(covar_mats, function(x) {tolower(gsub(" ", "", gsub("[[:punct:]]", "", strsplit(x, suffix)[[1]][1])))})) == tissue_name)
  covars = read.table(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  rownames(covars) = covars$ID
  covars = covars[,colnames(expression_matrix)]

# } else {

#   covars = read.table(paste(args[7], 'Covariates/GTEx/genotype_platform_gender_covariates.txt', sep=''), sep='\t', stringsAsFactors=FALSE, header=TRUE)
#   rownames(covars) = covars$X
#   covars = covars[colnames(expression_matrix),]
#   covars = covars[,(2:6)]
#   if (length(unique(covars$GENDER)) == 1) {
#     covars = subset(covars, select = -(GENDER) )
#   }
#   if (length(unique(covars$Platform)) == 1) {
#     covars = subset(covars, select = -(Platform) )
#   }
#   covars = as.data.frame(t(covars))

# }

# Import MatrixEQTL
library(MatrixEQTL)

cvrt = SlicedData$new()
cvrt$CreateFromMatrix(as.matrix(covars));

snps = SlicedData$new()
snps$CreateFromMatrix(as.matrix(genotype_matrix));

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1e-2

#for (i in c(1:2)) {
useModel = modelLINEAR;

gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(expression_matrix));

  # gene_qnorm = gene$copy();
  # # Quantile-normalize the expression values - ties broken by averaging
  # for( sl in 1:length(gene_qnorm) ) {
  #   mat = gene_qnorm[[sl]];
  #   mat = t(apply(mat, 1, rank, ties.method = "average"));
  #   mat = qnorm(mat / (ncol(gene_qnorm)+1));
  #   gene_qnorm[[sl]] = mat;
  # }
  # rm(sl, mat);

me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = NULL,
      # record cis only
      pvOutputThreshold = 0,
      useModel = useModel,
      verbose = TRUE,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = pvOutputThreshold,
      snpspos = snp_positions,
      genepos = gene_positions,
      cisDist = cis_dist,
      pvalue.hist = 100,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)

me_cis = me$cis

me_cis$eqtls$tissue_MAF = rowSums(genotype_matrix[as.character(me_cis$eqtls$snps),])/dim(genotype_matrix)[2]/2
me_cis$eqtls$tissue_MAF = sapply(me_cis$eqtls$tissue_MAF, function(x) {min(x, 1-x)})

# We also repeat with the permuted genotypes: however, for all-by-all, this is not for calibration but for simply showing
# that the permuted runs don't deviate much from the uniform assumption.

# A strict comparison is hard because the permutations are not tissue-specific and it can change the MAF for particular variants
geno_option = 'continuous_permute'
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_Chr', chr_number, '_Final.txt')
genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]
# if (multiplier==10) {
#   genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):dim(genotype_matrix)[1]),]
# } else {
#   genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):((multiplier)*snp_set_size)),]
# }

snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(snp_positions) = snp_positions$rsID
snp_positions = snp_positions[rownames(snp_positions) %in% rownames(genotype_matrix),]
genotype_matrix = genotype_matrix[rownames(genotype_matrix) %in% rownames(snp_positions),]
# Incorporate a tissue-specific MAF filter so that we are not testing variants that actually have a low MAF in our test set
MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
genotype_matrix = genotype_matrix[MAF>0.05 & MAF<0.95,]
snp_positions = snp_positions[MAF>0.05 & MAF<0.95,]

snps = SlicedData$new()
snps$CreateFromMatrix(as.matrix(genotype_matrix));

useModel = modelLINEAR;

gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(expression_matrix));

me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = NULL,
      pvOutputThreshold = 0,
      useModel = useModel,
      verbose = TRUE,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = pvOutputThreshold,
      snpspos = snp_positions,
      genepos = gene_positions,
      cisDist = cis_dist,
      pvalue.hist = 100,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)

me_cis_permute = me$cis

me_cis_permute$eqtls$tissue_MAF = rowSums(genotype_matrix[as.character(me_cis_permute$eqtls$snps),])/dim(genotype_matrix)[2]/2
me_cis_permute$eqtls$tissue_MAF = sapply(me_cis_permute$eqtls$tissue_MAF, function(x) {min(x, 1-x)})


# Save the Matrix EQTL object as RData
save(me_cis, me_cis_permute, file=paste0(out_file, '_part', chr_number, '_me.RData'))
