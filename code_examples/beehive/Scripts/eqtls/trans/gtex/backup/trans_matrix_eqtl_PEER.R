##################################################
#  trans_matrix_eqtl.R
#
#  $proj/Scripts/eqtls/trans/gtex/trans_matrix_eqtl.R
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
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/muscleskeletal_Analysis.covariates.txt'
# args[2] = '1'
# args[3] = '2.5e8'
# args[4] = 'muscleskeletal'
# args[5] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/cov_association/muscleskeletal/muscleskeletal_cov_association_MatrixEQTL'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
geno_option = 'continuous'
chr_number = as.numeric(args[2])
cis_dist = as.numeric(args[3])
tissue_name = args[4]
out_file = args[5]

# Take all SNPs for the trans runs
# genotype_file_name = paste(args[8], 'Genotype/GTEx/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_', disc_or_cont , '_Chr', chr_number, '_Final.txt', sep="")
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr_number, '_Final.txt', sep="")
SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', chr_number, '.txt', sep="")

# Need to check that the subjects read in from expression matrices match the subjects in genotype data
library(dplyr)
library(matrixStats)

expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(expression_matrix) = expression_matrix$ID
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]

genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

# Just need a pseudo gene position matrix
gene_positions = data.frame(gene_id = rownames(expression_matrix), chr = 'chrC', start = 1, end = 2)

snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(snp_positions) = snp_positions$rsID
snp_positions = snp_positions[rownames(snp_positions) %in% rownames(genotype_matrix),]
genotype_matrix = genotype_matrix[rownames(genotype_matrix) %in% rownames(snp_positions),]

# Incorporate a tissue-specific MAF filter so that we are not testing variants that actually have a low MAF in our test set
gene_set_size = dim(expression_matrix)[1]
snp_set_size_init = dim(genotype_matrix)[1]

MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
genotype_matrix = genotype_matrix[MAF>0.05 & MAF<0.95,]
snp_positions = snp_positions[MAF>0.05 & MAF<0.95,]

snp_set_size_final = dim(genotype_matrix)[1]

# Make sure each cov has mean 0 and sd 1:
expression_matrix = expression_matrix - rowMeans(expression_matrix)
expression_matrix = expression_matrix / rowSds(as.matrix(expression_matrix))

# Import MatrixEQTL
library(MatrixEQTL)

gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(expression_matrix));

cvrt = SlicedData$new()

all_n_tests_trans = 0
all_hist_vals_trans = rep(0,100)

trans_eqtl_placeholder = vector('list', 10)

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1e-2

for (i in c(1:10)) {
  print(i)
  # Only for all-by-all runs - divide the chromosomes into ten because the jobs because there are a lot of tests
  snp_set_size = ceiling(dim(genotype_matrix)[1]/10)
  if (i==10) {
    genotype_matrix_part = genotype_matrix[(((i-1)*snp_set_size+1):dim(genotype_matrix)[1]),]
    snp_positions_part = snp_positions[(((i-1)*snp_set_size+1):dim(snp_positions)[1]),]
  } else {
    genotype_matrix_part = genotype_matrix[(((i-1)*snp_set_size+1):((i)*snp_set_size)),]
    snp_positions_part = snp_positions[(((i-1)*snp_set_size+1):((i)*snp_set_size)),]
  }

  snps = SlicedData$new()
  snps$CreateFromMatrix(as.matrix(genotype_matrix_part));

  useModel = modelLINEAR;

  # Only for generating the P-value hist:
  me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = NULL,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      verbose = TRUE,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = pvOutputThreshold,
      snpspos = snp_positions_part,
      genepos = gene_positions,
      cisDist = cis_dist,
      pvalue.hist = 100,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)

  # record trans- pairs
  if (dim(me$trans$eqtls)[1] > 0) {
    snps_list = sapply(me$trans$eqtls$snps, as.character)
    snp_chr = snp_positions[snps_list,]$chr
    gene_list = sapply(me$trans$eqtls$gene, as.character)
    snp_loc = snp_positions[snps_list,]$pos

    temp_df = me$trans$eqtls
    temp_df$snp_chr = snp_chr
    temp_df$snp_pos = snp_loc

    trans_eqtl_placeholder[[i]] = temp_df
  }

  all_n_tests_trans = all_n_tests_trans + me$trans$ntests
  all_hist_vals_trans = all_hist_vals_trans + me$trans$hist.counts
  rm(me)
}

cumulative_me_null = data.frame(snps = character(0), gene = character(0), statistic = numeric(0), pvalue = numeric(0), FDR = numeric(0), beta = numeric(0), snp_chr = character(0), snp_pos = numeric(0))

cumulative_me_trans = do.call("rbind", trans_eqtl_placeholder)
if (!is.null(cumulative_me_trans)) {
  ordering = order(cumulative_me_trans$pvalue)
  ordered_pval = cumulative_me_trans$pvalue[ordering]
  cumulative_me_trans = cumulative_me_trans[ordering,]
} else {
  cumulative_me_trans = cumulative_me_null
}

cumulative_me_trans$tissue_MAF = rowSums(genotype_matrix[as.character(cumulative_me_trans$snps),])/dim(genotype_matrix)[2]/2
cumulative_me_trans$tissue_MAF = sapply(cumulative_me_trans$tissue_MAF, function(x) {min(x, 1-x)})

# Save the Matrix EQTL object as RData
save(cumulative_me_trans, all_n_tests_trans, all_hist_vals_trans, gene_set_size, snp_set_size_init, snp_set_size_final, file=paste(out_file, '_part', chr_number, '_me.RData', sep="") )
