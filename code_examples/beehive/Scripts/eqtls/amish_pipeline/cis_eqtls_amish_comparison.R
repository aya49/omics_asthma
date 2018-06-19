##################################################
#  cis_eqtls_amish_comparison.R
#
#  $proj/Scripts/eqtls/amish_pipeline/cis_eqtls_amish_comparison.R
#
#  Script for generating comparison stats for Amish and GTEx cis-eQTLs
#
#  Author: Brian Jo
#
##################################################

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')
library(dplyr)

# args = c(1:8)
# args[1] = '/tigress/BEE/amish/analyses/ciseqtl/genomewide/cis_eqtls_1mb_chr22.txt'
# args[2] = '22'
# args[3] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
# args[4] = 'cellsebvtransformedlymphocytes_v6p_consortium_autosomes_normalized.txt'
# args[5] = 'v6p'
# args[6] = '/tigress/BEE/RNAseq/Data/Genotype/gtex/amish_comparison/chr22.txt'
# args[7] = 'cellsebvtransformedlymphocytes'
# args[8] = '/tigress/BEE/RNAseq/Output/cis-mapping/amish/'

amish_cis_eqtl_f = args[1]
chrom = args[2]
expression_dir = args[3]
expression_file = args[4]
version = args[5]
genotype_matrix_f = args[6]
tissue_name = args[7]
out_dir = args[8]

amish_cis_eqtls = read.table(amish_cis_eqtl_f, sep = '\t', header = T, stringsAsFactors = F)
snp_list = unique(amish_cis_eqtls$rs)
print(length(snp_list))
gene_list = unique(amish_cis_eqtls$gene_id)
print(length(gene_list))

expression_matrix = read.csv(file = paste0(expression_dir, expression_file), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]

genotype_matrix = read.table(genotype_matrix_f, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt')
# Also obtain genes and snps that are in the expression and genotype files, respectively.
gene_positions = read.table(file = gene_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gene_positions = gene_positions[,c('gene_id', 'chr', 'start', 'end')]
gene_positions = mutate(gene_positions, chr = paste("chr", chr, sep = ""))
gene_positions = select(gene_positions, c(gene_id, chr, start, end))
rownames(gene_positions) = gene_positions$gene_id

# Leave the genes that are in the amish files only
inds = sapply(gene_list, function(x) {x %in% rownames(expression_matrix)})
expression_matrix = expression_matrix[gene_list[inds],]
gene_positions = gene_positions[gene_list[inds],]

# Also obtain the data frame with SNP positions
snp_positions = data.frame(rsID = rownames(genotype_matrix), chr = paste0('chr', chrom), pos = sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][2]}), stringsAsFactors = F)
snp_positions$pos = as.numeric(snp_positions$pos)
MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
MAF = sapply(MAF, function(x) {min(x, 1-x)})

# Load in the covariates
suffix = '_Analysis.covariates.txt'
covars = read.csv(paste0(expression_dir, 'covariates/', tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

# Import MatrixEQTL
library(MatrixEQTL)
source(paste0(proj_dir, '/Scripts/eqtls/MatrixEQTL_wrapper.R'))

# Only interested in cis-eQTLs - within 1mb
set.seed(42)
genotype_matrix_perm = genotype_matrix[,sample(colnames(genotype_matrix))]
colnames(genotype_matrix_perm) = colnames(genotype_matrix)
me = MatrixEQTL_wrapper(genotype_matrix, expression_matrix, snp_positions, gene_positions, pvThresh = 0, pvThresh_cis = 1, covariates = covars, cis_dist = 1000000)
me_perm = MatrixEQTL_wrapper(genotype_matrix_perm, expression_matrix, snp_positions, gene_positions, pvThresh = 0, pvThresh_cis = 1, covariates = covars, cis_dist = 1000000)

# construct comparison data frame
rownames(amish_cis_eqtls) = paste(amish_cis_eqtls$gene_id, amish_cis_eqtls$rs, sep='_')
comp_df = me$cis$eqtls
comp_df$MAF = MAF[comp_df$snps]
rownames(comp_df) = paste(comp_df$gene, comp_df$snps, sep='_')

comp_df_perm = me_perm$cis$eqtls
# comp_df_perm$MAF = MAF[comp_df_perm$snps]
rownames(comp_df_perm) = paste(comp_df_perm$gene, comp_df_perm$snps, sep='_')
comp_df_perm = comp_df_perm[,c('statistic', 'pvalue', 'beta')]
colnames(comp_df_perm) = c('statistic_permuted', 'pvalue_permuted', 'beta_permuted')

out_df = cbind(comp_df, comp_df_perm[rownames(comp_df),])

out_df = cbind(out_df, amish_cis_eqtls[rownames(out_df),c('af', 'allele1', 'beta','se','l_remle','l_mle','p_wald','p_lrt')])
colnames(out_df)[c(16:19)] = c('amish_af', 'amish_allele1', 'amish_beta', 'amish_se')

write.table(out_df, quote = F, sep = '\t', row.names = F, file = paste0(out_dir, tissue_name, '/gtex_amish_comparison_ciseqtls_chr', chrom, '.txt'))
