##################################################
#  trans_matrix_eqtl_null_geneFDR.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/trans_matrix_eqtl_null_geneFDR.R
# 
#  This version performs the all-by-all trans-mapping analysis for null genes, and saves only the best eSNP.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# # Example
# args = c(1:7)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/skinnotsunexposedsuprapubic_v6p_consortium_autosomes_normalized.txt'
# args[2] = '1'
# args[3] = '2.5e8'
# args[4] = 'skinnotsunexposedsuprapubic'
# args[5] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
# args[6] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/skinnotsunexposedsuprapubic/skinnotsunexposedsuprapubic_nonverlapping_certain_autosomes_normalized_MatrixEQTL'
# # Attach PEER number at the end
# args[7] = '10'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
geno_option = 'continuous'
part_number = as.numeric(args[2])
cis_dist = as.numeric(args[3])
tissue_name = args[4]
cov_dir = args[5]
out_file = args[6]
num_split = as.numeric(args[7])
chr_number = ceiling(part_number/num_split)

# Take all SNPs for the trans runs
# genotype_file_name = paste(args[8], 'Genotype/GTEx/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_', disc_or_cont , '_Chr', chr_number, '_Final.txt', sep="")
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr_number, '_Final.txt', sep="")
SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', chr_number, '.txt', sep="")
gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt')

# Need to check that the subjects read in from expression matrices match the subjects in genotype data
library(dplyr)

expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]

# Load in the mappability list:
mappability_cutoff = 0.8
mappability_list = read.table(paste0(proj_dir, '/Data/Resources/annotations/avg_mappability_Exon_UTR.txt'), stringsAsFactors=F)
rownames(mappability_list) = mappability_list$V1
# Arbitrary 0.8 cutoff - can be modified for a different threshold
mappability_list = mappability_list[(mappability_list$V2>mappability_cutoff),]
# Filter out genes with low mappability
expression_matrix = expression_matrix[rownames(expression_matrix) %in% rownames(mappability_list),]

genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

# Only for all-by-all runs - divide the chromosomes into ten because the jobs because there are a lot of tests
snp_set_size = ceiling(dim(genotype_matrix)[1]/num_split)
multiplier = part_number - num_split*(chr_number-1)
if (multiplier==num_split) {
  genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):dim(genotype_matrix)[1]),]
} else {
  genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):((multiplier)*snp_set_size)),]
}

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
gene_positions = gene_positions[gene_positions$gene_id %in% rownames(expression_matrix),]
rownames(gene_positions) = gene_positions$gene_id
gene_positions = gene_positions[rownames(expression_matrix),]
gene_positions$start = as.numeric(gene_positions$start)
gene_positions$end = as.numeric(gene_positions$end)
# gene_positions = gene_positions[sapply(rownames(expression_matrix), function(x) {match(x, gene_positions$gene_id)}),]
# rownames(gene_positions) = gene_positions$gene_id

# Incorporate a tissue-specific MAF filter so that we are not testing variants that actually have a low MAF in our test set
gene_set_size = dim(expression_matrix)[1]
snp_set_size_init = dim(genotype_matrix)[1]

MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
genotype_matrix = genotype_matrix[MAF>0.05 & MAF<0.95,]
snp_positions = snp_positions[MAF>0.05 & MAF<0.95,]

snp_set_size_final = dim(genotype_matrix)[1]

# Also prepare permuted matrix
# geno_option = 'continuous_permute'
# genotype_permute_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr_number, '_Final.txt', sep="")
# genotype_matrix_permute = read.table(genotype_permute_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
# genotype_matrix_permute = genotype_matrix_permute[rownames(genotype_matrix),colnames(genotype_matrix)]

# Load in the covariates
suffix = '_Analysis.covariates.txt'
covars = read.csv(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

# Import MatrixEQTL
library(MatrixEQTL)
source(paste0(proj_dir, '/Scripts/eqtls/MatrixEQTL_wrapper.R'))

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1e-4

# Permute the matrix:
expression_matrix_permute = expression_matrix
n = 111
for (i in c(1:nrow(expression_matrix))) {
	set.seed(n*i)
	order = sample(1:ncol(expression_matrix))
	expression_matrix_permute[i,] = expression_matrix[i,order]
}

me_permute = MatrixEQTL_wrapper(genotype_matrix, expression_matrix_permute, snp_positions, gene_positions, pvThresh = pvOutputThreshold, pvThresh_cis = 0, covariates = covars)

# select the strongest SNP
if (nrow(me_permute$trans$eqtls) > 0) {
	me_permute$trans$eqtls = me_permute$trans$eqtls[!duplicated(me_permute$trans$eqtls$gene),]
}

# From the original expression: 
genes_represented_file = paste0(proj_dir, '/Output/trans-mapping/gtex/MatrixEQTL/gene_level_FDR/genes_represented/', tissue_name, '_genes_present.txt')
genes_represented = read.table(genes_represented_file, header=T)

total_list = c(as.character(genes_represented$gene), rownames(expression_matrix))
unrepresented = (!duplicated(total_list))[(nrow(genes_represented)+1):length(total_list)]

if (sum(unrepresented) > 0) {
	expression_matrix_part = expression_matrix[unrepresented,]
	gene_positions_part = gene_positions[unrepresented,]
	me = MatrixEQTL_wrapper(genotype_matrix, expression_matrix_part, snp_positions, gene_positions_part, pvThresh = pvOutputThreshold, pvThresh_cis = 0, covariates = covars)

	# select the strongest SNP
	if (nrow(me$trans$eqtls) > 0) {
		me$trans$eqtls = me$trans$eqtls[!duplicated(me$trans$eqtls$gene),]
	}

	save(me, me_permute, gene_set_size, snp_set_size_init, snp_set_size_final, file=paste0(out_file, '_part', sprintf("%03d", part_number), '.RData') )
} else {
	save(me_permute, gene_set_size, snp_set_size_init, snp_set_size_final, file=paste0(out_file, '_part', sprintf("%03d", part_number), '.RData') )
}
