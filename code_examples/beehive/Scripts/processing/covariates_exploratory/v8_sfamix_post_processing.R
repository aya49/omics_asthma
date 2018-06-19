##################################################
#  v8_sfamix_post_processing.R
#
#  $proj/Scripts/processing/covariates_exploratory/v8_sfamix_post_processing.R
# 
#  This script calculates the MR stats (both frequentist and Bayesian) for SNP-cis-trans trios for the GTEx v8 data.
#
#  Author: Brian Jo
#
##################################################

sfamix_out_dir = '/tigress/BEE/RNAseq/Output/processing/exploratory/v8/sfamix/GTEx_Analysis_v8/'
tissue_list = list.files(path = sfamix_out_dir)

n_factor_df = list()

tissue_list = c('Adipose_Subcutaneous')

for (tissue in tissue_list) {
	n_steps = sapply(list.files(path = paste0(sfamix_out_dir, tissue), pattern = 'Z_'), function(x) {strsplit(x, '_')[[1]][2]})
	n_steps = as.numeric(n_steps)[order(as.numeric(n_steps))]

	n_factors = sapply(n_steps[2:length(n_steps)], function(steps) {
		length(as.numeric(read.csv(file = paste0(sfamix_out_dir, tissue, '/Z_', steps), sep = ' ', header = F)[1,]))
	})
	n_dense_factors = sapply(n_steps[2:length(n_steps)], function(steps) {
		sum(as.numeric(read.csv(file = paste0(sfamix_out_dir, tissue, '/Z_', steps), sep = ' ', header = F)[1,]))
	})
	n_factor_df[[tissue]] = data.frame(n_steps = n_steps[2:length(n_steps)], n_factors = n_factors, n_dense_factors = n_dense_factors, tissue = tissue)

	peer_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
	PEER = read.csv(paste0(peer_dir, tissue, '.v8.covariates.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	rownames(PEER) = PEER$ID
	PEER = PEER[,c(2:ncol(PEER))]

	exp_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
	suffix = '.v8.normalized_expression.bed.gz'
	expression_file_location = paste0(exp_dir, tissue, '.v8.normalized_expression.bed.gz')
	header = readLines(gzfile(expression_file_location), n = 1)
	header = strsplit(header, '\t')[[1]]

	for (steps in n_steps[2:length(n_steps)]) {
		# Read in LAM and PEER for comparison:
		LAM = read.table(paste0(sfamix_out_dir, tissue, '/LAM_', steps))
		r = cor(LAM, t(PEER))



		# If there are any sparse factors, read in EX for the gene list

	}
}

n_factor_df = do.call('rbind', n_factor_df)
save(n_factor_df, file = '/tigress/BEE/RNAseq/Analysis/Figures/exploratory/sfamix/n_factors.RData')

library(ggplot2)

n_factor_df = n_factor_df[n_factor_df$n_factors <= 500,]
ggplot(n_factor_df, aes(x = n_steps, y = n_factors, col = tissue)) + geom_line() + guides(col = F)

n_factor_df$n_sparse_factors = n_factor_df$n_factors - n_factor_df$n_dense_factors

ggplot(n_factor_df, aes(x = n_steps, y = n_sparse_factors, col = tissue)) + geom_line() + guides(col = F)
