##################################################
#  GGR_validation_specific_pairs.R
#
#  $proj/Scripts/causality/GGR/gtex/GGR_validation_specific_pairs.R
# 
#  GGR validation pipeline with trans-eQTLs and MR implementations.
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

graph_dir = '/tigress/BEE/RNAseq/Analysis/Figures/GGR_validation_plots/eqtl_cor_plots/'
library(ggplot2)

plot_function_cor = function(plot_pair, covars, gene_pair_type = '') {
	chr = strsplit(plot_pair$snp_chr, 'chr')[[1]][2]
	genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr, '_Final.txt', sep="")
	genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	genotype_matrix = genotype_matrix[,colnames(expression_matrix)]
	x = genotype_matrix[plot_pair$snps,]
	y = expression_matrix[plot_pair$cis_gene,]
	z = expression_matrix[plot_pair$gene,]

	X = t(as.matrix(covars) - rowMeans(covars))
	geno_center = x - mean(as.numeric(x))

	res_y = y - t(X %*% solve(t(X) %*% X) %*% t(X) %*% as.numeric(y))
	res_z = z - t(X %*% solve(t(X) %*% X) %*% t(X) %*% as.numeric(z))
	res_x = geno_center - t(X %*% solve(t(X) %*% X) %*% t(X) %*% as.numeric(geno_center))

	plot_df = data.frame(cis = as.numeric(y), trans = as.numeric(z), genotype = as.numeric(round(x)))
	g = ggplot(plot_df, aes(cis, trans, col = factor(genotype))) + geom_point() + xlab(plot_pair$cis_gene) + ylab(plot_pair$gene) + theme_classic()
	ggsave(filename = paste0(graph_dir, 'expression_cor/', condition, '_', regression, '_', gene_pair_type, '_expression_cor.pdf'), plot = g)
	g = ggplot(plot_df) + geom_boxplot(aes(x = factor(genotype), y=trans), outlier.shape=NA) + geom_jitter(data = plot_df, aes(x = genotype+1, y = trans), width = 0.2) + xlab(plot_pair$snps) + ylab(plot_pair$gene) + theme_classic()
	ggsave(filename = paste0(graph_dir, 'trans_eqtl/', condition, '_', regression, '_', gene_pair_type, '_trans_eqtl.pdf'), plot = g)

	plot_df_res = data.frame(cis = as.numeric(res_y), trans = as.numeric(res_z), genotype = as.numeric(round(res_x)))
	g = ggplot(plot_df_res, aes(cis, trans, col = factor(genotype))) + geom_point() + xlab(plot_pair$cis_gene) + ylab(plot_pair$gene) + theme_classic()
	ggsave(filename = paste0(graph_dir, 'expression_cor/', condition, '_', regression, '_', gene_pair_type, '_expression_cor_corrected.pdf'), plot = g)
	g = ggplot(plot_df_res) + geom_boxplot(aes(x = factor(genotype), y=trans), outlier.shape=NA) + geom_jitter(data = plot_df_res, aes(x = genotype+1, y = trans), width = 0.2) + xlab(plot_pair$snps) + ylab(plot_pair$gene) + theme_classic()
	ggsave(filename = paste0(graph_dir, 'trans_eqtl/', condition, '_', regression, '_', gene_pair_type, '_trans_eqtl_corrected.pdf'), plot = g)

}

# FDR Control of p-values from the null distribution
# condition = '0mean-unnormalized_g-null_l-fdr'
# original_dir = paste0('/tigress/BEE/RNAseq/Output/causality/GGR/output/', condition, '/')
# regression = 'enet-2'
# write_dir = paste0('/tigress/BEE/RNAseq/Output/causality/GGR/final_list/', condition, '/')

trans_eqtl_dir = paste0(proj_dir, '/Output/causality/GGR/final_list/all_by_all/')
gene_pair_dir = paste0(proj_dir, '/Data/Networks/GGR/priority/analysis/tf-interesting-edge-networks/')

# regressions = c('enet-1', 'enet-2')
conditions = c('0mean-1var', '0mean-unnormalized')
regressions = c('enet-1', 'enet-2')
gene_pair_types = c('tf-imm', 'tf-metab', 'tf-tf')

# Load in expression values
expression_file_location = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/lung_v6p_consortium_autosomes_normalized.txt'
expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
rownames(expression_matrix) = sapply(rownames(expression_matrix), function(x) {strsplit(x, '\\.')[[1]][1]})

# Load in the covariates
suffix = '_Analysis.covariates.txt'
cov_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
covars = read.csv(paste0(cov_dir, 'lung', suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

for (condition in conditions) {

	for (regression in regressions) {
		trans_eqtl_table = read.table(paste0(trans_eqtl_dir, condition, '_g-null_l-fdr/lung_', regression, '_trans_eQTL_table.txt'), header=T, stringsAsFactors=F, sep='\t')
		trans_eqtl_table$identifier = paste(trans_eqtl_table$cis_gene, trans_eqtl_table$gene, sep='_')
		# print(condition)
		# print(regression)
		# print(trans_eqtl_table[1,])
		plot_pair = trans_eqtl_table[1,]
		plot_function_cor(plot_pair, covars)

		gene_pair_file = paste0('/tigress/BEE/RNAseq/Data/Networks/GGR/retrofitted/', condition, '_g-null_l-fdr/prot2TPM-er-reps_', condition, '_g-null_l-fdr-0.05_', regression, '-union-network.txt')
		input_gene_pair = read.csv(file = gene_pair_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

		for (gene_pair_type in gene_pair_types) {
			gene_pair_table = read.table(paste0(gene_pair_dir, gene_pair_type, '/', condition, '_', regression, '-', gene_pair_type, '.txt'), header=T, stringsAsFactors=F, sep='\t')
			gene_pair_table$Cause_v19 = as.character(sapply(gene_pair_table$Cause, function(x) {input_gene_pair$Cause_v19[which(x == input_gene_pair$Cause)][1]}))
			gene_pair_table$Effect_v19 = as.character(sapply(gene_pair_table$Effect, function(x) {input_gene_pair$Effect_v19[which(x == input_gene_pair$Effect)][1]}))
			rownames(gene_pair_table) = paste(gene_pair_table$Cause_v19, gene_pair_table$Effect_v19, sep='_')
			print(condition)
			print(regression)
			print(gene_pair_type)
			plot_pair = trans_eqtl_table[which(sapply(trans_eqtl_table$identifier, function(x) {x %in% rownames(gene_pair_table)}))[1],]
			print(plot_pair)
			plot_function_cor(plot_pair, covars, gene_pair_type = gene_pair_type)
		}
	}
}


# trans_eqtl_dir = paste0(proj_dir, '/Output/causality/GGR/final_list/best_cis_only/')

# for (condition in conditions) {
# 	for (regression in regressions) {
# 		trans_eqtl_table = read.table(paste0(trans_eqtl_dir, condition, '_g-null_l-fdr/lung_', regression, '_trans_eQTL_table.txt'), header=T, stringsAsFactors=F, sep='\t')

# 		for (gene_pair_type in gene_pair_types) {
# 			gene_pair_table = read.table(paste0(gene_pair_dir, gene_pair_type, '/', condition, '_', regression, '-', gene_pair_type, '.txt'), header=T, stringsAsFactors=F, sep='\t')

# 			print(condition)
# 			print(regression)
# 			print(gene_pair_type)
# 			print(trans_eqtl_table[head(which(sapply(trans_eqtl_table$Cause.Effect, function(x) {x %in% gene_pair_table$Cause.Effect}))),])
# 		}
# 	}
# }

library(ggplot2)
pdf("rs1867277_tmem253.pdf", useDingbats = F)
ggplot(datplot) + 
geom_boxplot(aes(x = as.factor(round(rs1867277)), y=expression.subset), outlier.shape=NA)+
#geom_jitter(data = datplot[which(!datplot$rs1867277 %in% c(0,1,2)),], aes(x = rs1867277+1, y = expression.subset), color = "gray30") +
geom_jitter(data = datplot, aes(x = rs1867277+1, y = expression.subset), color = "gray30", width = 0.2) +
theme_classic() + theme(legend.position="none")
dev.off()
