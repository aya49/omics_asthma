##################################################
#  MR_PPA_simulation_prep.R
#
#  $proj/Scripts/causality/bayes_MR/main_pipeline/MR_PPA_simulation_prep.R
# 
#  Prepares the simulation expression dataset and calculates MR stats
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
# args = c(1:3)
# args[1] = 'Whole_Blood'
# args[2] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/prep_files/'
# args[3] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/MR_stats/'

tissue_name = args[1]
# Make sure that the chromosome and part numbers are the same as the eQTL run
prep_dir = args[2]
MR_stats_dir = paste0(args[3], args[1], '/')

if (!file.exists(paste0(prep_dir, tissue_name, '_simulation_prep_data.RData'))) {

load(paste0(prep_dir, tissue_name, '.RData'))

MR_result_files = list.files(MR_stats_dir)
MR_stats_comb = list()
counter = 1
for (f in MR_result_files) {
	load(paste0(MR_stats_dir, f))
	MR_stats_comb[[counter]] = MR_stats
	counter = counter + 1
}

# Which gene pairs had cis-trans causal relationships?
MR_stats_comb = do.call('rbind', MR_stats_comb)
gene_pairs = data.frame(cis_gene = MR_stats_comb$cis_gene, trans_gene = MR_stats_comb$trans_gene, gene_pair = paste0(MR_stats_comb$cis_gene, '_', MR_stats_comb$trans_gene), stringsAsFactors = F)
gene_pairs = gene_pairs[!duplicated(gene_pairs$gene_pair),]

# Randomly select gene set while making sure not to include any cis-trans causal relationships
gene_set = c()

while (length(gene_set) <= 2000) {
	set.seed(42)
	gene_set = c(gene_set, sample(setdiff(rownames(expression_matrix), gene_set), size = 500, replace = F))
	trans_targets = lapply(gene_set, function(x) {gene_pairs[which(gene_pairs$cis_gene == x),]$trans_gene})
	trans_targets = unique(unlist(trans_targets))
	gene_set = setdiff(gene_set, trans_targets)
	print(length(gene_set))
}

# Final gene set
gene_set = sample(gene_set, size = 2000, replace = F)

expression_matrix = expression_matrix[gene_set,]
gene_positions = gene_positions[gene_set,]
# Select gene pairs that will receive a causal signal

set.seed(111)
src = sample(rownames(expression_matrix), size = 100000, replace = T)
set.seed(1111)
dest = sample(rownames(expression_matrix), size = 100000, replace = T)
gene_pair_list = data.frame(src = src, dest = dest, stringsAsFactors = F)
# Make sure we don't have self- and reverse-causation events
gene_pair_list = gene_pair_list[(gene_pair_list$src != gene_pair_list$dest),]
gene_pair_list$id = sapply(c(1:nrow(gene_pair_list)), function(x) {paste0(gene_pair_list[x,][order(gene_pair_list[x,])[1]], '_', gene_pair_list[x,][order(gene_pair_list[x,])[2]])})
gene_pair_list = gene_pair_list[!duplicated(gene_pair_list$id),]
gene_pair_list = gene_pair_list[c(1:50000),]

mult_mat = matrix(0, nrow(expression_matrix), nrow(expression_matrix))
rownames(mult_mat) = rownames(expression_matrix)
colnames(mult_mat) = rownames(expression_matrix)
l = ncol(expression_matrix)

for (i in c(1:10000)) {
	if (i%%100 == 0) {print(i)}
	src = gene_pair_list$src[i]
	dest = gene_pair_list$dest[i]
	if (cor(as.numeric(expression_matrix[src,]), as.numeric(expression_matrix[dest,]))>0) { mult_mat[dest,src] = 0.5 } else { mult_mat[dest,src] = -0.5 }
}
for (i in c(10001:20000)) {
	if (i%%100 == 0) {print(i)}
	src = gene_pair_list$src[i]
	dest = gene_pair_list$dest[i]
	if (cor(as.numeric(expression_matrix[src,]), as.numeric(expression_matrix[dest,]))>0) { mult_mat[dest,src] = 0.4 } else { mult_mat[dest,src] = -0.4 }
}
for (i in c(20001:30000)) {
	if (i%%100 == 0) {print(i)}
	src = gene_pair_list$src[i]
	dest = gene_pair_list$dest[i]
	if (cor(as.numeric(expression_matrix[src,]), as.numeric(expression_matrix[dest,]))>0) { mult_mat[dest,src] = 0.3 } else { mult_mat[dest,src] = -0.3 }
}
for (i in c(30001:40000)) {
	if (i%%100 == 0) {print(i)}
	src = gene_pair_list$src[i]
	dest = gene_pair_list$dest[i]
	if (cor(as.numeric(expression_matrix[src,]), as.numeric(expression_matrix[dest,]))>0) { mult_mat[dest,src] = 0.2 } else { mult_mat[dest,src] = -0.2 }
}
for (i in c(40001:50000)) {
	if (i%%100 == 0) {print(i)}
	src = gene_pair_list$src[i]
	dest = gene_pair_list$dest[i]
	if (cor(as.numeric(expression_matrix[src,]), as.numeric(expression_matrix[dest,]))>0) { mult_mat[dest,src] = 0.1 } else { mult_mat[dest,src] = -0.1 }
}

expression_matrix_sim = expression_matrix + (mult_mat %*% as.matrix(expression_matrix))

# Add some noise and re-quantile-normalize
expression_matrix_sim = expression_matrix_sim + do.call('rbind', lapply(c(1:nrow(expression_matrix)), function(x) {rnorm(ncol(expression_matrix), 0, 0.1)}))
expression_matrix_sim = t(apply(expression_matrix_sim, 1, rank, ties.method = "average"))
expression_matrix_sim = qnorm(expression_matrix_sim / (ncol(expression_matrix)+1))

save(expression_matrix, mult_mat, expression_matrix_sim, gene_pair_list, gene_positions, file = paste0(prep_dir, tissue_name, '_simulation_prep_data.RData'))

}