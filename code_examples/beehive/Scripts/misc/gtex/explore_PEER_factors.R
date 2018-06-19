##################################################
#  explore_PEER_factors.R
#
#  $proj/Scripts/misc/gtex/explore_PEER_factors.R
# 
#  PEER factor exploartory analysis - currently only calculates the FVU (Fraction of Variance Unexplained)
#  on the expression data. 
#
#  Author: Brian Jo
#
##################################################

args = c(1:3)
args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
args[2] = '_nonverlapping_certain_autosomes_normalized.txt'
args[3] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'

exp_dir = args[1]
suffix = args[2]
cov_dir = args[3]

exp_mats = list.files(exp_dir, pattern = suffix)

for (item in exp_mats) {
	tissue = strsplit(item, '_')[[1]][1]
	print(tissue)
	
	expression_matrix = read.csv(file = paste0(exp_dir, item), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	rownames(expression_matrix) = expression_matrix$gene_id
	expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]

	covars = read.table(paste0(cov_dir, tissue, '_Analysis.covariates.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	rownames(covars) = covars$ID
	covars = covars[,colnames(expression_matrix)]
	if ('gender' %in% rownames(covars)) {
		N_peer = dim(covars)[1] - 5
		initial_list = c('C1', 'C2', 'C3', 'gender', 'Platform')
	} else {
		N_peer = dim(covars)[1] - 4
		initial_list = c('C1', 'C2', 'C3', 'Platform')
	}
	rem_list = rownames(covars)[!(rownames(covars) %in% initial_list)]

	# expr.svd = svd(t(expression_matrix))
	# print(diag(cor(expr.svd$u[,1:N_peer], t(covars[4:(N_peer+3),]), method = 'spearman')))

	expression_matrix = expression_matrix - rowMeans(expression_matrix)

	calc_FVU = function(y, X) {
		y_hat = (X %*% solve(t(X) %*% X) %*% t(X) %*% y)
		return(1 - (sum(y_hat*y_hat) / sum(y*y)))
	}
	FVU = calc_FVU(as.matrix(t(expression_matrix)), as.matrix(t(covars[initial_list,])))

	var_unexp_df = c('init', FVU)

	for (i in c(1:N_peer)) {
		print(i)
		next_FVU = as.numeric(sapply(rem_list, function(x) {calc_FVU(as.matrix(t(expression_matrix)), as.matrix(t(covars[c(initial_list, x),])))}))
		min = which(next_FVU == min(next_FVU))

		print(rem_list[min])
		print(next_FVU[min])
		var_unexp_df = rbind(var_unexp_df, c(rem_list[min], next_FVU[min]))
		initial_list = c(initial_list, rem_list[min])
		rem_list = rownames(covars)[!(rownames(covars) %in% initial_list)]
	}

}