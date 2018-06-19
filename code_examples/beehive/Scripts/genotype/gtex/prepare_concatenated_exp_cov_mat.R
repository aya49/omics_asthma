# Preparing concatenated run:
expression_in_dir = '/tigress/BEE/eQTLs/Data/Expression/Expression_matrices/GTEx/hg19/filtered/normalized_gtex_gcp_v6p/'
expression_in_suffix = '_nonoverlapping_certain_biotype_normalized_GTEx_gct_RNA-SeQCv1.1.8_gene_rpkm.txt'

covar_in_dir = '/tigress/BEE/eQTLs/Data/References/GTEx/v6p_All_Tissues_covariates_FOR_QC_ONLY/'
covar_in_suffix = '_Analysis.covariates.txt'

expression_out_dir = '/tigress/BEE/eQTLs/Data/Expression/Expression_matrices/GTEx/hg19/filtered/normalized_gtex_gcp_v6p_concat/'

exp_mats = list.files(path = expression_in_dir, pattern = expression_in_suffix)
cov_mats = list.files(path = covar_in_dir, pattern = covar_in_suffix)

for (i in c(1:44)) {
	exp_file = exp_mats[i]
	cov_file = cov_mats[i]
	tissue_name = strsplit(exp_file, expression_in_suffix)[[1]][1]
	print(tissue_name)
	# Process concatenated expression matrix
	exp_mat = read.table(paste(expression_in_dir, exp_file, sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')
	colnames(exp_mat) = as.character(sapply(colnames(exp_mat), function(x) {paste(x, '.', tissue_name, sep='')}))
	if (i==1) {
		expression_concat = exp_mat
	} else {
		# First add new columns
		dummy_df = data.frame(matrix(0L, nrow = dim(expression_concat)[1], ncol = dim(exp_mat)[2]))
		rownames(dummy_df) = rownames(expression_concat)
		colnames(dummy_df) = colnames(exp_mat)
		expression_concat = cbind(expression_concat, dummy_df)
		# Then add new rows
		print(sum(!(rownames(exp_mat) %in% rownames(expression_concat))))
		if (sum(!(rownames(exp_mat) %in% rownames(expression_concat))) > 0) {
			dummy_df = data.frame(matrix(0L, nrow = sum(!(rownames(exp_mat) %in% rownames(expression_concat))), ncol = dim(expression_concat)[2]))
			rownames(dummy_df) = rownames(exp_mat)[which(!(rownames(exp_mat) %in% rownames(expression_concat)))]
			colnames(dummy_df) = colnames(expression_concat)
			expression_concat = rbind(expression_concat, dummy_df)
		}
		# Fill in new expression_values
		expression_concat[rownames(exp_mat),colnames(exp_mat)] = exp_mat[rownames(exp_mat),colnames(exp_mat)]
	}

	# Process concatenated covariate matrix
	cov_mat = read.table(paste(covar_in_dir, cov_file, sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')
	rownames(cov_mat) = cov_mat$ID
	cov_mat = cov_mat[,c(2:dim(cov_mat)[2])]
	colnames(cov_mat) = as.character(sapply(colnames(cov_mat), function(x) {paste(x, '.', tissue_name, sep='')}))
	common_cov = c('C1','C2','C3','gender','Platform')
	inds = as.numeric(sapply(common_cov[common_cov %in% rownames(cov_mat)], function(x) {which(rownames(cov_mat) == x)}))
	rem_mat = cov_mat[-inds,]
	rownames(rem_mat) = sapply(c(1:dim(rem_mat)[1]), function(x) {paste(tissue_name, '_', x, sep='')})
	cov_mat = rbind(cov_mat[inds,], rem_mat)
	if (i==1) {
		covar_concat = cov_mat
	} else {
		# First add new columns
		dummy_df = data.frame(matrix(0L, nrow = dim(covar_concat)[1], ncol = dim(cov_mat)[2]))
		rownames(dummy_df) = rownames(covar_concat)
		colnames(dummy_df) = colnames(cov_mat)
		covar_concat = cbind(covar_concat, dummy_df)
		# Then add new rows
		print(sum(!(rownames(cov_mat) %in% rownames(covar_concat))))
		if (sum(!(rownames(cov_mat) %in% rownames(covar_concat))) > 0) {
			dummy_df = data.frame(matrix(0L, nrow = sum(!(rownames(cov_mat) %in% rownames(covar_concat))), ncol = dim(covar_concat)[2]))
			rownames(dummy_df) = rownames(cov_mat)[which(!(rownames(cov_mat) %in% rownames(covar_concat)))]
			colnames(dummy_df) = colnames(covar_concat)
			covar_concat = rbind(covar_concat, dummy_df)
		}
		# Fill in new expression_values
		covar_concat[rownames(cov_mat),colnames(cov_mat)] = cov_mat[rownames(cov_mat),colnames(cov_mat)]
	}
}
#expression_concat[is.na(expression_concat)] = 0

write.table(expression_concat, file=paste(expression_out_dir, 'total', expression_in_suffix, sep=""), quote = FALSE, sep = "\t")
# Also divide into partitions, since the total matrix is large:

num_genes = ceiling(dim(expression_concat)[1] / 10)
for (i in c(1:10)) {
  gene_range = c( ((i-1)*num_genes + 1) : (i*num_genes) )
  if (i == 10) {
    gene_range = c( ((i-1)*num_genes + 1) : (dim(expression_concat)[1]) )
  }
  write.table(expression_concat[gene_range,], file=paste(expression_out_dir, 'total_part', i, expression_in_suffix, sep=""), quote = FALSE, sep = "\t")
  print(i)
}

write.table(covar_concat, file=paste(covar_in_dir, 'total_concatenated_covariates.txt', sep=""), quote = FALSE, sep = "\t")