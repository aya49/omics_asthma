##################################################
#  cis_trans_eqtl_MR_analysis.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/cis_trans_eqtl_MR_analysis.R
# 
#  Running MR stats on cis-trans-eQTLs
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

library(MatrixEQTL)
source('/tigress/BEE/RNAseq/Scripts/eqtls/MatrixEQTL_wrapper.R')

# Testing - how many trans are also cis?

cis_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/cis_eqtls_from_Stanford/'
trans_list = read.table('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/trans-eQTLs_FDR-0.1.txt', header=T, stringsAsFactors=F, sep='\t')

# sum = 0
# for (tissue in unique(trans_list$tissue)) {
# 	print(tissue)
# 	tissue_subset = trans_list[trans_list$tissue == tissue,]
# 	cis_list = read.table(paste0(cis_dir, tissue, '_eqtl_list.txt'), header=T, stringsAsFactors=F, sep='\t')
# 	in_cis = which(sapply(tissue_subset$snps, function(x) {x %in% cis_list$rsID}))
# 	sum = sum + length(in_cis)
# 	print(length(in_cis))

# 	# Running at pval cutoff of 1e-5 - 290
# }

calculate_t_MR = function(z, x, y, beta_xz, beta_yz) {
  beta_mr = beta_yz/beta_xz
  sig_sq = (y - x*beta_mr) %*% (y - x*beta_mr) / (length(z) - 3)
  var_beta = sig_sq * (z %*% z) / (x %*% z)^2
  return(beta_mr*beta_mr/var_beta)
}

return_row = function(snp, cis_gene, trans_gene, expression_matrix, genotype_matrix, beta_xz, beta_yz, covars, cis_pval, trans_pval) {
	genotype_matrix_part = genotype_matrix[snp,]
	expression_matrix_part = expression_matrix[c(cis_gene, trans_gene),]

    # Regress out covariates - Need this step for MR
    X = t(as.matrix(covars) - rowMeans(covars))
    y = t(as.matrix(expression_matrix_part) - rowMeans(expression_matrix_part))
    g = t(as.matrix(genotype_matrix_part) - rowMeans(genotype_matrix_part))

    res_y = y - (X %*% solve(t(X) %*% X) %*% t(X) %*% y)
    res_g = g - (X %*% solve(t(X) %*% X) %*% t(X) %*% g)

    geno_res = data.frame(t(res_g))
    exp_res = data.frame(t(res_y))

    t_MR = calculate_t_MR(as.numeric(geno_res[snp,]), as.numeric(exp_res[cis_gene,]), as.numeric(exp_res[trans_gene,]), beta_xz, beta_yz)

    return(data.frame(snp = snp, cis_gene = cis_gene, trans_gene = trans_gene, beta_xz = beta_xz, beta_yz = beta_yz, t_MR = t_MR, cis_pval = cis_pval, trans_pval = trans_pval))
}

t_MR_counter = 1
t_MR_counter_permute = 1
t_MR_list = list()
t_MR_permute_list = list()

for (tissue in c('brainputamenbasalganglia','cellstransformedfibroblasts','muscleskeletal','pancreas','testis','thyroid')) {
	# Run MR
	print(tissue)
	tissue_subset = trans_list[trans_list$tissue == tissue,]
	cis_list = read.table(paste0(cis_dir, tissue, '_eqtl_list.txt'), header=T, stringsAsFactors=F, sep='\t')
	in_cis = which(sapply(tissue_subset$snps, function(x) {x %in% cis_list$rsID}))

	prev_chr = ''
	if (length(in_cis) > 0) {
		test_set = tissue_subset[in_cis,]
		test_set = test_set[order(test_set$snp_chr),]

		# load expression matrix
		expression_file_location = paste0('/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/', tissue, '_v6p_consortium_autosomes_normalized.txt')
		expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
		rownames(expression_matrix) = expression_matrix$gene_id
		expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
		# load covariate matrix
		covars = read.csv(paste0('/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/', tissue, '_Analysis.covariates.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
		rownames(covars) = covars$ID
		covars = covars[,colnames(expression_matrix)]

		print(nrow(test_set))
		for (i in c(1:nrow(test_set))) {
			print(i)
			current_chr = strsplit(test_set$snp_chr[i], 'chr')[[1]][2]
			if (prev_chr != current_chr) {
				prev_chr = current_chr
				genotype_file_name = paste0('/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', current_chr, '_Final.txt', sep="")
				genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
				genotype_matrix = genotype_matrix[,colnames(expression_matrix)]
			}

			inds = which(test_set[i,]$snps == cis_list$rsID)
			cis_list_part = cis_list[inds,]
			for (j in c(1:length(inds))) {
				cis_gene = cis_list_part[j,]$gene_id
				if (cis_gene %in% rownames(expression_matrix)) {
					beta_xz = cis_list_part[j,]$slope
					cis_pval = cis_list_part[j,]$pval_nom
					t_MR_list[[t_MR_counter]] = return_row(test_set[i,]$snps, cis_gene, test_set[i,]$gene, expression_matrix, genotype_matrix, beta_xz, test_set[i,]$beta, covars, cis_pval, test_set[i,]$pvalue)
					t_MR_list[[t_MR_counter]]$tissue = tissue

					snp_pos = data.frame(rsID = test_set[i,]$snps, chr = 'chr1', pos = 1)
					expression_matrix_permute = expression_matrix[c(1:100),]
					gene_pos = data.frame(gene_id = paste(test_set[i,]$gene, c(1:100), sep='_'), chr = 'chr2', start = 1, end = 2)

					for (k in c(1:100)) {
						expression_matrix_permute[k,] = expression_matrix[test_set[i,]$gene, sample(c(1:ncol(expression_matrix)), replace = FALSE)]
					}

					rownames(expression_matrix_permute) = gene_pos$gene_id
					me = MatrixEQTL_wrapper(genotype_matrix[test_set[i,]$snps,], expression_matrix_permute, snp_pos, gene_pos, pvThresh = 1, pvThresh_cis = 0, covariates = covars)
					t_MRs = lapply(c(1:100), function(x) {return_row(test_set[i,]$snps, cis_gene, as.character(me$trans$eqtls$gene[x]), rbind(expression_matrix_permute[as.character(me$trans$eqtls$gene[x]),], expression_matrix[cis_gene,]), genotype_matrix, beta_xz, me$trans$eqtls$beta[x], covars, cis_pval, me$trans$eqtls$pvalue[x])})
					t_MR_permute_list[[t_MR_counter]] = do.call('rbind', t_MRs)
					t_MR_permute_list[[t_MR_counter]]$tissue = tissue
					
					t_MR_counter = t_MR_counter + 1
				}
			}
		}
	}
}

# t_MR_list = do.call('rbind', t_MR_list)
# t_MR_permute_list = do.call('rbind', t_MR_permute_list)

save(t_MR_list, t_MR_permute_list, file = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MR_eQTL/cis_trans_eQTLs/MR_out.RData')

# Number of cis variants

# variant_list = c()
# for (tissue in c('brainputamenbasalganglia','cellstransformedfibroblasts','muscleskeletal','pancreas','testis','thyroid')) {
# 	print(tissue)
# 	cis_list = read.table(paste0(cis_dir, tissue, '_eqtl_list.txt'), header=T, stringsAsFactors=F, sep='\t')
# 	cis_list = cis_list[cis_list$pval_nom < 1e-5,]
# 	variant_list = c(variant_list, unique(cis_list$snp))
# 	variant_list = unique(variant_list)
# }
# 	cis_list = read.table(paste0(cis_dir, tissue, '_eqtl_list.txt'), header=T, stringsAsFactors=F, sep='\t')


# for (file in z) {
# 	print(file)
# 	a = read.table(file, sep='\t', stringsAsFactors=F, header=T)
# 	a = a[a$pval_nom < 1e-5,]
# 	cis_list = c(cis_list, unique(a$snp))
# 	cis_list = unique(cis_list)
# }
