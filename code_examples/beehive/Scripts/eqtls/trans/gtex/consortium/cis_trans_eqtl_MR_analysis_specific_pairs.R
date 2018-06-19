##################################################
#  cis_trans_eqtl_MR_analysis_specific_pairs.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/cis_trans_eqtl_MR_analysis_specific_pairs.R
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

return_row_uncor = function(snp, cis_gene, trans_gene, expression_matrix, genotype_matrix, beta_xz, beta_yz, cis_pval, trans_pval) {
	genotype_matrix_part = genotype_matrix[snp,]
	expression_matrix_part = expression_matrix[c(cis_gene, trans_gene),]

    t_MR = calculate_t_MR(as.numeric(genotype_matrix_part[snp,]), as.numeric(expression_matrix_part[cis_gene,]), as.numeric(expression_matrix_part[trans_gene,]), beta_xz, beta_yz)

    return(data.frame(snp = snp, cis_gene = cis_gene, trans_gene = trans_gene, beta_xz = beta_xz, beta_yz = beta_yz, t_MR = t_MR, cis_pval = cis_pval, trans_pval = trans_pval))
}

# Specific examples for paper:
tissue = 'thyroid'
expression_file_location = paste0('/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/', tissue, '_v6p_consortium_autosomes_normalized.txt')
expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
# load covariate matrix
covars = read.csv(paste0('/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/', tissue, '_Analysis.covariates.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

chr = '9'
genotype_file_name = paste0('/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr, '_Final.txt', sep="")
genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
genotype_matrix = genotype_matrix[,colnames(expression_matrix)]

variants = c('rs7037324', 'rs1867277')
# TMEM253, ARFGEF3
trans_genes = c('ENSG00000232070.4', 'ENSG00000112379.8')
# C9org156, FOXE1
cis_genes = c('ENSG00000136932.9', 'ENSG00000178919.7')

snp_pos = data.frame(rsID = variants, chr = 'chr1', pos = 1)
rownames(snp_pos) = variants
gene_pos = data.frame(gene_id = c(cis_genes, trans_genes), chr = c('chr1', 'chr1', 'chr2', 'chr2'), start = 1, end = 2)

# get cis stats
me_cis = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[cis_genes,], snp_pos, gene_pos, pvThresh = 0, pvThresh_cis = 1, covariates = covars)
me_trans = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[trans_genes,], snp_pos, gene_pos, pvThresh = 1, pvThresh_cis = 0, covariates = covars)
me_cis_uncor = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[cis_genes,], snp_pos, gene_pos, pvThresh = 0, pvThresh_cis = 1)
me_trans_uncor = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[trans_genes,], snp_pos, gene_pos, pvThresh = 1, pvThresh_cis = 0)

t_MR_counter = 1
t_MR_list = list()
t_MR_uncor_list = list()


variant = 'rs7037324'
trans_gene = 'ENSG00000232070.4'
cis_gene = 'ENSG00000136932.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs7037324'
trans_gene = 'ENSG00000232070.4'
cis_gene = 'ENSG00000178919.7'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs7037324'
trans_gene = 'ENSG00000112379.8'
cis_gene = 'ENSG00000136932.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs7037324'
trans_gene = 'ENSG00000112379.8'
cis_gene = 'ENSG00000178919.7'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs1867277'
trans_gene = 'ENSG00000232070.4'
cis_gene = 'ENSG00000136932.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs1867277'
trans_gene = 'ENSG00000232070.4'
cis_gene = 'ENSG00000178919.7'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs1867277'
trans_gene = 'ENSG00000112379.8'
cis_gene = 'ENSG00000136932.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs1867277'
trans_gene = 'ENSG00000112379.8'
cis_gene = 'ENSG00000178919.7'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

# t_MR_list[[8]] = return_row(variants[2], cis_genes[2], trans_genes[2], expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[4], me_trans$trans$eqtls$beta[4], covars, me_cis$cis$eqtls$pvalue[4], me_trans$trans$eqtls$pvalue[4])

# variant = variants[2]
# cis_gene = cis_genes[2]
# trans_gene = trans_genes[2]
# beta_xz = me_cis$cis$eqtls$beta[4]
# beta_yz = me_trans$trans$eqtls$beta[4]
# cis_pval = me_cis$cis$eqtls$pvalue[4]
# trans_pval = me_trans$trans$eqtls$pvalue[4]

# t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, beta_xz, beta_yz, covars, cis_pval, trans_pval)
# expression_matrix_permute = expression_matrix[c(1:100),]
# for (k in c(1:100)) {
# 	expression_matrix_permute[k,] = expression_matrix[trans_gene, sample(c(1:ncol(expression_matrix)), replace = FALSE)]
# }
# gene_pos = data.frame(gene_id = paste(trans_gene, c(1:100), sep='_'), chr = 'chr2', start = 1, end = 2)
# rownames(expression_matrix_permute) = gene_pos$gene_id
# me = MatrixEQTL_wrapper(genotype_matrix[variant,], expression_matrix_permute, snp_pos[variant,], gene_pos, pvThresh = 1, pvThresh_cis = 0, covariates = covars)
# t_MRs = lapply(c(1:100), function(x) {return_row(variant, cis_gene, as.character(me$trans$eqtls$gene[x]), rbind(expression_matrix_permute[as.character(me$trans$eqtls$gene[x]),], expression_matrix[cis_gene,]), genotype_matrix, beta_xz, me$trans$eqtls$beta[x], covars, cis_pval, me$trans$eqtls$pvalue[x])})
# t_MR_permute_list[[t_MR_counter]] = do.call('rbind', t_MRs)

# rm(me)

t_MR = do.call('rbind', t_MR_list)
t_MR_uncor = do.call('rbind', t_MR_uncor_list)
t_MR$t_MR = sqrt(t_MR$t_MR)
t_MR_uncor$t_MR = sqrt(t_MR_uncor$t_MR)

save(t_MR, t_MR_uncor, me_cis, me_trans, me_cis_uncor, me_trans_uncor, file = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MR_eQTL/cis_trans_eQTLs/MR_out_thyroid.RData')


# skeletalmuscle example
tissue = 'muscleskeletal'
expression_file_location = paste0('/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/', tissue, '_v6p_consortium_autosomes_normalized.txt')
expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
# load covariate matrix
covars = read.csv(paste0('/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/', tissue, '_Analysis.covariates.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

chr = '5'
genotype_file_name = paste0('/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr, '_Final.txt', sep="")
genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
genotype_matrix = genotype_matrix[,colnames(expression_matrix)]

variants = c('rs2706381', 'rs1012793')
# PSME1, ARTD10
trans_genes = c('ENSG00000092010.10', 'ENSG00000178685.9')
# IRF1
cis_genes = c('ENSG00000125347.9')

snp_pos = data.frame(rsID = variants, chr = 'chr1', pos = 1)
rownames(snp_pos) = variants
gene_pos = data.frame(gene_id = c(cis_genes, trans_genes), chr = c('chr1', 'chr2', 'chr2'), start = 1, end = 2)

# get cis stats
me_cis = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[cis_genes,], snp_pos, gene_pos, pvThresh = 0, pvThresh_cis = 1, covariates = covars)
me_trans = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[trans_genes,], snp_pos, gene_pos, pvThresh = 1, pvThresh_cis = 0, covariates = covars)
me_cis_uncor = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[cis_genes,], snp_pos, gene_pos, pvThresh = 0, pvThresh_cis = 1)
me_trans_uncor = MatrixEQTL_wrapper(genotype_matrix[variants,], expression_matrix[trans_genes,], snp_pos, gene_pos, pvThresh = 1, pvThresh_cis = 0)

t_MR_counter = 1
t_MR_list = list()
t_MR_uncor_list = list()


variant = 'rs2706381'
trans_gene = 'ENSG00000092010.10'
cis_gene = 'ENSG00000125347.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs2706381'
trans_gene = 'ENSG00000178685.9'
cis_gene = 'ENSG00000125347.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs1012793'
trans_gene = 'ENSG00000092010.10'
cis_gene = 'ENSG00000125347.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

variant = 'rs1012793'
trans_gene = 'ENSG00000178685.9'
cis_gene = 'ENSG00000125347.9'
cis_ind = intersect(which(me_cis$cis$eqtls$snps == variant), which(me_cis$cis$eqtls$gene == cis_gene))
trans_ind = intersect(which(me_trans$trans$eqtls$snps == variant), which(me_trans$trans$eqtls$gene == trans_gene))
cis_uncor_ind = intersect(which(me_cis_uncor$cis$eqtls$snps == variant), which(me_cis_uncor$cis$eqtls$gene == cis_gene))
trans_uncor_ind = intersect(which(me_trans_uncor$trans$eqtls$snps == variant), which(me_trans_uncor$trans$eqtls$gene == trans_gene))
t_MR_list[[t_MR_counter]] = return_row(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis$cis$eqtls$beta[cis_ind], me_trans$trans$eqtls$beta[trans_ind], covars, me_cis$cis$eqtls$pvalue[cis_ind], me_trans$trans$eqtls$pvalue[trans_ind])
t_MR_uncor_list[[t_MR_counter]] = return_row_uncor(variant, cis_gene, trans_gene, expression_matrix, genotype_matrix, me_cis_uncor$cis$eqtls$beta[cis_uncor_ind], me_trans_uncor$trans$eqtls$beta[trans_uncor_ind], me_cis_uncor$cis$eqtls$pvalue[cis_uncor_ind], me_trans_uncor$trans$eqtls$pvalue[trans_uncor_ind])
t_MR_counter = t_MR_counter + 1

t_MR = do.call('rbind', t_MR_list)
t_MR_uncor = do.call('rbind', t_MR_uncor_list)
t_MR$t_MR = sqrt(t_MR$t_MR)
t_MR_uncor$t_MR = sqrt(t_MR_uncor$t_MR)

save(t_MR, t_MR_uncor, me_cis, me_trans, me_cis_uncor, me_trans_uncor, file = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MR_eQTL/cis_trans_eQTLs/MR_out_muscleskeletal.RData')
