##################################################
#  calculate_bayes_freq_MR_stats.R
#
#  $proj/Scripts/causality/bayes_MR/calculate_bayes_freq_MR_stats.R
# 
#  This script calculates the MR stats (both frequentist and Bayesian) for SNP-cis-trans trios for the GTEx v8 data.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

# The first part of the script is identical to the GTEx v8 frequentist eQTL mapping pipeline

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

library(MASS)
# Example
# args = c(1:6)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz'
# args[2] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/WABF/raw/Whole_Blood/wabf_raw_output_chr10_part1_risk_1.5_pi1_1e-05.RData'
# args[3] = '10'
# args[4] = 'Whole_Blood'
# args[5] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
# args[6] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/raw/Whole_Blood/bayes_freq_MR_stats_chr10_part1_risk_1.5_pi1_1e-05_'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
WABF_raw_output_location = args[2]
chr_number = args[3]
tissue_name = args[4]
cov_dir = args[5]
out_file = args[6]

# Read in the list of cis-eQTLs
load(WABF_raw_output_location)

# Read in the expression files and the gene positions
header = readLines(gzfile(expression_file_location), n = 1)
header = strsplit(header, '\t')[[1]]
expression_matrix = read.csv(gzfile(expression_file_location, 'r'), sep = '\t', stringsAsFactors = FALSE)

colnames(expression_matrix) = header
rownames(expression_matrix) = expression_matrix$gene_id

gene_positions = expression_matrix[,c(4,1,2,3)]
colnames(gene_positions) = c('gene_id', 'chr', 'start', 'end')
expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]

# Read in the genotype positions
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_05_ld_pruned.RData')
load(genotype_file_name)
# This loads in the data frame "genotype_matrix_master"

# Make sure the columns are the same
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix_master))]
genotype_matrix = genotype_matrix_master[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix_master))})]

# Fix the data type
genotype_matrix_temp = data.frame(lapply(genotype_matrix,as.numeric))
rownames(genotype_matrix_temp) = rownames(genotype_matrix)
genotype_matrix = genotype_matrix_temp

# Get the SNP positions
snp_positions = data.frame(ID = rownames(genotype_matrix), chr = sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][1]}), pos = as.numeric(sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][2]})))

# Load in the covariates
suffix = '.v8.covariates.txt'
covars = read.csv(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
colnames(covars) = as.character(sapply(colnames(covars), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

# invert and mean-center covars
inv_cov = t(covars)
# center with 'colMeans()'
center_colmeans = function(x) {
    xcenter = colMeans(x)
    x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}
# apply it
inv_cov = center_colmeans(inv_cov)
cov = data.frame(inv_cov)

# Mark the gene locations
wabf_df = out_df

wabf_df$gene_chr = gene_positions[as.character(out_df$gene),]$chr
wabf_df$gene_start = gene_positions[as.character(out_df$gene),]$start
wabf_df$gene_end = gene_positions[as.character(out_df$gene),]$end
wabf_df = wabf_df[wabf_df$gene_chr == paste0('chr', chr_number),]

diff = abs(as.numeric(sapply(as.character(wabf_df$snp), function(x) {strsplit(x, '_')[[1]][2]})) - wabf_df$gene_start)
# Say our cis-threshold is 150 kb
threshold = 150000
cis_eqtl_list = wabf_df[diff <= threshold,]

# For each row, do the following:
# Get the list of genes to test for trans-effects - say the threshold is 1 mb
# calculate the betas for SNP-trans gene
# calculate frequentist MR Wald stats
# calculate the MR-ABF (Baye factors)
# Output

# Helper functions
calc_WABF = function(exp, temp_cov, W, PO, thresh) {
	# calculate the univariate WABF
	X = as.matrix(temp_cov)
	y = as.matrix(exp)
	Z = solve(t(X) %*% X) %*% t(X) %*% y
	# strength of association
	betas = Z[nrow(Z),]
	# empirical variance V
	V = var(betas)
	# shrinkage factor r
	r = W/(V+W)
	# Z-statistic
	Z = betas / sqrt(V)
	# Approximate Bayes Factor
	ABF = (1/sqrt(1-r)) * exp(-r*(Z^2)/2)
	# Posterior probability of association
	PPA = 1 / (ABF*PO + 1)

	return_frame = data.frame(betas, Z, ABF, PPA)
	# Let's relax the criteria a bit to 0.1 - corresponding to roughly ABF 9e-5
	return_frame = return_frame[return_frame$PPA >= thresh,]
	return(return_frame)
}

calc_freq_MR = function(exp_trans_cands, exp_trans_cands_perm, temp_cov, exp_cis, beta_xz, trans_candidates) {
	# calculate the frequentist MR Wald stats
	X = as.matrix(temp_cov[,c(1:ncol(temp_cov)-1)])
	g = as.matrix(temp_cov[,ncol(temp_cov)])
	# orthogonalize w.r.t. covariates
	res_y = exp_trans_cands - (X %*% solve(t(X) %*% X) %*% t(X) %*% exp_trans_cands)
	res_y_perm = exp_trans_cands_perm - (X %*% solve(t(X) %*% X) %*% t(X) %*% exp_trans_cands_perm)
	res_y_cis = exp_cis - (X %*% solve(t(X) %*% X) %*% t(X) %*% exp_cis)
	res_g = g - (X %*% solve(t(X) %*% X) %*% t(X) %*% g)
	# Calculate the Wald stats
	beta_mr = trans_candidates$betas/beta_xz
	# res_y_cis is the cis-gene expression value - essentially (y-beta*x)^T * (y-beta*x) / (n-3)
	sig_sq = sapply(c(1:length(beta_mr)), function(n) {t(res_y[,n] - res_y_cis*beta_mr[n]) %*% (res_y[,n] - res_y_cis*beta_mr[n]) / (length(res_g) - 3)})
	var_beta = sig_sq * as.numeric((t(res_g) %*% res_g) / (t(res_y_cis) %*% res_g)^2)
	trans_candidates_MR_stats = trans_candidates
	trans_candidates_MR_stats$MR_Wald = (beta_mr)^2/var_beta
	# Calculate the Wald stats for permuted
	beta_mr = trans_candidates$betas_perm/beta_xz
	sig_sq = sapply(c(1:length(beta_mr)), function(n) {t(res_y_perm[,n] - res_y_cis*beta_mr[n]) %*% (res_y_perm[,n] - res_y_cis*beta_mr[n]) / (length(res_g) - 3)})
	var_beta = sig_sq * as.numeric((t(res_g) %*% res_g) / (t(res_y_cis) %*% res_g)^2)
	trans_candidates_MR_stats$MR_Wald_perm = (beta_mr)^2/var_beta

	return(trans_candidates_MR_stats)
}

# TODO: also report the betas
calc_MR_ABF = function(exp_cis, exp_trans, temp_cov, gene_list, W_MR, snp_pi_1, gene_pi_1) {
	# Calculate the Bayesian MR-ABF
	# orthogonalize X with respect to genotype
	X = as.matrix(temp_cov)
	# X[,ncol(X)] is the genotype
	ortho_covs = X[,c(1:(ncol(X)-1))] - (X[,ncol(X)] %*% solve(t(X[,ncol(X)]) %*% X[,ncol(X)]) %*% t(X[,ncol(X)]) %*% X[,c(1:(ncol(X)-1))])
	ortho_covs = cbind(ortho_covs, X[,'SNP'])
	# orthogonalize exp_cis with respect to all covariates
	exp_cis_ortho = exp_cis - (ortho_covs %*% solve(t(ortho_covs) %*% ortho_covs) %*% t(ortho_covs) %*% exp_cis)
	ortho_covs = cbind(ortho_covs, exp_cis_ortho)
	# Add the appropriate column names
	colnames(ortho_covs) = c(colnames(ortho_covs)[c(1:(ncol(ortho_covs)-2))], c('SNP', 'exp_cis'))
	# Now solve fit exp_trans jointly with respect to all covariates, genotype, and exp_cis
	Z = solve(t(ortho_covs) %*% ortho_covs) %*% t(ortho_covs) %*% exp_trans
	# strength of association
	betas = data.frame(beta_trans = Z[(nrow(Z)-1),], theta = Z[nrow(Z),])
	# empirical variance V
	total_V = var(as.matrix(betas))
	# take out the empirical mean? currently just take the zero vector
	# total_means = colMeans(as.matrix(betas))
	total_means = as.matrix(c(0,0))

	H00_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V) %*% t(as.matrix(betas[x,] - total_means))))})
	total_V_01 = total_V
	total_V_01[1,1] = total_V_01[1,1] + W_MR
	H01_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V_01))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V_01) %*% t(as.matrix(betas[x,] - total_means))))})
	total_V_10 = total_V
	total_V_10[2,2] = total_V_10[2,2] + W_MR
	H10_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V_10))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V_10) %*% t(as.matrix(betas[x,] - total_means))))})
	total_V_11 = total_V
	total_V_11[1,1] = total_V_11[1,1] + W_MR
	total_V_11[2,2] = total_V_11[2,2] + W_MR
	H11_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V_11))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V_11) %*% t(as.matrix(betas[x,] - total_means))))})

	return_frame = data.frame(H00_ABF, H01_ABF, H10_ABF, H11_ABF)
	return_frame$MR_PPA = ((gene_pi_1 * H10_ABF) + (snp_pi_1 * gene_pi_1 * H11_ABF)) / (H00_ABF + (snp_pi_1 * H01_ABF) + (gene_pi_1 * H10_ABF) + (snp_pi_1 * gene_pi_1 * H11_ABF))
	return(return_frame)
}

# For now, let use the W value of 0.15^2 - roughly translating to 95% chance that the relative risk is between 2/3 and 3/2
cis_risk = 1.5
W = (log(cis_risk)/1.96)^2
trans_risk = 1.5
W_MR = (log(trans_risk)/1.96)^2
# Say we expect roughly one out of 1e5 SNPs to be eQTLs
pi_1 = 1e-5
PO = (1-pi_1)/pi_1
# parameters for MR-ABF
snp_pi_1 = 1e-4 # expecting 1 out of ~10000 true trans-eQTLs among chosen cis-eQTLs
gene_pi_1 = 5e-3 # expecting nonzero contribution in ~100 trans genes
# trans distance threshold
threshold = 1000000

MR_stats_list = list()

for (i in c(1:nrow(cis_eqtl_list))) {
	print(i)
	snp = as.character(cis_eqtl_list$snp[i])
	cis_gene = as.character(cis_eqtl_list$gene[i])
	# Which genotypes are not NA?
	inds = !is.na(genotype_matrix[snp,])
	if (length(unique(as.numeric(genotype_matrix[snp,inds]))) == 1) {next}
	# Get the expression subset corresponding to available genotypes
	exp = t(expression_matrix)[as.logical(inds),]
	# which genes to include for trans analysis?
	snp_chr = strsplit(snp,'_')[[1]][1]
	snp_pos = as.numeric(strsplit(snp,'_')[[1]][2])
	# Get the set of trans genes to test for
	gene_inds = sapply(c(1:nrow(gene_positions)), function(x) {!(gene_positions[x,'chr'] == snp_chr && abs(gene_positions[x,'start'] - snp_pos) <= threshold)})
	exp = exp[,gene_inds]
	exp = center_colmeans(exp)
	# Permuted version for comparison
	set.seed(111)
	exp_perm = exp[sample(nrow(exp)),]

	# Calculate the trans betas and ABF
	temp_cov = cov[inds,]
	temp_cov$SNP = as.numeric(genotype_matrix[snp,inds])
	# Mean center
	temp_cov$SNP = temp_cov$SNP - mean(temp_cov$SNP)

	# The set of trios to test for - get all trios that have trans-eQTL PPA >= 5 percent
	trans_candidates = calc_WABF(exp, temp_cov, W, PO, 0.05)
	if (nrow(trans_candidates) == 0) {next}

	# Calculate frequentist MR Wald stats
	exp_cis = t(expression_matrix)[as.logical(inds),cis_gene]
	exp_trans_cands = exp[,rownames(trans_candidates)]
	exp_trans_cands_perm = exp_perm[,rownames(trans_candidates)]
	# Get the trans-betas for MR stats
	trans_candidates_perm = calc_WABF(exp_trans_cands_perm, temp_cov, W, PO, 0)
	trans_candidates$betas_perm = sapply(rownames(trans_candidates), function(x) {trans_candidates_perm[x,'betas']})
	beta_xz = cis_eqtl_list[i,]$beta
	trans_candidates_MR_stats = calc_freq_MR(exp_trans_cands, exp_trans_cands_perm, temp_cov, exp_cis, beta_xz, trans_candidates)

	# Calculate Bayesian MR-ABF
	MR_ABF_df = calc_MR_ABF(exp_cis, exp, temp_cov, rownames(trans_candidates), W_MR, snp_pi_1, gene_pi_1)
	MR_ABF_df_perm = calc_MR_ABF(exp_cis, exp_perm, temp_cov, rownames(trans_candidates), W_MR, snp_pi_1, gene_pi_1)

	total_df = cbind(trans_candidates_MR_stats, MR_ABF_df)
	colnames(MR_ABF_df_perm) = c('H00_ABF_perm', 'H01_ABF_perm', 'H10_ABF_perm', 'H11_ABF_perm', 'MR_PPA_perm')
	total_df = cbind(total_df, MR_ABF_df_perm)
	# Add in additional information
	total_df = cbind(gene_positions[rownames(total_df),], total_df)
	total_df = cbind(gene_positions[cis_gene,], total_df)
	total_df = cbind(snp, total_df)
	colnames(total_df)[c(1:9)] = c('snp', 'cis_gene', 'cis_chr', 'cis_start', 'cis_end', 'trans_gene', 'trans_chr', 'trans_start', 'trans_end')

	MR_stats_list[[i]] = total_df
}

out_df = do.call('rbind', MR_stats_list)

save(out_df, file = paste0(out_file, '_trans_risk_', trans_risk, '_snp_effect_', snp_pi_1, '_gene_effect_', gene_pi_1, '.RData'))
