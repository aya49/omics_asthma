##################################################
#  wabf_calculate_bayes_factor_MR_stats.R
#
#  $proj/Scripts/causality/bayes_MR/wabf_calculate_bayes_factor_MR_stats.R
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
args = c(1:7)
args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Muscle_Skeletal.v8.normalized_expression.bed.gz'
args[2] = '10'
args[3] = '1'
args[4] = '20000'
args[5] = 'Muscle_Skeletal'
args[6] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
args[7] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/raw/Muscle_Skeletal/bayes_freq_MR_stats'

expression_file_location = args[1]
chr_number = args[2]
part_number = as.numeric(args[3])
partition_size = as.numeric(args[4])
tissue_name = args[5]
cov_dir = args[6]
out_file = args[7]

# Read in the expression files and the gene positions
header = readLines(gzfile(expression_file_location), n = 1)
header = strsplit(header, '\t')[[1]]
expression_matrix = read.csv(gzfile(expression_file_location, 'r'), sep = '\t', stringsAsFactors = FALSE)

colnames(expression_matrix) = header
rownames(expression_matrix) = expression_matrix$gene_id

# Filter out genes with low mappability
mappability_list = read.table('/tigress/BEE/RNAseq/Output/processing/mappability/annotation/hg38_gene_mappability.txt', col.names = c('gene', 'mappability'), stringsAsFactors = FALSE)
rownames(mappability_list) = mappability_list$gene
expression_matrix = expression_matrix[sapply(rownames(expression_matrix), function(x) {(x %in% mappability_list$gene) && (mappability_list[x,2] >= 0.8)}),]

# Only take the genes that are in the filtered list:
filtered_gene_list = read.table('/tigress/BEE/RNAseq/Data/Resources/annotations/silver/gene_filter/Bayes_MR_list.txt', stringsAsFactors = FALSE)
expression_matrix = expression_matrix[sapply(rownames(expression_matrix), function(x) {x %in% filtered_gene_list$V1}),]

gene_positions = expression_matrix[,c(4,1,2,3)]
colnames(gene_positions) = c('gene_id', 'chr', 'start', 'end')
expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]

# Read in the genotype positions - both for MAF 1% and 5%
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_05_ld_pruned.RData')
load(genotype_file_name)
genotypes_all = genotype_matrix_master
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_01_ld_pruned.RData')
load(genotype_file_name)
genotypes_all = rbind(genotypes_all, genotype_matrix_master)
# This loads in the data frame "genotype_matrix_master"

# Get the appropriate partition
num_parts = ceiling(nrow(genotypes_all) / partition_size)
num_inds = ceiling(nrow(genotypes_all) / num_parts)
if (part_number == num_parts) {
  genotype_matrix = genotypes_all[c((((part_number-1)*num_inds)+1):nrow(genotypes_all)),]
} else {
  genotype_matrix = genotypes_all[c((((part_number-1)*num_inds)+1):(part_number*num_inds)),]
}
# There should be 838 indivs for genotype

# Make sure the columns are the same
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

# Fix the data type
genotype_matrix_temp = data.frame(lapply(genotype_matrix,as.numeric))
rownames(genotype_matrix_temp) = rownames(genotype_matrix)
colnames(genotype_matrix_temp) = colnames(genotype_matrix)
genotype_matrix = genotype_matrix_temp
# Fix the genotype matrix colnames

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

# Hack for ASHG - only use v6p samples - remove later!
# v6p_subject_list = read.table('/tigress/BEE/RNAseq/Data/Resources/gtex/genotype/subjects_with_genotypes.txt', stringsAsFactors=F)
# v6p_subject_list = as.character(v6p_subject_list$V1)

# present = sapply(v6p_subject_list, function(x) {x %in% colnames(genotype_matrix)})

# genotype_matrix = genotype_matrix[,v6p_subject_list[present]]
# expression_matrix = expression_matrix[,v6p_subject_list[present]]
# cov = cov[v6p_subject_list[present],]
# cov = center_colmeans(cov)

# remove covs columns that are not unique
if (length(unique(cov$pcr)) == 1) {cov = subset(cov, select = -c(ncol(cov)-2))}
if (length(unique(cov$platform)) == 1) {cov = subset(cov, select = -c(ncol(cov)-1))}
if (length(unique(cov$sex)) == 1) {cov = subset(cov, select = -c(ncol(cov)))}

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

# Simpler version only for betas
calc_betas = function(exp, temp_cov) {
	# calculate the univariate WABF
	X = as.matrix(temp_cov)
	y = as.matrix(exp)
	Z = solve(t(X) %*% X) %*% t(X) %*% y
	# strength of association
	betas = Z[nrow(Z),]
	return(betas)
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
calc_MR_ABF = function(exp_cis, exp_trans, temp_cov, gene_list, beta_xz, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1) {
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
	# scale the trans beta with cis beta for beta IV
	betas[,1] = betas[,1]/beta_xz
	# empirical variance V
	total_V = var(as.matrix(betas))
	total_V[1,2] = total_V[2,1] = 0
	# take out the empirical mean? currently just take the zero vector
	# total_means = colMeans(as.matrix(betas))
	total_means = as.matrix(c(0,0))

	H00_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V) %*% t(as.matrix(betas[x,] - total_means))))})
	total_V_01 = total_V
	total_V_01[2,2] = total_V_01[2,2] + W_MR_2
	H01_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V_01))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V_01) %*% t(as.matrix(betas[x,] - total_means))))})
	total_V_10 = total_V
	total_V_10[1,1] = total_V_10[1,1] + W_MR_1
	H10_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V_10))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V_10) %*% t(as.matrix(betas[x,] - total_means))))})
	total_V_11 = total_V
	total_V_11[1,1] = total_V_11[1,1] + W_MR_1
	total_V_11[2,2] = total_V_11[2,2] + W_MR_2
	H11_ABF = sapply(gene_list, function(x) {(1 / sqrt(det(total_V_11))) * exp(-0.5 * (as.matrix(betas[x,] - total_means) %*% solve(total_V_11) %*% t(as.matrix(betas[x,] - total_means))))})

	return_frame = data.frame(H00_ABF, H01_ABF, H10_ABF, H11_ABF)
	return_frame$MR_PPA = ((snp_pi_1 * H10_ABF) + (snp_pi_1 * gene_pi_1 * H11_ABF)) / (H00_ABF + (snp_pi_1 * H10_ABF) + (gene_pi_1 * H01_ABF) + (snp_pi_1 * gene_pi_1 * H11_ABF))
	return_frame = cbind(betas[gene_list,], return_frame)
	return(return_frame)
}

# For now, let use the W value of 0.15^2 - roughly translating to 95% chance that the relative risk is between 2/3 and 3/2
cis_risk = 1.5
W = (log(cis_risk)/1.96)^2
trans_risk_1 = 1.5
W_MR_1 = (log(trans_risk_1)/1.96)^2
trans_risk_2 = 1.1
W_MR_2 = (log(trans_risk_2)/1.96)^2
# Say we expect roughly one out of 1e5 SNPs to be eQTLs
pi_1 = 1e-5
PO = (1-pi_1)/pi_1
# parameters for MR-ABF
snp_pi_1 = 1e-3 # expecting 1 out of ~10000 true trans-eQTLs among chosen cis-eQTLs
gene_pi_1 = 1e-3 # expecting nonzero contribution in ~100 trans genes
# trans distance threshold - will set 1Mb for both
cis_threshold = 1000000
trans_threshold = 1000000

# Parts of the pipeline that are subject to change
# Separate the Samples into AA and EA individuals
pop_clust_dir = '/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/support_files/pop_clust/'
AA_list_1 = read.table(paste0(pop_clust_dir, 'AA_1.txt'), stringsAsFactors = FALSE)
AA_list_2 = read.table(paste0(pop_clust_dir, 'AA_2.txt'), stringsAsFactors = FALSE)
EA_list_1 = read.table(paste0(pop_clust_dir, 'EA_1.txt'), stringsAsFactors = FALSE)
EA_list_2 = read.table(paste0(pop_clust_dir, 'EA_2.txt'), stringsAsFactors = FALSE)
EA_list_3 = read.table(paste0(pop_clust_dir, 'EA_3.txt'), stringsAsFactors = FALSE)
AA_list = c(AA_list_1$V1, AA_list_2$V1)
EA_list = c(EA_list_1$V1, EA_list_2$V1, EA_list_3$V1)
# Which ones are available?
AA_list = AA_list[sapply(AA_list, function(x) {x %in% colnames(expression_matrix)})]
EA_list = EA_list[sapply(EA_list, function(x) {x %in% colnames(expression_matrix)})]

abf_eqtl_list = list()
MR_stats_list = list()

# Let's first obtain the list of eQTLs - PPA threshold of 0.5 for cis and 0.05 for trans
# This will be divided into EA and AA individuals and processed separately
count = 1
for (i in c(1:nrow(genotype_matrix))) {
	snp = rownames(genotype_matrix)[i]


	# Which genotypes are not NA?
	inds = !is.na(genotype_matrix[snp,])
	if (length(unique(as.numeric(genotype_matrix[snp,inds]))) == 1) {next}
	# Get the SNP MAF
	snp_maf = sum(as.numeric(genotype_matrix[snp,inds]))/(2*sum(inds))
	snp_maf = min(snp_maf, 1-snp_maf)
	# Get the expression subset corresponding to available genotypes
	exp_all_genes = t(expression_matrix)[as.logical(inds),]
	exp_all_genes = center_colmeans(exp_all_genes)
	temp_cov = cov[inds,]
	temp_cov$SNP = as.numeric(genotype_matrix[i,inds])
	# Mean center
	temp_cov$SNP = temp_cov$SNP - mean(temp_cov$SNP)
	# Calculate the general WABF
	eqtl_candidates = calc_WABF(exp_all_genes, temp_cov, W, PO, 0.5)
  	# save the stats for any associations with PPA over 0.5
  	if (nrow(eqtl_candidates) == 0) {next}
    # Record the eQTL candidates
    eqtl_candidates = cbind(snp, rownames(eqtl_candidates), eqtl_candidates)
    colnames(eqtl_candidates)[2] = 'gene'
    eqtl_candidates$gene_chr = gene_positions[as.character(eqtl_candidates$gene),]$chr
	eqtl_candidates$gene_start = gene_positions[as.character(eqtl_candidates$gene),]$start
	eqtl_candidates$gene_end = gene_positions[as.character(eqtl_candidates$gene),]$end
	eqtl_candidates$snp_maf = snp_maf
    abf_eqtl_list[[i]] = eqtl_candidates

    # Mark the gene locations
	cis_eqtl_list = abf_eqtl_list[[i]]
	cis_eqtl_list = cis_eqtl_list[cis_eqtl_list$gene_chr == paste0('chr', chr_number),]

	diff = abs(as.numeric(sapply(as.character(cis_eqtl_list$snp), function(x) {strsplit(x, '_')[[1]][2]})) - cis_eqtl_list$gene_start)
	# Say our cis-threshold is 150 kb
	cis_eqtl_list = cis_eqtl_list[diff <= cis_threshold,]

	if (nrow(cis_eqtl_list) == 0) {next}
	print('cis_eqtl found')
	print(i)
	print(nrow(cis_eqtl_list))
	# If there are cis eqtls to test for:
	for (j in c(1:nrow(cis_eqtl_list))) {
		print(count)
		cis_gene = as.character(cis_eqtl_list$gene[j])
		# which genes to include for trans analysis?
		snp_chr = strsplit(snp,'_')[[1]][1]
		snp_pos = as.numeric(strsplit(snp,'_')[[1]][2])
		# Get the set of trans genes to test for
		gene_inds = sapply(c(1:nrow(gene_positions)), function(x) {!(gene_positions[x,'chr'] == snp_chr && abs(gene_positions[x,'start'] - snp_pos) <= trans_threshold)})
		exp = exp_all_genes[,gene_inds]
		# exp = center_colmeans(exp)
		# Permuted version for comparison
		set.seed(111)
		exp_perm = exp[sample(nrow(exp)),]

		# The set of trios to test for - get all trios that have trans-eQTL PPA >= 5 percent
		trans_candidates = calc_WABF(exp, temp_cov, W, PO, 0.05)
		if (nrow(trans_candidates) == 0) {next}

		# Calculate frequentist MR Wald stats
		exp_cis = exp_all_genes[,cis_gene]
		exp_trans_cands = exp[,rownames(trans_candidates)]
		exp_trans_cands_perm = exp_perm[,rownames(trans_candidates)]
		# Get the trans-betas for MR stats
		trans_candidates_perm = calc_betas(exp_trans_cands_perm, temp_cov)
		trans_candidates$betas_perm = as.numeric(trans_candidates_perm)
		beta_xz = cis_eqtl_list$betas[j]
		trans_candidates_MR_stats = calc_freq_MR(exp_trans_cands, exp_trans_cands_perm, temp_cov, exp_cis, beta_xz, trans_candidates)

		# Calculate Bayesian MR-ABF
		gene_list = rownames(trans_candidates)
		MR_ABF_df = calc_MR_ABF(exp_cis, exp, temp_cov, gene_list, beta_xz, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)
		MR_ABF_df_perm = calc_MR_ABF(exp_cis, exp_perm, temp_cov, gene_list, beta_xz, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)

		total_df = cbind(trans_candidates_MR_stats, MR_ABF_df)
		colnames(MR_ABF_df_perm) = c('beta_trans_perm', 'theta_perm', 'H00_ABF_perm', 'H01_ABF_perm', 'H10_ABF_perm', 'H11_ABF_perm', 'MR_PPA_perm')
		total_df = cbind(total_df, MR_ABF_df_perm)
		# Add in additional information
		total_df = cbind(gene_positions[rownames(total_df),], total_df)
		total_df = cbind(gene_positions[cis_gene,], total_df)
		total_df = cbind(snp, total_df)
		colnames(total_df)[c(1:9)] = c('snp', 'cis_gene', 'cis_chr', 'cis_start', 'cis_end', 'trans_gene', 'trans_chr', 'trans_start', 'trans_end')
		total_df$snp_maf = snp_maf

		MR_stats_list[[count]] = total_df
		count = count + 1
	}
}

eqtl_out_df = do.call('rbind', abf_eqtl_list)
mr_out_df = do.call('rbind', MR_stats_list)

# Save the resulting data frames with the parameter set as file name
save(eqtl_out_df, mr_out_df, file = paste0(out_file, '_chr_', chr_number, '_part_', part_number, '_cis_risk_', cis_risk, '_pi1_', pi_1,  '_trans_risk1_', trans_risk_1, '_trans_risk2_', trans_risk_2, '_snp_effect_', snp_pi_1, '_gene_effect_', gene_pi_1, '.RData'))
