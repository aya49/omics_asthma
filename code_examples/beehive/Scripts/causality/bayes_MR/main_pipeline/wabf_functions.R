##################################################
#  wabf_functions.R
#
#  $proj/Scripts/causality/bayes_MR/main_pipeline/wabf_functions.R
# 
#  This helper functions for the bayes MR pipeline.
#
#  Author: Brian Jo
#
##################################################

# Helper functions
calc_WABF_with_NAs = function(exp, genotype, cov, W, PO) {
	# Remove rows from genotypes that are NA
	inds = !is.na(genotype)
	# Return if genotype is uniform
	if (length(unique(as.numeric(genotype[,inds]))) == 1) {return(NULL)}
	# Get the SNP MAF
	snp_maf = sum(as.numeric(genotype[,inds]))/(2*sum(inds))
	snp_maf = min(snp_maf, 1-snp_maf)
	# Get the expression subset corresponding to available genotypes
	exp_all_genes = t(exp)[as.logical(inds),]
	exp_all_genes = center_colmeans(exp_all_genes)
	vars = cov[inds,]
	vars = center_colmeans(vars)
	# What is the proportion of variance explained by covs?
	y = as.matrix(genotype[,inds] - mean(as.numeric(genotype[,inds])) )
	x = as.matrix(vars)
	beta_x = tryCatch({(x %*% solve(t(x) %*% x) %*% t(x) %*% t(y))}, error = function(e) {c(NA)})
	if (is.na(beta_x[1])) {return(NULL)}
	var_exp = sum(beta_x * beta_x)/sum(y * y)
	if (var_exp > 0.9) {return(NULL)}
	# Now, append the SNP to the variable set
	vars$SNP = as.numeric(y)
	# Set up the linear equation
	X = as.matrix(vars)
	y = as.matrix(exp_all_genes)
	fit = solve(t(X) %*% X) %*% t(X) %*% y
	betas = fit[nrow(fit),]
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

	return_frame = data.frame(row.names = NULL, snp = rownames(genotype), gene = colnames(exp_all_genes), beta = betas, stat = Z, abf = ABF, ppa = PPA, maf = snp_maf, snp_var_exp = var_exp)
	return(return_frame)
}

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
calc_betas = function(exp, genotype, cov) {
	# calculate the univariate WABF
	vars = cbind(cov, as.numeric(genotype))
	X = as.matrix(vars)
	y = as.matrix(exp)
	Z = solve(t(X) %*% X) %*% t(X) %*% y
	# strength of association
	betas = Z[nrow(Z),]
	return(as.numeric(betas))
}

calc_freq_MR_with_NAs = function(exp_cis, exp_trans, exp_trans_perm, genotype, cov, beta_xz = NULL, beta_yz = NULL, beta_yz_perm = NULL) {
	# Remove rows from genotypes that are NA
	inds = !is.na(genotype)
	# Return if genotype is uniform
	if (length(unique(as.numeric(genotype[,inds]))) == 1) {return(NULL)}

	# Format the expression values
	exp_cis = t(exp_cis[,inds])
	exp_cis = center_colmeans(exp_cis)
	exp_trans = t(exp_trans[,inds])
	exp_trans = center_colmeans(exp_trans)
	exp_trans_perm = t(exp_trans_perm[,inds])
	exp_trans_perm = center_colmeans(exp_trans_perm)

	# Get the beta values if necessary
	if (is.null(beta_xz)) {beta_xz = calc_betas(exp_cis, genotype[,inds], cov[inds,])}
	if (is.null(beta_yz)) {beta_yz = calc_betas(exp_trans, genotype[,inds], cov[inds,])}
	if (is.null(beta_yz_perm)) {beta_yz_perm = calc_betas(exp_trans_perm, genotype[,inds], cov[inds,])}

	# Set up the linear equation
	X = as.matrix(cov[inds,])
	X = center_colmeans(X)
	g = as.matrix(t(genotype[,inds]))
	g = center_colmeans(g)

	fit_space = tryCatch({X %*% solve(t(X) %*% X) %*% t(X)}, error = function(e) {c(NA)})
	if (is.na(fit_space[1])) {return(NULL)}

	# orthogonalize w.r.t. covariates
	res_y = exp_trans - (fit_space %*% exp_trans)
	res_y_perm = exp_trans_perm - (fit_space %*% exp_trans_perm)
	res_x = exp_cis - (fit_space %*% exp_cis)
	res_g = g - (fit_space %*% g)
	# Calculate the Wald stats
	beta_mr = beta_yz/beta_xz
	# res_x is the cis-gene expression value - essentially (y-beta*x)^T * (y-beta*x) / (n-3)
	sig_sq = sapply(c(1:length(beta_mr)), function(n) {t(res_y[,n] - res_x*beta_mr[n]) %*% (res_y[,n] - res_x*beta_mr[n]) / (length(res_g) - 3)})
	var_beta = sig_sq * as.numeric((t(res_g) %*% res_g) / (t(res_x) %*% res_g)^2)
	MR_Wald = (beta_mr)^2/var_beta
	# Calculate the Wald stats for permuted
	beta_mr_perm = beta_yz_perm/beta_xz
	sig_sq_perm = sapply(c(1:length(beta_mr_perm)), function(n) {t(res_y_perm[,n] - res_x*beta_mr_perm[n]) %*% (res_y_perm[,n] - res_x*beta_mr_perm[n]) / (length(res_g) - 3)})
	var_beta_perm = sig_sq_perm * as.numeric((t(res_g) %*% res_g) / (t(res_x) %*% res_g)^2)
	MR_Wald_perm = (beta_mr_perm)^2/var_beta_perm

	return(data.frame(MR_Wald = MR_Wald, MR_Wald_perm = MR_Wald_perm))
}

calc_freq_MR = function(exp_trans, exp_trans_perm, genotype, cov, exp_cis, beta_xz, beya_yz, trans_candidates) {
	# calculate the frequentist MR Wald stats
	X = as.matrix(temp_cov[,c(1:ncol(temp_cov)-1)])
	g = as.matrix(temp_cov[,ncol(temp_cov)])
	# orthogonalize w.r.t. covariates
	res_y = exp_trans - (X %*% solve(t(X) %*% X) %*% t(X) %*% exp_trans)
	res_y_perm = exp_trans_perm - (X %*% solve(t(X) %*% X) %*% t(X) %*% exp_trans_perm)
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

calc_MR_ABF_with_NAs = function(exp_cis, exp_trans, genotype, cov, beta_xz, gene_list, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1) {
	# Remove rows from genotypes that are NA
	inds = !is.na(genotype)
	# Return if genotype is uniform
	if (length(unique(as.numeric(genotype[,inds]))) == 1) {return(NULL)}

	# Format the expression values
	exp_cis = t(exp_cis[,inds])
	exp_cis = center_colmeans(exp_cis)
	exp_trans = t(exp_trans[,inds])
	exp_trans = center_colmeans(exp_trans)

	# Calculate the Bayesian MR-ABF
	X = as.matrix(cov[inds,])
	X = center_colmeans(X)
	g = as.matrix(t(genotype[,inds]))
	g = center_colmeans(g)
	# orthogonalize X with respect to genotype
	ortho_covs = X - (g %*% solve(t(g) %*% g) %*% t(g) %*% X)
	ortho_covs = cbind(ortho_covs, g)
	# orthogonalize exp_cis with respect to all covariates
	exp_cis_ortho = tryCatch({exp_cis - (ortho_covs %*% solve(t(ortho_covs) %*% ortho_covs) %*% t(ortho_covs) %*% exp_cis)}, error = function(e) {c(NA)})
	if (is.na(exp_cis_ortho[1])) {return(NULL)}
	ortho_covs = cbind(ortho_covs, exp_cis_ortho)
	# Add the appropriate column names
	# colnames(ortho_covs) = c(colnames(ortho_covs)[c(1:(ncol(ortho_covs)-2))], c('SNP', 'exp_cis'))
	# Now solve fit exp_trans jointly with respect to all covariates, genotype, and exp_cis
	Z = tryCatch({solve(t(ortho_covs) %*% ortho_covs) %*% t(ortho_covs) %*% exp_trans}, error = function(e) {c(NA)})
	if (is.na(Z[1])) {return(NULL)}
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

# center with 'colMeans()'
center_colmeans = function(x) {
    xcenter = colMeans(x)
    x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}