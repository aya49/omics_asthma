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

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# Example
args = c(1:7)
args[1] = 'Muscle_Skeletal'
args[2] = '10'
args[3] = '1'
args[4] = '20000'
args[5] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/prep_files/'
args[6] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/raw/Muscle_Skeletal/bayes_freq_MR_stats'

tissue_name = args[1]
chr_number = args[2]
part_number = as.numeric(args[3])
partition_size = as.numeric(args[4])
prep_dir = args[5]
out_file = args[6]

library(MASS)
source("/tigress/BEE/RNAseq/Scripts/causality/bayes_MR/wabf_functions.R")
load(paste0(prep_dir, tissue_name, '.RData'))

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

# Fix the data type
genotype_matrix_temp = data.frame(lapply(genotype_matrix,as.numeric))
rownames(genotype_matrix_temp) = rownames(genotype_matrix)
colnames(genotype_matrix_temp) = colnames(genotype_matrix)
genotype_matrix = genotype_matrix_temp

# match the colnames to expression matrix
genotype_matrix = genotype_matrix[,colnames(expression_matrix)]
genotype_matrix_AA = genotype_matrix[,colnames(expression_matrix_AA)]
genotype_matrix_EA = genotype_matrix[,colnames(expression_matrix_EA)]

# Get the SNP positions
snp_positions = data.frame(ID = rownames(genotype_matrix), chr = sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][1]}), pos = as.numeric(sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][2]})))

# Parameter settings
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

abf_eqtl_list_AA = list()
abf_eqtl_list_EA = list()
abf_eqtl_list = list()

MR_stats_list = list()

# Let's first obtain the list of cis-eQTLs - PPA threshold of 0.5 for cis
# This task will be divided into EA and AA individuals and processed separately
count = 1
for (i in c(1:nrow(genotype_matrix))) {
	snp = rownames(genotype_matrix)[i]
	# which genes are cis?
	cis_gene_pos = gene_positions[gene_positions$chr == paste0('chr', chr_number),]
	cis_genes = names(which((sapply(rownames(cis_gene_pos), function(x) {cis_gene_pos[x, 'start'] - cis_threshold <= snp_positions[i,'pos'] && gene_positions[x, 'end'] + cis_threshold >= snp_positions[i,'pos']}))))

	abf_eqtl_list[[i]] = calc_WABF_with_NAs(expression_matrix, cis_genes, genotype_matrix[i,], cov, W, PO)
	abf_eqtl_list_AA[[i]] = calc_WABF_with_NAs(expression_matrix_AA, cis_genes, genotype_matrix_AA[i,], cov_AA, W, PO)
	abf_eqtl_list_EA[[i]] = calc_WABF_with_NAs(expression_matrix_EA, cis_genes, genotype_matrix_EA[i,], cov_EA, W, PO)


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
