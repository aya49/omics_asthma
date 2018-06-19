##################################################
#  wabf_calculate_MR_stats.R
#
#  $proj/Scripts/causality/bayes_MR/main_pipeline/wabf_calculate_MR_stats.R
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
# args = c(1:5)
# args[1] = 'Whole_Blood'
# args[2] = '4'
# args[3] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/prep_files/'
# args[4] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/eqtls/Whole_Blood/'
# args[5] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/MR_stats/Whole_Blood/wabf_mrstats_'

tissue_name = args[1]
# Make sure that the chromosome and part numbers are the same as the eQTL run
chr_number = args[2]
prep_dir = args[3]
wabf_out_dir = args[4]
out_file = args[5]

library(MASS)
source("/tigress/BEE/RNAseq/Scripts/causality/bayes_MR/main_pipeline/wabf_functions.R")
load(paste0(prep_dir, tissue_name, '.RData'))
# Load the results of eQTL run
cis_PPA_threshold = 0.5
trans_PPA_threshold = 0.5

# Load in the eQTL run results - combine for each chromosome
abf_cis_eqtls_temp = data.frame()
abf_cis_eqtls_AA_temp = data.frame()
abf_cis_eqtls_EA_temp = data.frame()
abf_trans_eqtls_temp = data.frame()
abf_trans_eqtls_AA_temp = data.frame()
abf_trans_eqtls_EA_temp = data.frame()

search_str = paste0('*chr_', chr_number, '_*cisPPA_', cis_PPA_threshold, '_transPPA_', trans_PPA_threshold, '.RData')
files = list.files(wabf_out_dir, pattern = glob2rx(search_str))
for (f in files) {
	load(paste0(wabf_out_dir, f))
	abf_cis_eqtls_temp = rbind(abf_cis_eqtls_temp, abf_cis_eqtls)
	abf_cis_eqtls_AA_temp = rbind(abf_cis_eqtls_AA_temp, abf_cis_eqtls_AA)
	abf_cis_eqtls_EA_temp = rbind(abf_cis_eqtls_EA_temp, abf_cis_eqtls_EA)
	abf_trans_eqtls_temp = rbind(abf_trans_eqtls_temp, abf_trans_eqtls)
	abf_trans_eqtls_AA_temp = rbind(abf_trans_eqtls_AA_temp, abf_trans_eqtls_AA)
	abf_trans_eqtls_EA_temp = rbind(abf_trans_eqtls_EA_temp, abf_trans_eqtls_EA)
}

abf_cis_eqtls = abf_cis_eqtls_temp
abf_cis_eqtls_AA = abf_cis_eqtls_AA_temp
abf_cis_eqtls_EA = abf_cis_eqtls_EA_temp
abf_trans_eqtls = abf_trans_eqtls_temp
abf_trans_eqtls_AA = abf_trans_eqtls_AA_temp
abf_trans_eqtls_EA = abf_trans_eqtls_EA_temp

# Read in the genotype positions - both for MAF 1% and 5%
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_05_ld_pruned.RData')
load(genotype_file_name)
genotypes_all = genotype_matrix_master
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_01_ld_pruned.RData')
load(genotype_file_name)
genotypes_all = rbind(genotypes_all, genotype_matrix_master)
# This loads in the data frame "genotype_matrix_master"

# Get the appropriate partition
genotype_matrix = genotypes_all[unique(as.character(abf_cis_eqtls$snp)),]
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

# Get the permuted versions of these matrices
set.seed(111)
expression_matrix_perm = expression_matrix[,sample(ncol(expression_matrix))]
set.seed(111)
expression_matrix_AA_perm = expression_matrix_AA[,sample(ncol(expression_matrix_AA))]
set.seed(111)
expression_matrix_EA_perm = expression_matrix_EA[,sample(ncol(expression_matrix_EA))]

# Compile the list of trios to test
rownames(abf_cis_eqtls) = paste(abf_cis_eqtls$snp, abf_cis_eqtls$gene, sep='_')
rownames(abf_cis_eqtls_AA) = paste(abf_cis_eqtls_AA$snp, abf_cis_eqtls_AA$gene, sep='_')
rownames(abf_cis_eqtls_EA) = paste(abf_cis_eqtls_EA$snp, abf_cis_eqtls_EA$gene, sep='_')
rownames(abf_trans_eqtls) = paste(abf_trans_eqtls$snp, abf_trans_eqtls$gene, sep='_')
rownames(abf_trans_eqtls_AA) = paste(abf_trans_eqtls_AA$snp, abf_trans_eqtls_AA$gene, sep='_')
rownames(abf_trans_eqtls_EA) = paste(abf_trans_eqtls_EA$snp, abf_trans_eqtls_EA$gene, sep='_')

MR_stats_list = list()
MR_stats_list_AA = list()
MR_stats_list_EA = list()
print(dim(abf_cis_eqtls))
for (i in c(1:nrow(abf_cis_eqtls))) {
	snp = as.character(abf_cis_eqtls$snp)[i]
	cis_gene = as.character(abf_cis_eqtls$gene)[i]
	trans_genes = as.character(abf_trans_eqtls$gene[which(as.character(abf_trans_eqtls$snp) == abf_cis_eqtls$snp[i])])
	trans_beta = abf_trans_eqtls[paste(snp, trans_genes, sep='_'),]$beta
	if (length(trans_genes) > 0) {
		MR_stats_list[[i]] = data.frame(snp = snp, cis_gene = cis_gene, trans_gene = trans_genes, cis_beta = abf_cis_eqtls$beta[i], trans_beta = trans_beta, stringsAsFactors = F)
	}
	if (paste(snp, cis_gene, sep='_') %in% rownames(abf_cis_eqtls_AA)) {
		trans_genes_AA = as.character(abf_trans_eqtls_AA[paste(snp, trans_genes, sep='_'),]$gene)
		trans_genes_AA = trans_genes_AA[which(!is.na(trans_genes_AA))]
		if (length(trans_genes_AA) > 0) {
			trans_beta = abf_trans_eqtls_AA[paste(snp, trans_genes, sep='_'),]$beta
			MR_stats_list_AA[[i]] = data.frame(snp = snp, cis_gene = cis_gene, trans_gene = trans_genes_AA, cis_beta = abf_cis_eqtls_AA$beta[i], trans_beta = trans_beta, stringsAsFactors = F)
		}
	}
	if (paste(snp, cis_gene, sep='_') %in% rownames(abf_cis_eqtls_EA)) {
		trans_genes_EA = as.character(abf_trans_eqtls_EA[paste(snp, trans_genes, sep='_'),]$gene)
		trans_genes_EA = trans_genes_EA[which(!is.na(trans_genes_EA))]
		if (length(trans_genes_EA) > 0) {
			trans_beta = abf_trans_eqtls_EA[paste(snp, trans_genes, sep='_'),]$beta
			MR_stats_list_EA[[i]] = data.frame(snp = snp, cis_gene = cis_gene, trans_gene = trans_genes_EA, cis_beta = abf_cis_eqtls_EA$beta[i], trans_beta = trans_beta, stringsAsFactors = F)
		}
	}
}

# Make the list lengths the same
if (length(MR_stats_list) > length(MR_stats_list_AA)) {MR_stats_list_AA = MR_stats_list_AA[c(1:length(MR_stats_list))]}
if (length(MR_stats_list) > length(MR_stats_list_EA)) {MR_stats_list_EA = MR_stats_list_EA[c(1:length(MR_stats_list))]}
# MR_stats_df = do.call('rbind', MR_stats_list)

# Parameter settings
# For now, let use the W value of 0.15^2 - roughly translating to 95% chance that the relative risk is between 2/3 and 3/2
trans_risk_1 = 1.5
W_MR_1 = (log(trans_risk_1)/1.96)^2
trans_risk_2 = 1.1
W_MR_2 = (log(trans_risk_2)/1.96)^2
# parameters for MR-ABF
snp_pi_1 = 1e-3 # expecting 1 out of ~1000 true trans-eQTLs among chosen cis-eQTLs
gene_pi_1 = 1e-3 # expecting nonzero contribution in ~1000 trans genes

print(length(MR_stats_list))
print(length(MR_stats_list_AA))
print(length(MR_stats_list_EA))

MR_stats_list_result = list()
MR_stats_list_AA_result = list()
MR_stats_list_EA_result = list()

# for (j in c(1:10)) {
for (j in c(1:length(MR_stats_list))) {
	print(j)
	if (is.null(MR_stats_list[[j]])) {next}
	snp = MR_stats_list[[j]]$snp[1]
	cis_gene = MR_stats_list[[j]]$cis_gene[1]

	freq_MR = calc_freq_MR_with_NAs(expression_matrix[cis_gene,], expression_matrix[MR_stats_list[[j]]$trans_gene,], expression_matrix_perm[MR_stats_list[[j]]$trans_gene,], genotype_matrix[snp,], cov, beta_xz = MR_stats_list[[j]]$cis_beta[1], beta_yz = MR_stats_list[[j]]$trans_beta)
	df_0 = cbind(MR_stats_list[[j]], freq_MR)
	if(!is.null(MR_stats_list_AA[[j]])) {
		freq_MR = calc_freq_MR_with_NAs(expression_matrix_AA[cis_gene,], expression_matrix_AA[MR_stats_list_AA[[j]]$trans_gene,], expression_matrix_AA_perm[MR_stats_list_AA[[j]]$trans_gene,], genotype_matrix_AA[snp,], cov_AA, beta_xz = MR_stats_list_AA[[j]]$cis_beta[1], beta_yz = MR_stats_list_AA[[j]]$trans_beta)
		df_1 = cbind(MR_stats_list_AA[[j]], freq_MR)
	}
	if(!is.null(MR_stats_list_EA[[j]])) {
		freq_MR = calc_freq_MR_with_NAs(expression_matrix_EA[cis_gene,], expression_matrix_EA[MR_stats_list_EA[[j]]$trans_gene,], expression_matrix_EA_perm[MR_stats_list_EA[[j]]$trans_gene,], genotype_matrix_EA[snp,], cov_EA, beta_xz = MR_stats_list_EA[[j]]$cis_beta[1], beta_yz = MR_stats_list_EA[[j]]$trans_beta)
		df_2 = cbind(MR_stats_list_EA[[j]], freq_MR)
	}

	# Calculate Bayesian MR-ABF
	MR_ABF_df = calc_MR_ABF_with_NAs(expression_matrix[cis_gene,], expression_matrix, genotype_matrix[snp,], cov, MR_stats_list[[j]]$cis_beta[1], MR_stats_list[[j]]$trans_gene, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)
	MR_ABF_df_perm = calc_MR_ABF_with_NAs(expression_matrix[cis_gene,], expression_matrix_perm, genotype_matrix[snp,], cov, MR_stats_list[[j]]$cis_beta[1], MR_stats_list[[j]]$trans_gene, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)
	if (!is.null(MR_ABF_df) && !is.null(MR_ABF_df_perm)) {
		colnames(MR_ABF_df_perm) = c('beta_trans_perm', 'theta_perm', 'H00_ABF_perm', 'H01_ABF_perm', 'H10_ABF_perm', 'H11_ABF_perm', 'MR_PPA_perm')
		total_df = cbind(MR_ABF_df, MR_ABF_df_perm)
		df_0 = cbind(df_0, total_df)
		MR_stats_list_result[[j]] = df_0
	}

	if(!is.null(MR_stats_list_AA[[j]])) {
		MR_ABF_df = calc_MR_ABF_with_NAs(expression_matrix_AA[cis_gene,], expression_matrix_AA, genotype_matrix_AA[snp,], cov_AA, MR_stats_list_AA[[j]]$cis_beta[1], MR_stats_list_AA[[j]]$trans_gene, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)
		MR_ABF_df_perm = calc_MR_ABF_with_NAs(expression_matrix_AA[cis_gene,], expression_matrix_AA_perm, genotype_matrix_AA[snp,], cov_AA, MR_stats_list_AA[[j]]$cis_beta[1], MR_stats_list_AA[[j]]$trans_gene, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)
		if (!is.null(MR_ABF_df) && !is.null(MR_ABF_df_perm)) {
			colnames(MR_ABF_df_perm) = c('beta_trans_perm', 'theta_perm', 'H00_ABF_perm', 'H01_ABF_perm', 'H10_ABF_perm', 'H11_ABF_perm', 'MR_PPA_perm')
			df_1 = cbind(df_1, total_df)
			MR_stats_list_AA_result[[j]] = df_1
		}
	}
	if(!is.null(MR_stats_list_EA[[j]])) {
		MR_ABF_df = calc_MR_ABF_with_NAs(expression_matrix_EA[cis_gene,], expression_matrix_EA, genotype_matrix_EA[snp,], cov_EA, MR_stats_list_EA[[j]]$cis_beta[1], MR_stats_list_EA[[j]]$trans_gene, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)
		MR_ABF_df_perm = calc_MR_ABF_with_NAs(expression_matrix_EA[cis_gene,], expression_matrix_EA_perm, genotype_matrix_EA[snp,], cov_EA, MR_stats_list_EA[[j]]$cis_beta[1], MR_stats_list_EA[[j]]$trans_gene, W_MR_1, W_MR_2, snp_pi_1, gene_pi_1)
		if (!is.null(MR_ABF_df) && !is.null(MR_ABF_df_perm)) {
			colnames(MR_ABF_df_perm) = c('beta_trans_perm', 'theta_perm', 'H00_ABF_perm', 'H01_ABF_perm', 'H10_ABF_perm', 'H11_ABF_perm', 'MR_PPA_perm')
			df_2 = cbind(df_2, total_df)
			MR_stats_list_EA_result[[j]] = df_2
		}
	}
}

MR_stats = do.call('rbind', MR_stats_list_result)
rnames = paste(MR_stats$snp, MR_stats$trans_gene, sep='_')
MR_stats[,c('trans_stat', 'trans_abf', 'trans_ppa')] = abf_trans_eqtls[rnames, c('stat', 'abf', 'ppa')]
rnames = paste(MR_stats$snp, MR_stats$cis_gene, sep='_')
MR_stats[,c('cis_stat', 'cis_abf', 'cis_ppa', 'maf', 'snp_var_exp')] = abf_cis_eqtls[rnames, c('stat', 'abf', 'ppa', 'maf', 'snp_var_exp')]

MR_stats_AA = do.call('rbind', MR_stats_list_AA_result)
rnames = paste(MR_stats_AA$snp, MR_stats_AA$trans_gene, sep='_')
MR_stats_AA[,c('trans_stat', 'trans_abf', 'trans_ppa')] = abf_trans_eqtls_AA[rnames, c('stat', 'abf', 'ppa')]
rnames = paste(MR_stats_AA$snp, MR_stats_AA$cis_gene, sep='_')
MR_stats_AA[,c('cis_stat', 'cis_abf', 'cis_ppa', 'maf', 'snp_var_exp')] = abf_cis_eqtls_AA[rnames, c('stat', 'abf', 'ppa', 'maf', 'snp_var_exp')]

MR_stats_EA = do.call('rbind', MR_stats_list_EA_result)
rnames = paste(MR_stats_EA$snp, MR_stats_EA$trans_gene, sep='_')
MR_stats_EA[,c('trans_stat', 'trans_abf', 'trans_ppa')] = abf_trans_eqtls_EA[rnames, c('stat', 'abf', 'ppa')]
rnames = paste(MR_stats_EA$snp, MR_stats_EA$cis_gene, sep='_')
MR_stats_EA[,c('cis_stat', 'cis_abf', 'cis_ppa', 'maf', 'snp_var_exp')] = abf_cis_eqtls_EA[rnames, c('stat', 'abf', 'ppa', 'maf', 'snp_var_exp')]

sav = paste0(out_file, 'chr_', chr_number, '_cisPPA_', cis_PPA_threshold, '_transPPA_', trans_PPA_threshold, '_trans_risk_1_', trans_risk_1, '_trans_risk_2_', trans_risk_2, '_snp_pi1_', snp_pi_1, '_gene_pi1_', gene_pi_1, '.RData')
save(MR_stats, MR_stats_AA, MR_stats_EA, file = sav)
