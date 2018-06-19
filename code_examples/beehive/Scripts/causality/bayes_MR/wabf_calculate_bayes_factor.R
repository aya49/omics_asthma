##################################################
#  wabf_calculate_bayes_factor.R
#
#  $proj/Scripts/causality/bayes_MR/wabf_calculate_bayes_factor.R
# 
#  This script calculates the Wakefield Approximate Bayes Factors for the GTEx v8 data.
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
# args = c(1:7)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz'
# args[2] = '16'
# args[3] = '1'
# args[4] = '10000'
# args[5] = 'Whole_Blood'
# args[6] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
# args[7] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/WABF/raw/Whole_Blood/wabf_raw_output_'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
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

gene_positions = expression_matrix[,c(4,1,2,3)]
colnames(gene_positions) = c('gene_id', 'chr', 'start', 'end')
expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]

# Read in the genotype positions
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_05_ld_pruned.RData')
load(genotype_file_name)
# This loads in the data frame "genotype_matrix_master"

# Get the appropriate partition
num_parts = ceiling(nrow(genotype_matrix_master) / partition_size)
num_inds = ceiling(nrow(genotype_matrix_master) / num_parts)
if (part_number == num_parts) {
  genotype_matrix = genotype_matrix_master[c((((part_number-1)*partition_size)+1):nrow(genotype_matrix_master)),]
} else {
  genotype_matrix = genotype_matrix_master[c((((part_number-1)*partition_size)+1):(part_number*partition_size)),]
}
# There should be 838 indivs for genotype

# Make sure the columns are the same
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

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

calc_beta = function(n, i, inds) {
  exp = as.numeric(expression_matrix[n,inds])
  exp = exp - mean(exp)
  temp_cov = cov[inds,]
  temp_cov$SNP = as.numeric(genotype_matrix[i,inds])
  # Mean center
  temp_cov$SNP = temp_cov$SNP - mean(temp_cov$SNP)

  z = lm(exp ~ ., data = temp_cov)
  return(as.numeric(z$coefficients['SNP']))
}

abf_eqtl_list = list()
# For now, let use the W value of 0.15^2 - roughly translating to 95% chance that the relative risk is between 2/3 and 3/2
risk = 1.5
W = (log(risk)/1.96)^2
# Say we expect roughtly one out of 1e5 SNPs to be eQTLs
pi_1 = 1e-5
PO = (1-pi_1)/pi_1

# Now that all the data has been imported, we can calculate the WABF:
for (i in c(120:nrow(genotype_matrix))) {
  print(i)
  # Which genotypes are not NA?
  inds = !is.na(genotype_matrix[i,])
  if (length(unique(as.numeric(genotype_matrix[i,inds]))) == 1) {next}
  # Get the expression subset corresponding to available genotypes
  exp = t(expression_matrix)[as.logical(inds),]
  exp = center_colmeans(exp)
  temp_cov = cov[inds,]
  temp_cov$SNP = as.numeric(genotype_matrix[i,inds])
  # Mean center
  temp_cov$SNP = temp_cov$SNP - mean(temp_cov$SNP)
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

  # save the stats for any associations with PPA over 0.5
  if (sum(PPA >= 0.9) > 0) {
    assoc = which(PPA >= 0.9)
    abf_eqtl_list[[i]] = data.frame(snp = rownames(genotype_matrix)[i], gene = names(assoc), beta = as.numeric(betas[assoc]), z = as.numeric(Z[assoc]), shrink = r, BF = as.numeric(ABF[assoc]), PPA = as.numeric(PPA[assoc]))
  }
}

out_df = do.call('rbind', abf_eqtl_list)

save(out_df, file = paste0(out_file, 'chr', chr_number, '_part', part_number, '_risk_', risk, '_pi1_', pi_1, '.RData'))
