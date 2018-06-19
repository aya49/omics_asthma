##################################################
#  GGR_validation_pipeline.R
#
#  $proj/Scripts/causality/GGR/gtex/GGR_validation_pipeline.R
# 
#  GGR validation pipeline with trans-eQTLs and MR implementations.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# # Example
# args = c(1:8)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/arterytibial_v6p_consortium_autosomes_normalized.txt'
# args[2] = '47'
# args[3] = '2.5e8'
# args[4] = 'arterytibial'
# args[5] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
# args[6] = '/tigress/BEE/RNAseq/Output/causality/GGR/output/0mean-unnormalized_g-null_l-fdr/'
# # Number of pairs to process per job
# args[7] = '/tigress/BEE/RNAseq/Data/Networks/GGR/retrofitted/0mean-unnormalized_g-null_l-fdr/prot2TPM-er-reps_0mean-unnormalized_g-null_l-fdr-0.05_enet-2-union-network.txt'
# args[8] = '56'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
geno_option = 'continuous'
part_number = as.numeric(args[2])
cis_dist = as.numeric(args[3])
tissue_name = args[4]
cov_dir = args[5]
output_dir = args[6]
gene_pair_file = args[7]
total_parts = as.numeric(args[8])

library(dplyr)
# Import MatrixEQTL
library(MatrixEQTL)
source(paste0(proj_dir, '/Scripts/eqtls/MatrixEQTL_wrapper.R'))

# MR test statistic function
calculate_t_MR = function(z, x, y, beta_xz, beta_yz) {
  beta_mr = beta_yz/beta_xz
  sig_sq = (y - x*beta_mr) %*% (y - x*beta_mr) / (length(z) - 3)
  var_beta = sig_sq * (z %*% z) / (x %*% z)^2
  return(beta_mr/var_beta)
}

# Read in the gene pairs
input_gene_pair = read.csv(file = gene_pair_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt')

gene_positions = read.table(file = gene_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gene_positions = gene_positions[,c('gene_id', 'chr', 'start', 'end')]
gene_positions = mutate(gene_positions, chr = paste("chr", chr, sep = ""))
gene_positions = select(gene_positions, c(gene_id, chr, start, end))

# Since we're working with different versions of gene annotations, we need an additional step to only take the ENSG number
gene_positions$gene_id = sapply(gene_positions$gene_id, function(x) {strsplit(x, '\\.')[[1]][1]})
rownames(gene_positions) = gene_positions$gene_id

# input_gene_pair$Cause = sapply(input_gene_pair$Cause, function(x) {strsplit(x, '\\.')[[1]][1]})
# input_gene_pair$Effect = sapply(input_gene_pair$Effect, function(x) {strsplit(x, '\\.')[[1]][1]})
# input_gene_pair = mutate(input_gene_pair, Cause.Effect = paste(Cause, Effect, sep = "-"))

# v19_v22_map = read.table(paste0(proj_dir, '/Data/Resources/annotations/gencode.v19.v22.map.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# input_gene_pair$Cause_v19 = input_gene_pair$Cause
# input_gene_pair$Effect_v19 = input_gene_pair$Effect

# input_gene_pair$geneChr_cis = as.character(sapply(input_gene_pair$Cause_v19, function(x) {gene_positions$chr[which(gene_positions$gene_id == x)]}))
# input_gene_pair$startPos_cis = as.numeric(sapply(input_gene_pair$Cause_v19, function(x) {gene_positions$start[which(gene_positions$gene_id == x)]}))


# input_gene_pair$Cause_v19 = as.character(sapply(input_gene_pair$Cause, function(x) {v19_v22_map$v19[which(v19_v22_map$v22 == x)]}))
# input_gene_pair$Effect_v19 = as.character(sapply(input_gene_pair$Effect, function(x) {v19_v22_map$v19[which(v19_v22_map$v22 == x)]}))

# input_gene_pair$geneChr = as.character(sapply(input_gene_pair$Cause_v19, function(x) {gene_positions$chr[which(gene_positions$gene_id == x)]}))
# input_gene_pair$startPos = as.numeric(sapply(input_gene_pair$Cause_v19, function(x) {gene_positions$start[which(gene_positions$gene_id == x)]}))
input_gene_pair = input_gene_pair[order(input_gene_pair$geneChr_cis, input_gene_pair$startPos_cis),]
# only take autosomes
input_gene_pair = input_gene_pair[(input_gene_pair$geneChr_cis!='X'),]
input_gene_pair = input_gene_pair[(input_gene_pair$geneChr_trans!='X'),]
input_gene_pair = input_gene_pair[(input_gene_pair$geneChr_cis!='Y'),]
input_gene_pair = input_gene_pair[(input_gene_pair$geneChr_trans!='Y'),]

# Select the subset to work on:
num_pairs = ceiling(nrow(input_gene_pair)/total_parts)

if (part_number==total_parts) {
  input_gene_pair = input_gene_pair[(((part_number-1)*num_pairs+1):nrow(input_gene_pair)),]
} else {
  input_gene_pair = input_gene_pair[(((part_number-1)*num_pairs+1):((part_number)*num_pairs)),]
}

expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
rownames(expression_matrix) = sapply(rownames(expression_matrix), function(x) {strsplit(x, '\\.')[[1]][1]})

# Load in the covariates
suffix = '_Analysis.covariates.txt'
covars = read.csv(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

# Not all genes are in the expression matrix, since it had its filters - only take the ones that are in the matrix
input_gene_pair_test = input_gene_pair[(input_gene_pair$Cause_v19 %in% rownames(expression_matrix)),]
input_gene_pair_test = input_gene_pair_test[(input_gene_pair_test$Effect_v19 %in% rownames(expression_matrix)),]
input_gene_pair_test = input_gene_pair_test[(input_gene_pair_test$geneChr_cis != input_gene_pair_test$geneChr_trans),]

# Get the list of genes
cis_gene_list = unique(input_gene_pair_test$Cause_v19)
# trans_gene_list = unique(input_gene_pair_test$Effect_v19)

# Load up the genotype matrix
current_chr = input_gene_pair_test$geneChr_cis[1]
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', current_chr, '_Final.txt', sep="")
SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', current_chr, '.txt', sep="")

genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(snp_positions) = snp_positions$rsID

genotype_matrix = genotype_matrix[,colnames(expression_matrix)]
snp_positions = snp_positions[rownames(genotype_matrix),]

# Initialize data frames
MR_table_master = data.frame()

# Only need to save trans eQTLs in this case
empty_df = data.frame(snps = character(), gene = character(), statistic = numeric(), pvalue = numeric(), FDR = numeric(), beta = numeric(), snp_chr = character(), snp_pos = numeric(), gene_chr = character(), gene_start = numeric(), gene_end = numeric(), cis_gene = character())
cumulative_me_total = list(trans = list(ntests = 0, hist.counts = rep(0, 100), eqtls = empty_df))
cumulative_me_permute_total = list(trans = list(ntests = 0, hist.counts = rep(0, 100), eqtls = empty_df))

for (cis_gene in cis_gene_list) {
  print(cis_gene)
  rows = which(input_gene_pair_test$Cause_v19 == cis_gene)
  trans_gene_list = input_gene_pair_test$Effect_v19[rows]
  expression_matrix_part = expression_matrix[c(cis_gene, trans_gene_list),]
  gene_positions_part = gene_positions[rownames(expression_matrix_part),]

  new_chr = input_gene_pair_test$geneChr_cis[rows[1]]
  if (new_chr != current_chr) {
    # Load the new genotype matrix
    current_chr = new_chr
    genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', current_chr, '_Final.txt', sep="")
    SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', current_chr, '.txt', sep="")

    genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    rownames(snp_positions) = snp_positions$rsID

    genotype_matrix = genotype_matrix[,colnames(expression_matrix)]
    snp_positions = snp_positions[rownames(genotype_matrix),]
  }

  TSS = input_gene_pair_test$startPos_cis[rows[1]]
  snp_inds = abs(snp_positions$pos - TSS) <= 20000

  if (sum(snp_inds) > 0) {

    genotype_matrix_part = genotype_matrix[snp_inds,]
    snp_positions_part = snp_positions[snp_inds,]

    # Regress out covariates - Need this step for MR
    X = t(as.matrix(covars) - rowMeans(covars))
    y = t(as.matrix(expression_matrix_part) - rowMeans(expression_matrix_part))
    g = t(as.matrix(genotype_matrix_part) - rowMeans(genotype_matrix_part))

    res_y = y - (X %*% solve(t(X) %*% X) %*% t(X) %*% y)
    res_g = g - (X %*% solve(t(X) %*% X) %*% t(X) %*% g)

    geno_res = data.frame(t(res_g))
    exp_res = data.frame(t(res_y))

    cumulative_me = MatrixEQTL_wrapper(geno_res, exp_res, snp_positions_part, gene_positions_part, pvThresh = 1, pvThresh_cis = 1)
    # Run permuted version
    set.seed(111)
    exp_res_perm = exp_res[,sample(ncol(exp_res))]
    cumulative_me_perm = MatrixEQTL_wrapper(geno_res, exp_res_perm, snp_positions_part, gene_positions_part, pvThresh = 1, pvThresh_cis = 1)

    snps = as.character(cumulative_me$cis$eqtls$snps)

    # cumulative_me_total$cis$ntests = cumulative_me_total$cis$ntests + cumulative_me$cis$ntests
    # cumulative_me_total$cis$hist.counts = cumulative_me_total$cis$hist.counts + cumulative_me$cis$hist.counts
    cumulative_me_total$trans$ntests = cumulative_me_total$trans$ntests + cumulative_me$trans$ntests
    cumulative_me_total$trans$hist.counts = cumulative_me_total$trans$hist.counts + cumulative_me$trans$hist.counts
    # cumulative_me_permute_total$cis$ntests = cumulative_me_permute_total$cis$ntests + cumulative_me_perm$cis$ntests
    # cumulative_me_permute_total$cis$hist.counts = cumulative_me_permute_total$cis$hist.counts + cumulative_me_perm$cis$hist.counts
    cumulative_me_permute_total$trans$ntests = cumulative_me_permute_total$trans$ntests + cumulative_me_perm$trans$ntests
    cumulative_me_permute_total$trans$hist.counts = cumulative_me_permute_total$trans$hist.counts + cumulative_me_perm$trans$hist.counts

    for (i in rows) {
      print(i)
      # Collect the best trans eQTLs for true and permuted tests
      trans_gene = input_gene_pair_test$Effect_v19[i]
      eqtl_row = cumulative_me$trans$eqtls[which(cumulative_me$trans$eqtls$gene == trans_gene)[1],]
      eqtl_row$cis_gene = cis_gene
      cumulative_me_total$trans$eqtls = rbind(cumulative_me_total$trans$eqtls, eqtl_row)

      eqtl_row = cumulative_me_perm$trans$eqtls[which(cumulative_me_perm$trans$eqtls$gene == trans_gene)[1],]
      eqtl_row$cis_gene = cis_gene
      cumulative_me_permute_total$trans$eqtls = rbind(cumulative_me_permute_total$trans$eqtls, eqtl_row)

      # Collect the MR statistics
      # Take snps that are associated with trans-gene, and all cis-genes with it
      MR_table = data.frame(snp = snps)
      rownames(MR_table) = snps
      MR_table$r_sq = cor(as.numeric(exp_res[cis_gene,]), as.numeric(exp_res[trans_gene,]))
      temp_df = cumulative_me$cis$eqtls
      rownames(temp_df) = temp_df$snps
      MR_table$beta_xz = as.numeric(temp_df[snps,'beta'])
      temp_df = cumulative_me$trans$eqtls[(as.character(cumulative_me$trans$eqtls$gene) == trans_gene),]
      rownames(temp_df) = temp_df$snps
      MR_table$beta_yz = as.numeric(temp_df[snps,'beta'])
      temp_df = cumulative_me_perm$trans$eqtls[(as.character(cumulative_me_perm$trans$eqtls$gene) == trans_gene),]
      rownames(temp_df) = temp_df$snps
      MR_table$beta_yz_perm = as.numeric(temp_df[snps,'beta'])

      MR_table$t_MR = sapply(as.character(MR_table$snp), function(x) {calculate_t_MR(as.numeric(geno_res[x,]), as.numeric(exp_res[cis_gene,]), as.numeric(exp_res[trans_gene,]), MR_table[x,'beta_xz'], MR_table[x,'beta_yz'])})
      MR_table$t_MR_perm = sapply(as.character(MR_table$snp), function(x) {calculate_t_MR(as.numeric(geno_res[x,]), as.numeric(exp_res[cis_gene,]), as.numeric(exp_res_perm[trans_gene,]), MR_table[x,'beta_xz'], MR_table[x,'beta_yz_perm'])})

      max_row = MR_table[which(abs(MR_table$t_MR) == max(abs(MR_table$t_MR)))[1],]
      max_row_perm = MR_table[which(abs(MR_table$t_MR_perm) == max(abs(MR_table$t_MR_perm)))[1],]
      max_row$t_MR_perm = max_row_perm$t_MR_perm
      max_row$beta_yz_perm = max_row_perm$beta_yz_perm
      max_row$beta_xz_perm = max_row_perm$beta_xz

      if (nrow(MR_table_master) == 0) {
        MR_table_master = cbind(input_gene_pair_test[i,], max_row)
      } else {
        MR_table_master = rbind(MR_table_master, cbind(input_gene_pair_test[i,], max_row))
      }
    }
  }
}

# Accumulate and save the final results
split = strsplit(gene_pair_file, '/')[[1]]
filename = paste0(output_dir, strsplit(split[length(split)], '\\.txt')[[1]][1], '_part_', part_number, '_', tissue_name, '.RData')
save(cumulative_me_total, cumulative_me_permute_total, MR_table_master, file = filename)
