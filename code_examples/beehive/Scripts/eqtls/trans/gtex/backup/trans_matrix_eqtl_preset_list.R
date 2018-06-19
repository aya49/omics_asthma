##################################################
#  trans_matrix_eqtl_preset_list.R
#
#  $proj/Scripts/eqtls/trans/gtex/trans_matrix_eqtl_preset_list.R
# 
#  This version is the most up-to-date version for trans- pipeline.
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
# args = c(1:7)
# args[1] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/Final_trans_eQTL_list_0.5.txt'
# args[2] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/seg_runs/normalized/arteryaorta_nonverlapping_certain_autosomes_normalized.txt.aa'
# args[3] = '1'
# args[4] = '2.5e8'
# args[5] = 'arteryaorta'
# args[6] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/seg_runs/normalized/covariates/'
# args[7] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/seg_runs/EA/arteryaorta/arteryaorta_nonverlapping_certain_autosomes_normalized_MatrixEQTL.ea'

input_trans_list_file = args[1]
expression_file_location = args[2]
# changed feature - geno_option is always continuous by default
geno_option = 'continuous'
num_split = as.numeric(args[3])
cis_dist = as.numeric(args[4])
tissue_name = args[5]
cov_dir = args[6]
out_file = args[7]

input_trans_list = read.csv(file = input_trans_list_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
input_trans_list = input_trans_list[input_trans_list$tissue == tissue_name,]

if (dim(input_trans_list)[1] == 0) {
  quit(save="no")
}

input_trans_list$identifier = sapply(c(1:dim(input_trans_list)[1]), function(x) {paste0(input_trans_list$snps[x], '_', input_trans_list$gene[x])})

# Take the set of SNPS for this run
snp_chrs = unique(input_trans_list$snp_chr)

i = snp_chrs[1]
chr_num = strsplit(i, "chr")[[1]][2]
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr_num, '_Final.txt', sep="")
SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', chr_num, '.txt', sep="")

genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(snp_positions) = snp_positions$rsID

genotype_matrix = genotype_matrix[input_trans_list$snps[input_trans_list$snps %in% rownames(genotype_matrix)],]
snp_positions = snp_positions[rownames(snp_positions) %in% rownames(genotype_matrix),]
genotype_matrix = genotype_matrix[rownames(genotype_matrix) %in% rownames(snp_positions),]

if (length(snp_chrs) > 1) {
  for (i in snp_chrs[2:length(snp_chrs)]) {
    chr_num = strsplit(i, "chr")[[1]][2]
    print(chr_num)
    genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr_num, '_Final.txt', sep="")
    SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', chr_num, '.txt', sep="")

    genotype_matrix_chr = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    snp_positions_chr = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    rownames(snp_positions_chr) = snp_positions_chr$rsID

    genotype_matrix_chr = genotype_matrix_chr[input_trans_list$snps[input_trans_list$snps %in% rownames(genotype_matrix_chr)],]
    snp_positions_chr = snp_positions_chr[rownames(snp_positions_chr) %in% rownames(genotype_matrix_chr),]
    genotype_matrix_chr = genotype_matrix_chr[rownames(genotype_matrix_chr) %in% rownames(snp_positions_chr),]

    genotype_matrix = rbind(genotype_matrix, genotype_matrix_chr)
    snp_positions = rbind(snp_positions, snp_positions_chr)
  }
}

# Need to check that the subjects read in from expression matrices match the subjects in genotype data
library(dplyr)

gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt')
expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]

# # Load in the mappability list:
# mappability_cutoff = 0.8
mappability_list = read.table(paste0(proj_dir, '/Data/Resources/annotations/avg_mappability_Exon_UTR.txt'), stringsAsFactors=F)
rownames(mappability_list) = mappability_list$V1
# # Arbitrary 0.8 cutoff - can be modified for a different threshold
# mappability_list = mappability_list[(mappability_list$V2>mappability_cutoff),]
# # Filter out genes with low mappability
# expression_matrix = expression_matrix[rownames(expression_matrix) %in% rownames(mappability_list),]

expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
expression_matrix = expression_matrix[input_trans_list$gene,]
genotype_matrix = genotype_matrix[,colnames(expression_matrix)]

# # Only for all-by-all runs - divide the chromosomes into ten because the jobs because there are a lot of tests
# snp_set_size = ceiling(dim(genotype_matrix)[1]/10)
# multiplier = as.numeric(part_number) - 10*(chr_number-1)
# if (multiplier==10) {
#   genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):dim(genotype_matrix)[1]),]
# } else {
#   genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):((multiplier)*snp_set_size)),]
# }

# Also obtain genes and snps that are in the expression and genotype files, respectively.
gene_positions = read.table(file = gene_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gene_positions = gene_positions[,c('gene_id', 'chr', 'start', 'end')]
gene_positions = mutate(gene_positions, chr = paste("chr", chr, sep = ""))
gene_positions = select(gene_positions, c(gene_id, chr, start, end))

# Make sure the list of genes and their ordering are the same
# We are also filtering out non-autosomal genes
expression_matrix = expression_matrix[rownames(expression_matrix) %in% gene_positions$gene_id,]
gene_positions = gene_positions[gene_positions$gene_id %in% rownames(expression_matrix),]
rownames(gene_positions) = gene_positions$gene_id
gene_positions = gene_positions[rownames(expression_matrix),]
gene_positions$start = as.numeric(gene_positions$start)
gene_positions$end = as.numeric(gene_positions$end)
# gene_positions = gene_positions[sapply(rownames(expression_matrix), function(x) {match(x, gene_positions$gene_id)}),]
# rownames(gene_positions) = gene_positions$gene_id

# Incorporate a tissue-specific MAF filter so that we are not testing variants that actually have a low MAF in our test set
gene_set_size = dim(expression_matrix)[1]
snp_set_size_init = dim(genotype_matrix)[1]

# MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
# genotype_matrix = genotype_matrix[MAF>0.05 & MAF<0.95,]
# snp_positions = snp_positions[MAF>0.05 & MAF<0.95,]

snp_set_size_final = dim(genotype_matrix)[1]

# if(method == 'gtex_gct_v6p') {

#   suffix = '_Analysis.covariates.txt'
#   # covar_mats = list.files(path = cov_dir, pattern = suffix)
#   # ind = which(as.character(sapply(covar_mats, function(x) {tolower(gsub(" ", "", gsub("[[:punct:]]", "", strsplit(x, suffix)[[1]][1])))})) == tissue_name)
#   covars = read.table(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#   rownames(covars) = covars$ID
#   covars = covars[,colnames(expression_matrix)]

# if (dim(covars)[1] >= dim(covars)[2]) {
#   quit(save="no")
# }

# } else {

#   covars = read.table(paste(args[7], 'Covariates/GTEx/genotype_platform_gender_covariates.txt', sep=''), sep='\t', stringsAsFactors=FALSE, header=TRUE)
#   rownames(covars) = covars$X
#   covars = covars[colnames(expression_matrix),]
#   covars = covars[,(2:6)]
#   if (length(unique(covars$GENDER)) == 1) {
#     covars = subset(covars, select = -(GENDER) )
#   }
#   if (length(unique(covars$Platform)) == 1) {
#     covars = subset(covars, select = -(Platform) )
#   }
#   covars = as.data.frame(t(covars))

# }

# Import MatrixEQTL
library(MatrixEQTL)

cvrt = SlicedData$new()
# cvrt$CreateFromMatrix(as.matrix(covars));

snps = SlicedData$new()
snps$CreateFromMatrix(as.matrix(genotype_matrix));

# all_n_tests_cis = 0
all_n_tests_trans = 0
# all_hist_vals_cis = rep(0,100)
all_hist_vals_trans = rep(0,100)

# cis_eqtl_placeholder = vector('list', num_split)
trans_eqtl_placeholder = vector('list', num_split)

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1

#for (i in c(1:2)) {
for (i in c(1:num_split)) {
  num_genes = ceiling(dim(expression_matrix)[1] / as.numeric(num_split))
  gene_range = c( ((i-1)*num_genes + 1) : (i*num_genes) )
  if (i == num_split) {
    gene_range = c( ((i-1)*num_genes + 1) : (dim(expression_matrix)[1]) )
  }

  print(gene_range[1])
  print(gene_range[length(gene_range)])

  expression_matrix_part = expression_matrix[gene_range,]
  gene_positions_part = gene_positions[gene_range,]

  useModel = modelLINEAR;

  gene = SlicedData$new()
  gene$CreateFromMatrix(as.matrix(expression_matrix_part));

  gene_qnorm = gene$copy();
  # Quantile-normalize the expression values - ties broken by averaging
  for( sl in 1:length(gene_qnorm) ) {
    mat = gene_qnorm[[sl]];
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(gene_qnorm)+1));
    gene_qnorm[[sl]] = mat;
  }
  rm(sl, mat);

  # Only for generating the P-value hist:
  me = Matrix_eQTL_main(
      snps = snps,
      gene = gene_qnorm,
      cvrt = cvrt,
      output_file_name = NULL,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      verbose = TRUE,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = pvOutputThreshold,
      snpspos = snp_positions,
      genepos = gene_positions,
      cisDist = cis_dist,
      pvalue.hist = 100,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)

  # record trans- pairs
  if (dim(me$trans$eqtls)[1] > 0) {
    snps_list = sapply(me$trans$eqtls$snps, as.character)
    snp_chr = snp_positions[snps_list,]$chr
    gene_list = sapply(me$trans$eqtls$gene, as.character)
    gene_chr = gene_positions[gene_list,]$chr

    snp_loc = snp_positions[snps_list,]$pos
    gene_start_loc = gene_positions[gene_list,]$start
    gene_end_loc = gene_positions[gene_list,]$end

    temp_df = me$trans$eqtls
    temp_df$snp_chr = snp_chr
    temp_df$snp_pos = snp_loc
    temp_df$gene_chr = gene_chr
    temp_df$gene_start = gene_start_loc
    temp_df$gene_end = gene_end_loc

    temp_df$identifier = sapply(c(1:dim(temp_df)[1]), function(x) {paste0(temp_df$snps[x], '_', temp_df$gene[x])})
    temp_df = temp_df[(temp_df$identifier %in% input_trans_list$identifier),]

    trans_eqtl_placeholder[[i]] = temp_df
  }

  # # record cis- pairs
  # if (dim(me$cis$eqtls)[1] > 0) {
  #   snps_list = sapply(me$cis$eqtls$snps, as.character)
  #   snp_chr = snp_positions[snps_list,]$chr
  #   gene_list = sapply(me$cis$eqtls$gene, as.character)
  #   gene_chr = gene_positions[gene_list,]$chr

  #   snp_loc = snp_positions[snps_list,]$pos
  #   gene_start_loc = gene_positions[gene_list,]$start
  #   gene_end_loc = gene_positions[gene_list,]$end

  #   temp_df = me$cis$eqtls
  #   temp_df$snp_chr = snp_chr
  #   temp_df$snp_pos = snp_loc
  #   temp_df$gene_chr = gene_chr
  #   temp_df$gene_start = gene_start_loc
  #   temp_df$gene_end = gene_end_loc

  #   cis_eqtl_placeholder[[i]] = temp_df
  # }

  # all_n_tests_cis = all_n_tests_cis + me$cis$ntests
  all_n_tests_trans = all_n_tests_trans + me$trans$ntests

  # all_hist_vals_cis = all_hist_vals_cis + me$cis$hist.counts
  all_hist_vals_trans = all_hist_vals_trans + me$trans$hist.counts
  rm(me)
}

cumulative_me_null = data.frame(snps = character(0), gene = character(0), statistic = numeric(0), pvalue = numeric(0), FDR = numeric(0), beta = numeric(0), snp_chr = character(0), snp_pos = numeric(0), gene_chr = character(0), gene_start = numeric(0), gene_end = numeric(0))

cumulative_me_trans = do.call("rbind", trans_eqtl_placeholder)
if (!is.null(cumulative_me_trans)) {
  ordering = order(cumulative_me_trans$pvalue)
  ordered_pval = cumulative_me_trans$pvalue[ordering]
  cumulative_me_trans = cumulative_me_trans[ordering,]
} else {
  cumulative_me_trans = cumulative_me_null
}

cumulative_me_trans$mappability = mappability_list[as.character(cumulative_me_trans$gene),'V2']
cumulative_me_trans$tissue_MAF = rowSums(genotype_matrix[as.character(cumulative_me_trans$snps),])/dim(genotype_matrix)[2]/2
cumulative_me_trans$tissue_MAF = sapply(cumulative_me_trans$tissue_MAF, function(x) {min(x, 1-x)})

# cumulative_me_cis = do.call("rbind", cis_eqtl_placeholder)
# if (!is.null(cumulative_me_cis)) {
#   ordering = order(cumulative_me_cis$pvalue)
#   ordered_pval = cumulative_me_cis$pvalue[ordering]
#   cumulative_me_cis = cumulative_me_cis[ordering,]
# } else {
#   cumulative_me_cis = cumulative_me_null
# }

# cumulative_me_cis$mappability = mappability_list[as.character(cumulative_me_cis$gene),'V2']
# cumulative_me_cis$tissue_MAF = rowSums(genotype_matrix[as.character(cumulative_me_cis$snps),])/dim(genotype_matrix)[2]/2
# cumulative_me_cis$tissue_MAF = sapply(cumulative_me_cis$tissue_MAF, function(x) {min(x, 1-x)})

# We also repeat with the permuted genotypes: however, for all-by-all, this is not for calibration but for simply showing
# that the permuted runs don't deviate much from the uniform assumption.

# # A strict comparison is hard because the permutations are not tissue-specific and it can change the MAF for particular variants
# geno_option = 'continuous_permute'
# genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr_number, '_Final.txt', sep="")
# genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
# genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]
# if (multiplier==10) {
#   genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):dim(genotype_matrix)[1]),]
# } else {
#   genotype_matrix = genotype_matrix[(((multiplier-1)*snp_set_size+1):((multiplier)*snp_set_size)),]
# }

# snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
# rownames(snp_positions) = snp_positions$rsID
# snp_positions = snp_positions[rownames(snp_positions) %in% rownames(genotype_matrix),]
# genotype_matrix = genotype_matrix[rownames(genotype_matrix) %in% rownames(snp_positions),]
# # Incorporate a tissue-specific MAF filter so that we are not testing variants that actually have a low MAF in our test set
# MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
# genotype_matrix = genotype_matrix[MAF>0.05 & MAF<0.95,]
# snp_positions = snp_positions[MAF>0.05 & MAF<0.95,]

# snps = SlicedData$new()
# snps$CreateFromMatrix(as.matrix(genotype_matrix));

# # Repeat for permuted matrices
# all_n_tests_cis_permute = 0
# all_n_tests_trans_permute = 0
# all_hist_vals_cis_permute = rep(0,100)
# all_hist_vals_trans_permute = rep(0,100)

# cis_eqtl_placeholder_permute = vector('list', num_split)
# trans_eqtl_placeholder_permute = vector('list', num_split)

# for (i in c(1:num_split)) {
#   num_genes = ceiling(dim(expression_matrix)[1] / as.numeric(num_split))
#   gene_range = c( ((i-1)*num_genes + 1) : (i*num_genes) )
#   if (i == num_split) {
#     gene_range = c( ((i-1)*num_genes + 1) : (dim(expression_matrix)[1]) )
#   }

#   print(gene_range[1])
#   print(gene_range[length(gene_range)])

#   expression_matrix_part = expression_matrix[gene_range,]
#   gene_positions_part = gene_positions[gene_range,]

#   useModel = modelLINEAR;

#   gene = SlicedData$new()
#   gene$CreateFromMatrix(as.matrix(expression_matrix_part));

#   gene_qnorm = gene$copy();
#   # Quantile-normalize the expression values - ties broken by averaging
#   for( sl in 1:length(gene_qnorm) ) {
#     mat = gene_qnorm[[sl]];
#     mat = t(apply(mat, 1, rank, ties.method = "average"));
#     mat = qnorm(mat / (ncol(gene_qnorm)+1));
#     gene_qnorm[[sl]] = mat;
#   }
#   rm(sl, mat);

#   # Only for generating the P-value hist:
#   me = Matrix_eQTL_main(
#       snps = snps,
#       gene = gene_qnorm,
#       cvrt = cvrt,
#       output_file_name = NULL,
#       pvOutputThreshold = pvOutputThreshold,
#       useModel = useModel,
#       verbose = TRUE,
#       output_file_name.cis = NULL,
#       pvOutputThreshold.cis = pvOutputThreshold,
#       snpspos = snp_positions,
#       genepos = gene_positions,
#       cisDist = cis_dist,
#       pvalue.hist = 100,
#       min.pv.by.genesnp = FALSE,
#       noFDRsaveMemory = FALSE)

#   # record trans- pairs
#   if (dim(me$trans$eqtls)[1] > 0) {
#     snps_list = sapply(me$trans$eqtls$snps, as.character)
#     snp_chr = snp_positions[snps_list,]$chr
#     gene_list = sapply(me$trans$eqtls$gene, as.character)
#     gene_chr = gene_positions[gene_list,]$chr

#     snp_loc = snp_positions[snps_list,]$pos
#     gene_start_loc = gene_positions[gene_list,]$start
#     gene_end_loc = gene_positions[gene_list,]$end

#     temp_df = me$trans$eqtls
#     temp_df$snp_chr = snp_chr
#     temp_df$snp_pos = snp_loc
#     temp_df$gene_chr = gene_chr
#     temp_df$gene_start = gene_start_loc
#     temp_df$gene_end = gene_end_loc

#     trans_eqtl_placeholder_permute[[i]] = temp_df
#   }

#   # record cis- pairs
#   if (dim(me$cis$eqtls)[1] > 0) {
#     snps_list = sapply(me$cis$eqtls$snps, as.character)
#     snp_chr = snp_positions[snps_list,]$chr
#     gene_list = sapply(me$cis$eqtls$gene, as.character)
#     gene_chr = gene_positions[gene_list,]$chr

#     snp_loc = snp_positions[snps_list,]$pos
#     gene_start_loc = gene_positions[gene_list,]$start
#     gene_end_loc = gene_positions[gene_list,]$end

#     temp_df = me$cis$eqtls
#     temp_df$snp_chr = snp_chr
#     temp_df$snp_pos = snp_loc
#     temp_df$gene_chr = gene_chr
#     temp_df$gene_start = gene_start_loc
#     temp_df$gene_end = gene_end_loc

#     cis_eqtl_placeholder_permute[[i]] = temp_df
#   }

#   all_n_tests_cis_permute = all_n_tests_cis_permute + me$cis$ntests
#   all_n_tests_trans_permute = all_n_tests_trans_permute + me$trans$ntests

#   all_hist_vals_cis_permute = all_hist_vals_cis_permute + me$cis$hist.counts
#   all_hist_vals_trans_permute = all_hist_vals_trans_permute + me$trans$hist.counts
#   rm(me)
# }

# cumulative_me_trans_permute = do.call("rbind", trans_eqtl_placeholder_permute)
# if (!is.null(cumulative_me_trans_permute)) {
#   ordering = order(cumulative_me_trans_permute$pvalue)
#   ordered_pval = cumulative_me_trans_permute$pvalue[ordering]
#   cumulative_me_trans_permute = cumulative_me_trans_permute[ordering,]
# } else {
#   cumulative_me_trans_permute = cumulative_me_null
# }

# cumulative_me_trans_permute$mappability = mappability_list[as.character(cumulative_me_trans_permute$gene),'V2']
# cumulative_me_trans_permute$tissue_MAF = rowSums(genotype_matrix[as.character(cumulative_me_trans_permute$snps),])/dim(genotype_matrix)[2]/2
# cumulative_me_trans_permute$tissue_MAF = sapply(cumulative_me_trans_permute$tissue_MAF, function(x) {min(x, 1-x)})

# cumulative_me_cis_permute = do.call("rbind", cis_eqtl_placeholder_permute)
# if (!is.null(cumulative_me_cis_permute)) {
#   ordering = order(cumulative_me_cis_permute$pvalue)
#   ordered_pval = cumulative_me_cis_permute$pvalue[ordering]
#   cumulative_me_cis_permute = cumulative_me_cis_permute[ordering,]
# } else {
#   cumulative_me_cis_permute = cumulative_me_null
# }

# cumulative_me_cis_permute$mappability = mappability_list[as.character(cumulative_me_cis_permute$gene),'V2']
# cumulative_me_cis_permute$tissue_MAF = rowSums(genotype_matrix[as.character(cumulative_me_cis_permute$snps),])/dim(genotype_matrix)[2]/2
# cumulative_me_cis_permute$tissue_MAF = sapply(cumulative_me_cis_permute$tissue_MAF, function(x) {min(x, 1-x)})

# Save the Matrix EQTL object as table
write.table(cumulative_me_trans, file=paste0(out_file, '.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
