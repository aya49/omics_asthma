##################################################
#  trans_matrix_eqtl_seg_runs.R
#
#  $proj/Scripts/eqtls/trans/gtex/trans_matrix_eqtl_seg_runs.R
# 
#  Split up the AA/EA populations when performing eQTL analysis, to make sure this effect doesn't seriously coufound our results
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
# args = c(1:7)
# args[1] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/Final_trans_eQTL_list_0.1.txt'
# args[2] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/brainputamenbasalganglia_nonverlapping_certain_autosomes_normalized.txt'
# args[3] = '1'
# args[4] = '2.5e8'
# args[5] = 'brainputamenbasalganglia'
# args[6] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
# args[7] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/seg_runs_orthogonalize/brainputamenbasalganglia/brainputamenbasalganglia_nonverlapping_certain_autosomes_normalized_MatrixEQTL'

input_trans_list_file = args[1]
expression_file_location = args[2]
# changed feature - geno_option is always continuous by default
geno_option = 'continuous'
num_split = as.numeric(args[3])
cis_dist = as.numeric(args[4])
tissue_name = args[5]
cov_dir = args[6]
out_file = args[7]
subj_anno_file = paste0(proj_dir, '/Data/Resources/gtex/covariates/GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt')

# Import MatrixEQTL
library(MatrixEQTL)

# Need to have a custom MatrixEQTL wrapper to make sure not to qn the residuals
MatrixEQTL_wrapper = function(genotype_matrix, expression_matrix, snp_positions, gene_positions, pvThresh = 1e-2, pvThresh_cis = 1e-2, hist_bins = 100, cis_dist = 2.5e8, covariates = NULL) {
  library(MatrixEQTL)
  slice_size = 5000

  cvrt = SlicedData$new()
  if (!is.null(covariates)) {
    cvrt$CreateFromMatrix(as.matrix(covariates));
  }

  snps = SlicedData$new()
  snps$CreateFromMatrix(as.matrix(genotype_matrix));

  gene = SlicedData$new()
  gene$CreateFromMatrix(as.matrix(expression_matrix));

  useModel = modelLINEAR;

  # change pvThresh_cis to pvThresh if it is 0:
  pvThresh_cis_input = pvThresh_cis
  if (pvThresh_cis == 0) {
    pvThresh_cis_input = pvThresh
  }

  # Since this function will be primarily run in the cluster, memory constraints will be a serious consideration:
  # For now, the solution is to partition the genotypes into slices of 5000, and merge the results back into a single data frame.

  all_n_tests_cis = 0
  all_n_tests_trans = 0
  all_hist_vals_cis = rep(0,hist_bins)
  all_hist_vals_trans = rep(0,hist_bins)

  num_split = ceiling(dim(genotype_matrix)[1] / slice_size)

  cis_eqtl_placeholder = vector('list', num_split)
  trans_eqtl_placeholder = vector('list', num_split)

  for (i in c(1:num_split)) {

    num_snps = ceiling(dim(genotype_matrix)[1] / num_split)
    snp_range = c( ((i-1)*num_snps + 1) : (i*num_snps) )
    if (i == num_split) {
      snp_range = c( ((i-1)*num_snps + 1) : (dim(genotype_matrix)[1]) )
    }

    genotype_matrix_part = genotype_matrix[snp_range,]
    snp_positions_part = snp_positions[snp_range,]

    snps = SlicedData$new()
    snps$CreateFromMatrix(as.matrix(genotype_matrix_part));

    me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = NULL,
        pvOutputThreshold = pvThresh,
        useModel = useModel,
        verbose = TRUE,
        output_file_name.cis = NULL,
        pvOutputThreshold.cis = pvThresh_cis_input,
        snpspos = snp_positions_part,
        genepos = gene_positions,
        cisDist = cis_dist,
        pvalue.hist = hist_bins,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)

    # record trans- pairs
    if (pvThresh > 0) {
      if (dim(me$trans$eqtls)[1] > 0) {
        snps_list = sapply(me$trans$eqtls$snps, as.character)
        gene_list = sapply(me$trans$eqtls$gene, as.character)

        temp_df = me$trans$eqtls
        temp_df$snp_chr = snp_positions[snps_list,]$chr
        temp_df$snp_pos = snp_positions[snps_list,]$pos
        temp_df$gene_chr = gene_positions[gene_list,]$chr
        temp_df$gene_start = gene_positions[gene_list,]$start
        temp_df$gene_end = gene_positions[gene_list,]$end

        trans_eqtl_placeholder[[i]] = temp_df
        all_n_tests_trans = all_n_tests_trans + me$trans$ntests
        all_hist_vals_trans = all_hist_vals_trans + me$trans$hist.counts
      }
    }

    # record cis- pairs
    if (pvThresh_cis > 0) {
      if (dim(me$cis$eqtls)[1] > 0) {
        snps_list = sapply(me$cis$eqtls$snps, as.character)
        gene_list = sapply(me$cis$eqtls$gene, as.character)

        temp_df = me$cis$eqtls
        temp_df$snp_chr = snp_positions[snps_list,]$chr
        temp_df$snp_pos = snp_positions[snps_list,]$pos
        temp_df$gene_chr = gene_positions[gene_list,]$chr
        temp_df$gene_start = gene_positions[gene_list,]$start
        temp_df$gene_end = gene_positions[gene_list,]$end

        cis_eqtl_placeholder[[i]] = temp_df
        all_n_tests_cis = all_n_tests_cis + me$cis$ntests
        all_hist_vals_cis = all_hist_vals_cis + me$cis$hist.counts
      }
    }

  }

  me = list()
  me$trans = list(eqtls = do.call("rbind", trans_eqtl_placeholder), ntests = all_n_tests_trans, hist.counts = all_hist_vals_trans)

  return(me)
}


# Read in the list of trans eqtls
input_trans_list = read.csv(file = input_trans_list_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
input_trans_list = input_trans_list[input_trans_list$tissue == tissue_name,]

if (dim(input_trans_list)[1] == 0) {
  quit(save="no")
}

input_trans_list$identifier = sapply(c(1:dim(input_trans_list)[1]), function(x) {paste0(input_trans_list$snps[x], '_', input_trans_list$gene[x])})

# Take the set of SNPS for this run
snp_chrs = unique(input_trans_list$snp_chr)

geno_mat_list = vector('list', length(snp_chrs))
snp_pos_list = vector('list', length(snp_chrs))

ind = 1
for (i in snp_chrs) {
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

  geno_mat_list[[ind]] = genotype_matrix_chr
  snp_pos_list[[ind]] = snp_positions_chr
  ind = ind + 1
}

genotype_matrix = do.call("rbind", geno_mat_list)
snp_positions = do.call("rbind", snp_pos_list)

# Need to check that the subjects read in from expression matrices match the subjects in genotype data
library(dplyr)

gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt')
expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
expression_matrix = expression_matrix[input_trans_list$gene,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,colnames(expression_matrix)]

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

# First, run the analysis without any covariates for comparison:

me_trans = MatrixEQTL_wrapper(genotype_matrix, expression_matrix, snp_positions, gene_positions, pvThresh = 1, pvThresh_cis = 0)

cumulative_me_trans = me_trans$trans$eqtls
rownames(cumulative_me_trans) = sapply(c(1:dim(cumulative_me_trans)[1]), function(x) {paste(cumulative_me_trans$snps[x], cumulative_me_trans$gene[x], sep='_')})
cumulative_me_trans = cumulative_me_trans[input_trans_list$identifier, !(names(cumulative_me_trans) == 'FDR')]

cumulative_me_trans$tissue_MAF = rowSums(genotype_matrix[as.character(cumulative_me_trans$snps),])/dim(genotype_matrix)[2]/2
cumulative_me_trans$tissue_MAF = sapply(cumulative_me_trans$tissue_MAF, function(x) {min(x, 1-x)})

write.table(cumulative_me_trans, file=paste0(out_file, '.txt.nocov'), quote = FALSE, sep = "\t", row.names = FALSE)

# Load in the list of EA/AA individuals
subj_annotation = read.csv(file = subj_anno_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

EA_indiv = subj_annotation$SUBJID[subj_annotation$RACE == "3"]
EA_indiv = as.character(sapply(EA_indiv, function(x) {paste(strsplit(x, '-')[[1]][1], strsplit(x, '-')[[1]][2], sep='.')}))
AA_indiv = subj_annotation$SUBJID[subj_annotation$RACE == "2"]
AA_indiv = as.character(sapply(AA_indiv, function(x) {paste(strsplit(x, '-')[[1]][1], strsplit(x, '-')[[1]][2], sep='.')}))

suffix = '_Analysis.covariates.txt'
covars = read.table(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

X = t(as.matrix(covars) - rowMeans(covars))
y = t(as.matrix(expression_matrix) - rowMeans(expression_matrix))
g = t(as.matrix(genotype_matrix) - rowMeans(genotype_matrix))

# residualization version
# res_y = y - (X %*% solve(t(X) %*% X) %*% t(X) %*% y)
# res_g = g - (X %*% solve(t(X) %*% X) %*% t(X) %*% g)

# orthogonalization version
Q = qr.Q(qr(X))
res_y = y - Q %*% (t(Q) %*% y)
res_g = g - Q %*% (t(Q) %*% g)

geno_res = data.frame(t(res_g))
exp_res = data.frame(t(res_y))

me_trans_all = MatrixEQTL_wrapper( geno_res[,colnames(geno_res) %in% c(EA_indiv, AA_indiv)], exp_res[,colnames(exp_res) %in% c(EA_indiv, AA_indiv)], snp_positions, gene_positions, pvThresh = 1, pvThresh_cis = 0)
me_trans_ea = MatrixEQTL_wrapper( geno_res[,colnames(geno_res) %in% EA_indiv], exp_res[,colnames(exp_res) %in% EA_indiv], snp_positions, gene_positions, pvThresh = 1, pvThresh_cis = 0)
me_trans_aa = MatrixEQTL_wrapper( geno_res[,colnames(geno_res) %in% AA_indiv], exp_res[,colnames(exp_res) %in% AA_indiv], snp_positions, gene_positions, pvThresh = 1, pvThresh_cis = 0)

# AA+EA
cumulative_me_trans = me_trans_all$trans$eqtls
rownames(cumulative_me_trans) = sapply(c(1:dim(cumulative_me_trans)[1]), function(x) {paste(cumulative_me_trans$snps[x], cumulative_me_trans$gene[x], sep='_')})
cumulative_me_trans = cumulative_me_trans[input_trans_list$identifier, !(names(cumulative_me_trans) == 'FDR')]

temp_genotype_mat = genotype_matrix[as.character(cumulative_me_trans$snps),colnames(genotype_matrix) %in% c(EA_indiv, AA_indiv)]

cumulative_me_trans$tissue_MAF = rowSums(temp_genotype_mat)/dim(temp_genotype_mat)[2]/2
cumulative_me_trans$tissue_MAF = sapply(cumulative_me_trans$tissue_MAF, function(x) {min(x, 1-x)})

write.table(cumulative_me_trans, file=paste0(out_file, '.txt'), quote = FALSE, sep = "\t", row.names = FALSE)

# EA
cumulative_me_trans = me_trans_ea$trans$eqtls
rownames(cumulative_me_trans) = sapply(c(1:dim(cumulative_me_trans)[1]), function(x) {paste(cumulative_me_trans$snps[x], cumulative_me_trans$gene[x], sep='_')})
cumulative_me_trans = cumulative_me_trans[input_trans_list$identifier, !(names(cumulative_me_trans) == 'FDR')]

temp_genotype_mat = genotype_matrix[as.character(cumulative_me_trans$snps),colnames(genotype_matrix) %in% EA_indiv]

cumulative_me_trans$tissue_MAF = rowSums(temp_genotype_mat)/dim(temp_genotype_mat)[2]/2
cumulative_me_trans$tissue_MAF = sapply(cumulative_me_trans$tissue_MAF, function(x) {min(x, 1-x)})

write.table(cumulative_me_trans, file=paste0(out_file, '.txt.ea'), quote = FALSE, sep = "\t", row.names = FALSE)

# AA
cumulative_me_trans = me_trans_aa$trans$eqtls
rownames(cumulative_me_trans) = sapply(c(1:dim(cumulative_me_trans)[1]), function(x) {paste(cumulative_me_trans$snps[x], cumulative_me_trans$gene[x], sep='_')})
cumulative_me_trans = cumulative_me_trans[input_trans_list$identifier, !(names(cumulative_me_trans) == 'FDR')]

temp_genotype_mat = genotype_matrix[as.character(cumulative_me_trans$snps),colnames(genotype_matrix) %in% AA_indiv]

cumulative_me_trans$tissue_MAF = rowSums(temp_genotype_mat)/dim(temp_genotype_mat)[2]/2
cumulative_me_trans$tissue_MAF = sapply(cumulative_me_trans$tissue_MAF, function(x) {min(x, 1-x)})

write.table(cumulative_me_trans, file=paste0(out_file, '.txt.aa'), quote = FALSE, sep = "\t", row.names = FALSE)

fileConn<-file(paste0(out_file, '.nsample'))
line1 = paste('Total Run', dim(expression_matrix)[2], sep='\t')
line2 = paste('EA/AA', sum(colnames(geno_res) %in% c(EA_indiv, AA_indiv)), sep='\t')
line3 = paste('EA', sum(colnames(geno_res) %in% EA_indiv), sep='\t')
line4 = paste('AA', sum(colnames(geno_res) %in% AA_indiv), sep='\t')
writeLines(c(line1, line2, line3, line4), fileConn)
close(fileConn)
