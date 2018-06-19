##################################################
#  trans_matrix_eqtl_MR_prep.R
#
#  $proj/Scripts/causality/MR_eQTL/gtex/trans_eqtl/trans_matrix_eqtl_MR_prep.R
# 
#  Prepare the list of associations for cis and trans genes for each trans-eQTL pair
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
args[1] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/Final_trans_eQTL_list_0.5.txt'
args[2] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/thyroid_nonverlapping_certain_autosomes_normalized.txt'
args[3] = '1'
args[4] = '2.5e8'
args[5] = 'thyroid'
args[6] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
args[7] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MR_eQTL/FDR_0.5/thyroid_MR_results.txt'

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

# Simple MatrixEQTL wrapper function
MatrixEQTL_wrapper_function = function(genotype_matrix, expression_matrix, snp_positions_chr, gene_positions) {

    cvrt = SlicedData$new()
    # cvrt$CreateFromMatrix(as.matrix(covars));
    snps = SlicedData$new()
    snps$CreateFromMatrix(as.matrix(genotype_matrix));

    # Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
    pvOutputThreshold = 1

    gene = SlicedData$new()
    gene$CreateFromMatrix(as.matrix(expression_matrix));

    useModel = modelLINEAR;
    gene_qnorm = gene$copy();
    # Quantile-normalize the expression values - ties broken by averaging
    for( sl in 1:length(gene_qnorm) ) {
        mat = gene_qnorm[[sl]];
        mat = t(apply(mat, 1, rank, ties.method = "average"));
        mat = qnorm(mat / (ncol(gene_qnorm)+1));
        gene_qnorm[[sl]] = mat;
    }
    rm(sl, mat);

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
    }

    cumulative_me_null = data.frame(snps = character(0), gene = character(0), statistic = numeric(0), pvalue = numeric(0), FDR = numeric(0), beta = numeric(0), snp_chr = character(0), snp_pos = numeric(0), gene_chr = character(0), gene_start = numeric(0), gene_end = numeric(0))

    cumulative_me_trans = temp_df
    if (!is.null(cumulative_me_trans)) {
      ordering = order(cumulative_me_trans$pvalue)
      ordered_pval = cumulative_me_trans$pvalue[ordering]
      cumulative_me_trans = cumulative_me_trans[ordering,]
    } else {
      cumulative_me_trans = cumulative_me_null
    }

    # record cis- pairs
    if (dim(me$cis$eqtls)[1] > 0) {
        snps_list = sapply(me$cis$eqtls$snps, as.character)
        snp_chr = snp_positions[snps_list,]$chr
        gene_list = sapply(me$cis$eqtls$gene, as.character)
        gene_chr = gene_positions[gene_list,]$chr

        snp_loc = snp_positions[snps_list,]$pos
        gene_start_loc = gene_positions[gene_list,]$start
        gene_end_loc = gene_positions[gene_list,]$end

        temp_df = me$cis$eqtls
        temp_df$snp_chr = snp_chr
        temp_df$snp_pos = snp_loc
        temp_df$gene_chr = gene_chr
        temp_df$gene_start = gene_start_loc
        temp_df$gene_end = gene_end_loc

        temp_df$identifier = sapply(c(1:dim(temp_df)[1]), function(x) {paste0(temp_df$snps[x], '_', temp_df$gene[x])})
        temp_df = temp_df[(temp_df$identifier %in% cis_identifiers),]
    }

    cumulative_me_cis = temp_df
    if (!is.null(cumulative_me_cis)) {
      ordering = order(cumulative_me_cis$pvalue)
      ordered_pval = cumulative_me_cis$pvalue[ordering]
      cumulative_me_cis = cumulative_me_cis[ordering,]
    } else {
      cumulative_me_cis = cumulative_me_null
    }

    print(dim(cumulative_me_trans))
    return(rbind(cumulative_me_trans, cumulative_me_cis))
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
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,colnames(expression_matrix)]

# Also obtain genes and snps that are in the expression and genotype files, respectively.
gene_positions = read.table(file = gene_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gene_positions = gene_positions[,c('gene_id', 'chr', 'start', 'end')]
gene_positions = mutate(gene_positions, chr = paste("chr", chr, sep = ""))
gene_positions = select(gene_positions, c(gene_id, chr, start, end))

cis_gene_list = sapply(c(1:dim(snp_positions)[1]), function(x) {gene_positions$gene_id[(gene_positions$chr == snp_positions$chr[x]) & ((abs(gene_positions$start - snp_positions$pos[x])<150000) | (abs(gene_positions$end - snp_positions$pos[x])<150000))]})
cis_gene_array = unique(unlist(cis_gene_list))
cis_identifiers = unlist(sapply(c(1:length(cis_gene_list)), function(x) {paste(snp_positions$rsID[x], cis_gene_list[[x]], sep='_')}))

# length(unique(input_trans_list$gene))
# dim(expression_matrix[rownames(expression_matrix) %in% unique(input_trans_list$gene),])
# length(cis_gene_array)
# dim(expression_matrix[rownames(expression_matrix) %in% cis_gene_array,])
expression_matrix = expression_matrix[rownames(expression_matrix) %in% unique(c(input_trans_list$gene, cis_gene_array)),]
gene_positions = gene_positions[gene_positions$gene_id %in% rownames(expression_matrix),]
rownames(gene_positions) = gene_positions$gene_id
gene_positions = gene_positions[rownames(expression_matrix),]

suffix = '_Analysis.covariates.txt'
covars = read.table(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]
X = t(as.matrix(covars) - rowMeans(covars))
y = t(as.matrix(expression_matrix) - rowMeans(expression_matrix))
g = t(as.matrix(genotype_matrix) - rowMeans(genotype_matrix))

res_y = y - (X %*% solve(t(X) %*% X) %*% t(X) %*% y)
res_g = g - (X %*% solve(t(X) %*% X) %*% t(X) %*% g)

geno_res = data.frame(t(res_g))
exp_res = data.frame(t(res_y))

cumulative_me = MatrixEQTL_wrapper_function(geno_res, exp_res, snp_positions, gene_positions)

# Run permuted version
set.seed(111)
exp_res_perm = exp_res[,sample(dim(exp_res)[2])]
cumulative_me_perm = MatrixEQTL_wrapper_function(geno_res, exp_res_perm, snp_positions, gene_positions)

calculate_t_MR = function(z, x, y, beta_xz, beta_yz) {
  beta_mr = beta_yz/beta_xz
  sig_sq = (y - x*beta_mr) %*% (y - x*beta_mr) / (length(z) - 3)
  var_beta = sig_sq * (z %*% z) / (x %*% z)^2
  return(beta_mr/var_beta)
}

MR_table_list = vector('list', length(unique(input_trans_list$gene)))
i = 1
ind1 = cumulative_me$snp_chr == cumulative_me$gene_chr
for (gene in unique(input_trans_list$gene)) {
  print(i)
  # Take snps that are associated with trans-gene, and all cis-genes with it
  snp_set = as.character(input_trans_list[which(input_trans_list$gene == gene),]$snps)
  ind2 = sapply(cumulative_me$snps, function(x) {x %in% snp_set})
  ind = ind1 & ind2
  if (sum(ind) > 0) {
    MR_table = data.frame(snp = cumulative_me$snps[ind], cis_gene = cumulative_me$gene[ind], beta_xz = cumulative_me$beta[ind])
    MR_table$trans_gene = gene
    MR_table$beta_yz = sapply(c(1:dim(MR_table)[1]), function(x) {cumulative_me$beta[which(cumulative_me$identifier == paste(MR_table$snp[x], gene, sep='_'))]})
    MR_table$r_sq = sapply(c(1:dim(MR_table)[1]), function(x) {cor(as.numeric(exp_res[MR_table$cis_gene[x],]), as.numeric(exp_res[gene,]))})

    MR_table$beta_yz_perm = sapply(c(1:dim(MR_table)[1]), function(x) {cumulative_me_perm$beta[which(cumulative_me_perm$identifier == paste(MR_table$snp[x], gene, sep='_'))]})
    MR_table$r_sq_perm = sapply(c(1:dim(MR_table)[1]), function(x) {cor(as.numeric(exp_res[MR_table$cis_gene[x],]), as.numeric(exp_res_perm[gene,]))})

    MR_table$t_MR = sapply(c(1:dim(MR_table)[1]), function(x) {calculate_t_MR(as.numeric(geno_res[as.character(MR_table$snp)[x],]), as.numeric(exp_res[as.character(MR_table$cis_gene)[x],]), as.numeric(exp_res[gene,]), MR_table$beta_xz[x], MR_table$beta_yz[x])})
    MR_table$t_MR_perm = sapply(c(1:dim(MR_table)[1]), function(x) {calculate_t_MR(as.numeric(geno_res[as.character(MR_table$snp)[x],]), as.numeric(exp_res[as.character(MR_table$cis_gene)[x],]), as.numeric(exp_res_perm[gene,]), MR_table$beta_xz[x], MR_table$beta_yz_perm[x])})

    MR_table_list[[i]] = MR_table
  }
  i = i+1
}

MR_table_total = do.call("rbind", MR_table_list)

MR_table_total$identifier = sapply(c(1:dim(MR_table_total)[1]), function(x) {paste(MR_table_total$snp[x], MR_table_total$trans_gene[x], sep='_')})

order = order(abs(MR_table_total$t_MR), decreasing=TRUE)
MR_table_total = MR_table_total[order,]

library(qvalue)
qval = qvalue(p = empPvals(stat = abs(MR_table_total$t_MR), stat0 = abs(MR_table_total$t_MR_perm)))

MR_table_total$FDR = qval$qvalues

write.table(MR_table_total, file=out_file, quote = FALSE, sep = "\t", row.names = FALSE)



# ind = cumulative_me$snp_chr == cumulative_me$gene_chr
# MR_table = data.frame(snp = cumulative_me$snps[ind], cis_gene = cumulative_me$gene[ind], beta_xz = cumulative_me$beta[ind])
# MR_table$trans_gene = sapply(c(1:dim(MR_table)[1]), function(x) {input_trans_list$gene[which(input_trans_list$snps == MR_table$snp[x])]})
# MR_table$beta_yz = sapply(c(1:dim(MR_table)[1]), function(x) {cumulative_me$beta[which(cumulative_me$identifier == paste(MR_table$snp[x], MR_table$trans_gene[x], sep='_'))]})
# MR_table$r_sq = sapply(c(1:dim(MR_table)[1]), function(x) {cor(as.numeric(exp_res[MR_table$cis_gene[x],]), as.numeric(exp_res[MR_table$trans_gene[x],]))})


# MR_table$beta_yz_perm = sapply(c(1:dim(MR_table)[1]), function(x) {cumulative_me_perm$beta[which(cumulative_me_perm$identifier == paste(MR_table$snp[x], MR_table$trans_gene[x], sep='_'))]})
# MR_table$r_sq_perm = sapply(c(1:dim(MR_table)[1]), function(x) {cor(as.numeric(exp_res[MR_table$cis_gene[x],]), as.numeric(exp_res_perm[MR_table$trans_gene[x],]))})

# write.table(geno_res, file=paste0(out_file, '_geno_table.txt'), quote = FALSE, sep = "\t", row.names = TRUE)
# write.table(exp_res, file=paste0(out_file, '_exp_table.txt'), quote = FALSE, sep = "\t", row.names = TRUE)
# write.table(exp_res_perm, file=paste0(out_file, '_exp_table_perm.txt'), quote = FALSE, sep = "\t", row.names = TRUE)

# write.table(MR_table, file=paste0(out_file, '_MR_prep.txt'), quote = FALSE, sep = "\t", row.names = FALSE)

# write.table(MR_table, file='/Users/brian_jo/Desktop/Project/MR_eQTL/test/braincortex_MR_prep.txt' quote = FALSE, sep = "\t", row.names = FALSE)

# calculate_t_MR = function(z, x, y, beta_xz, beta_yz) {
#   beta_mr = beta_yz/beta_xz
#   sig_sq = (y - x*beta_mr) %*% (y - x*beta_mr) / (length(z) - 3)
#   var_beta = sig_sq * (z %*% z) / (x %*% z)^2
#   return(beta_mr/var_beta)
# }

# t_MR = sapply(c(1:dim(MR_table)[1]), function(x) {calculate_t_MR(as.numeric(geno_res[as.character(MR_table$snp)[x],]), as.numeric(exp_res[as.character(MR_table$cis_gene)[x],]), as.numeric(exp_res[as.character(MR_table$trans_gene)[x],]), MR_table$beta_xz[x], MR_table$beta_yz[x])})
# t_MR_perm = sapply(c(1:dim(MR_table)[1]), function(x) {calculate_t_MR(as.numeric(geno_res[as.character(MR_table$snp)[x],]), as.numeric(exp_res[as.character(MR_table$cis_gene)[x],]), as.numeric(exp_res_perm[as.character(MR_table$trans_gene)[x],]), MR_table$beta_xz[x], MR_table$beta_yz_perm[x])})

# MR_table$t_MR = t_MR
# MR_table$t_MR_perm = t_MR_perm
