##################################################
#  MatrixEQTL_wrapper.R
#
#  $proj/Scripts/eqtls/MatrixEQTL_wrapper.R
# 
#  This is a wrapper function for running MatrixEQTL
#  Constraints: You need to make sure the genotype matrix and snp positions have the same set of entries, and also that
#  the expression matrix and gene positions have the same set of entries. Also, make sure the column orders of
#  genotype, expression and covariates are all the same.
#  You can specify a different p-value threshold for cis and trans eqtls.
#  This function doesn't output the QQ plot stats, but the p-value histogram stats (number of pairs that have 
#  p-values from 0 to 1, in increments of 0.01)
#  Another thing to take note: 
#  MatrixEQTL requires that cis p-value threshold is larger than or equal to the overall threshold.
#  However, setting cis threshold to 0 is okay - in this case all eqtls are collapsed into a single data frame.
#  So I customized the algorithm such that setting cis threshold to 0 will only return trans eqtls, not collapse everything.
#
#  Author: Brian Jo
#
##################################################

# Five optional parameters - setting pvThresh to 0 only return cis, while pvThresh_cis to 0 only return trans
# You can also tweak the number of bins for the histogram.

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

  # By default, quantile-normalize the expression values - ties broken by averaging
  gene_qnorm = gene$copy();
  for( sl in 1:length(gene_qnorm) ) {
    mat = gene_qnorm[[sl]];
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(gene_qnorm)+1));
    gene_qnorm[[sl]] = mat;
  }
  rm(sl, mat);

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
        gene = gene_qnorm,
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

  # Combine the data frames
  cumulative_me_null = data.frame(snps = character(0), gene = character(0), statistic = numeric(0), pvalue = numeric(0), FDR = numeric(0), beta = numeric(0), snp_chr = character(0), snp_pos = numeric(0), gene_chr = character(0), gene_start = numeric(0), gene_end = numeric(0))

  me = list()

  if (pvThresh > 0) {
    cumulative_me_trans = do.call("rbind", trans_eqtl_placeholder)
    if (!is.null(cumulative_me_trans)) {
      cumulative_me_trans = cumulative_me_trans[order(cumulative_me_trans$pvalue),]
      cumulative_me_trans$FDR = BH_FDR(cumulative_me_trans$pvalue, all_n_tests_trans)
    } else {
      cumulative_me_trans = cumulative_me_null
    }

    me$trans = list(eqtls = cumulative_me_trans, ntests = all_n_tests_trans, hist.counts = all_hist_vals_trans)
  }

  if (pvThresh_cis > 0) {
    cumulative_me_cis = do.call("rbind", cis_eqtl_placeholder)
    if (!is.null(cumulative_me_cis)) {
      cumulative_me_cis = cumulative_me_cis[order(cumulative_me_cis$pvalue),]
      cumulative_me_cis$FDR = BH_FDR(cumulative_me_cis$pvalue, all_n_tests_cis)
    } else {
      cumulative_me_cis = cumulative_me_null
    }

    me$cis = list(eqtls = cumulative_me_cis, ntests = all_n_tests_cis, hist.counts = all_hist_vals_cis)
  }

  return(me)
}

# Perform Benjamini-Hochberg
BH_FDR = function(pvalues, ntests) {
  ordered_pvals = pvalues[order(pvalues)]
  alpha_hat = sapply(c(1:length(ordered_pvals)), function(x) {ordered_pvals[x] * ntests / x})
  alpha_hat[(alpha_hat >= 1.0)] = 1.0

  alpha_hat_copy = alpha_hat
  cur_index = 1
  if (alpha_hat_copy[1] < 1) {
    max_index = max(which(alpha_hat_copy<1))
    ptm <- proc.time()
    while (cur_index <= max_index) {
      temp = alpha_hat_copy[cur_index:length(alpha_hat_copy)]
      min_index = which.min(temp)
      min_alpha = min(temp)
      alpha_hat_copy[c(cur_index:(cur_index+min_index-1))] = min_alpha
      cur_index = cur_index + min_index
    }
    print("FDR calculation time: ")
    proc.time() - ptm
  }
  return(alpha_hat_copy)
}