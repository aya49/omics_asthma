##################################################
#  MatrixEQTL_FDR_control_intrachromosomal_final.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/MatrixEQTL_FDR_control_intrachromosomal_final.R
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# FDR Control of p-values from the null distribution
# original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/intrachromosomal/'
# tissue_name = 'vagina'
# write_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/intrachromosomal/'

original_dir = args[1]
tissue_name = args[2]
write_dir = args[3]

# library(ggplot2)
library(dplyr)

list_files = list.files(path=paste(original_dir, tissue_name, '/', sep=''))
print(length(list_files))

# Required data to aggregate
cumulative_me_trans_total = vector('list', length(list_files))
cumulative_me_cis_total = vector('list', length(list_files))

# cumulative_me_trans_permute_total = vector('list', length(list_files))
all_hist_vals_trans_total = rep(0,100)
all_hist_vals_cis_total = rep(0,100)

# all_hist_vals_trans_permute_total = rep(0,100)
all_n_tests_trans_total = 0
all_n_tests_cis_total = 0
snp_set_size_final_total = 0

i = 1
for (item in list_files) {
    load(paste(original_dir, tissue_name, '/', item, sep=""))
    
    cumulative_me_trans_total[[i]] = me$trans$eqtls
    cumulative_me_cis_total[[i]] = me$cis$eqtls

    all_hist_vals_trans_total = all_hist_vals_trans_total+me$trans$hist.counts
    all_hist_vals_cis_total = all_hist_vals_cis_total+me$cis$hist.counts

    all_n_tests_trans_total = all_n_tests_trans_total+me$trans$ntests
    all_n_tests_cis_total = all_n_tests_trans_total+me$cis$ntests
    snp_set_size_final_total = snp_set_size_final_total + snp_set_size_final

    i = i+1
}

cumulative_me_trans_total = do.call("rbind", cumulative_me_trans_total)
cumulative_me_cis_total = do.call("rbind", cumulative_me_cis_total)

# cumulative_me_trans_permute_total = do.call("rbind", cumulative_me_trans_permute_total)
print(dim(cumulative_me_trans_total))
print(dim(cumulative_me_cis_total))

BH_FDR = function(pvalues, ntests) {
    # Assumes that pvalues are ordered
    # ordered_pvals_list = pvalues[order(pvalues)]
    # alpha_hat = sapply(c(1:length(ordered_pvals_list)), function(x) {ordered_pvals_list[x] * ntests / x})

    alpha_hat = sapply(c(1:length(pvalues)), function(x) {pvalues[x] * ntests / x})
    alpha_hat[(alpha_hat >= 1.0)] = 1.0

    # Perform Benjamini-Hochberg
    alpha_hat_copy = alpha_hat
    cur_index = 1
    max_index = max(which(alpha_hat_copy<1))
    while (cur_index <= max_index) {
        temp = alpha_hat_copy[cur_index:length(alpha_hat_copy)]
        min_index = which.min(temp)
        min_alpha = min(temp)
        alpha_hat_copy[c(cur_index:(cur_index+min_index-1))] = min_alpha
        cur_index = cur_index + min_index
    }
    return(alpha_hat_copy)
}

cumulative_me_trans_total = cumulative_me_trans_total[order(cumulative_me_trans_total$pvalue),]
cumulative_me_trans_total$FDR = BH_FDR(cumulative_me_trans_total$pvalue, all_n_tests_trans_total)

cumulative_me_cis_total = cumulative_me_cis_total[order(cumulative_me_cis_total$pvalue),]
# temp - just keep pvalues < 1e-5
cumulative_me_cis_total = cumulative_me_cis_total[cumulative_me_cis_total$pvalue < 1e-5,]
cumulative_me_cis_total$FDR = BH_FDR(cumulative_me_cis_total$pvalue, all_n_tests_cis_total)

# Now we can finish the post-processing steps:
# Make necessary directories:
# dir.create(file.path(write_dir, 'pvalue_histograms'), showWarnings = FALSE)
dir.create(file.path(write_dir, 'eqtl_list'), showWarnings = FALSE)

## Generate p-value histograms
# print('Generateing p-value histograms:')
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans.png', sep=''), plot=g)
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans_permute.png', sep=''), plot=g)

# Save necessary files and write out summary
# As for actual eQTL list, write out only for pvalue<save_threshold, as requested by Yuan

save(cumulative_me_trans_total, cumulative_me_cis_total, all_hist_vals_trans_total, all_n_tests_trans_total, snp_set_size_final_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_intrachrom_list_by_pvalue.RData", sep=""))
