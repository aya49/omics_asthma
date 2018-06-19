##################################################
#  MatrixEQTL_FDR_control_cov.R
#
#  $proj/Scripts/eqtls/trans/gtex/MatrixEQTL_FDR_control_cov.R
# 
#  Just obtaining the list of trans-SNPs associated with known covariates
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)

# FDR Control of p-values from the null distribution

original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/cov_association/'
tissue_name = 'thyroid'
write_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/cov_association/'

original_dir = args[1]
tissue_name = args[2]
write_dir = args[3]

library(ggplot2)

list_files = list.files(path=paste(original_dir, tissue_name, '/', sep=''),pattern='_me.RData')
print(length(list_files))

# Required data to aggregate
cumulative_me_trans_total = vector('list', length(list_files))
all_hist_vals_trans_total = rep(0,100)
all_n_tests_trans_total = 0

i = 1
for (item in list_files) {
    load(paste(original_dir, tissue_name, '/', item, sep=""))
    
    cumulative_me_trans_total[[i]] = cumulative_me_trans
    all_hist_vals_trans_total = all_hist_vals_trans_total+all_hist_vals_trans
    all_n_tests_trans_total = all_n_tests_trans_total+all_n_tests_trans

    i = i+1
}

cumulative_me_trans_total = do.call("rbind", cumulative_me_trans_total)
print(dim(cumulative_me_trans_total))

# trans processing first
ordering = order(cumulative_me_trans_total$pvalue)
ordered_pvals_list = cumulative_me_trans_total$pvalue[ordering]

cumulative_me_trans_total = cumulative_me_trans_total[ordering,]
alpha_hat = sapply(c(1:length(ordered_pvals_list)), function(x) {ordered_pvals_list[x] * all_n_tests_trans_total / x})
alpha_hat[(alpha_hat >= 1.0)] = 1.0
cumulative_me_trans_total$FDR = alpha_hat

# Perform Benjamini-Hochberg
alpha_hat_copy = alpha_hat
cur_index = 1
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

cumulative_me_trans_total$FDR = alpha_hat_copy



# Now we can finish the post-processing steps:
# Make necessary directories:
dir.create(file.path(write_dir, 'pvalue_histograms'), showWarnings = FALSE)
dir.create(file.path(write_dir, 'eqtl_list'), showWarnings = FALSE)

## Generate p-value histograms
print('Generateing p-value histograms:')
g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_trans.png', sep=''), plot=g)

# Save necessary files and write out summary
# As for actual eQTL list, write out only for pvalue<save_threshold, as requested by Yuan

save_threshold = 1e-2
null_df = data.frame(snps = character(), gene = character(), statistic = numeric(), pvalue = numeric(), FDR = numeric(), beta = numeric(), snp_chr = character(), snp_pos = numeric(), gene_chr = character(), gene_start = numeric(), gene_end = numeric(), mappability = numeric(), tissue_MAF = numeric())

if (sum(cumulative_me_trans_total$pvalue<save_threshold) > 0) {
    write.table(cumulative_me_trans_total[c(1:sum(cumulative_me_trans_total$pvalue<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, "_trans_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
} else {
    write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, "_trans_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
}

save(all_n_tests_trans_total, all_hist_vals_trans_total, cumulative_me_trans_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_run_total_data.RData", sep=""))

# What proportion of trans SNPs are represented, both for known covariates and PEER factors?



