##################################################
#  MatrixEQTL_FDR_control.R
#
#  $proj/Scripts/eqtls/trans/gtex/MatrixEQTL_FDR_control.R
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)

# FDR Control of p-values from the null distribution

# original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/all-by-all/'
# tissue_name = 'thyroid'
# write_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/'

original_dir = args[1]
tissue_name = args[2]
write_dir = args[3]

library(ggplot2)

list_files = list.files(path=paste(original_dir, tissue_name, '/', sep=''),pattern='_me.RData')
print(length(list_files))

# Required data to aggregate
cumulative_me_trans_total = vector('list', length(list_files))
cumulative_me_cis_total = vector('list', length(list_files))
cumulative_me_trans_permute_total = vector('list', length(list_files))
cumulative_me_cis_permute_total = vector('list', length(list_files))

all_hist_vals_cis_total = rep(0,100)
all_hist_vals_trans_total = rep(0,100)
all_hist_vals_cis_permute_total = rep(0,100)
all_hist_vals_trans_permute_total = rep(0,100)

all_n_tests_cis_total = 0
all_n_tests_trans_total = 0
all_n_tests_cis_permute_total = 0
all_n_tests_trans_permute_total = 0

snp_set_size_init_total = 0
snp_set_size_final_total = 0

i = 1
for (item in list_files) {
    load(paste(original_dir, tissue_name, '/', item, sep=""))
    
    cumulative_me_trans_total[[i]] = cumulative_me_trans
    cumulative_me_cis_total[[i]] = cumulative_me_cis
    cumulative_me_trans_permute_total[[i]] = cumulative_me_trans_permute
    cumulative_me_cis_permute_total[[i]] = cumulative_me_cis_permute

    all_hist_vals_cis_total = all_hist_vals_cis_total+all_hist_vals_cis
    all_hist_vals_trans_total = all_hist_vals_trans_total+all_hist_vals_trans
    all_hist_vals_cis_permute_total = all_hist_vals_cis_permute_total+all_hist_vals_cis_permute
    all_hist_vals_trans_permute_total = all_hist_vals_trans_permute_total+all_hist_vals_trans_permute

    all_n_tests_cis_total = all_n_tests_cis_total+all_n_tests_cis
    all_n_tests_trans_total = all_n_tests_trans_total+all_n_tests_trans
    all_n_tests_cis_permute_total = all_n_tests_cis_permute_total+all_n_tests_cis_permute
    all_n_tests_trans_permute_total = all_n_tests_trans_permute_total+all_n_tests_trans_permute

    snp_set_size_init_total = snp_set_size_init_total+snp_set_size_init
    snp_set_size_final_total = snp_set_size_final_total+snp_set_size_final
    i = i+1
}

cumulative_me_trans_total = do.call("rbind", cumulative_me_trans_total)
cumulative_me_cis_total = do.call("rbind", cumulative_me_cis_total)
cumulative_me_trans_permute_total = do.call("rbind", cumulative_me_trans_permute_total)
cumulative_me_cis_permute_total = do.call("rbind", cumulative_me_cis_permute_total)
print(dim(cumulative_me_trans_total))
print(dim(cumulative_me_cis_total))
print(dim(cumulative_me_trans_permute_total))
print(dim(cumulative_me_cis_permute_total))

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

# trans processing, permuted
ordering = order(cumulative_me_trans_permute_total$pvalue)
ordered_pvals_list = cumulative_me_trans_permute_total$pvalue[ordering]

cumulative_me_trans_permute_total = cumulative_me_trans_permute_total[ordering,]
alpha_hat = sapply(c(1:length(ordered_pvals_list)), function(x) {ordered_pvals_list[x] * all_n_tests_trans_permute_total / x})
alpha_hat[(alpha_hat >= 1.0)] = 1.0
cumulative_me_trans_permute_total$FDR = alpha_hat

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

cumulative_me_trans_permute_total$FDR = alpha_hat_copy

# also do cis processing now:
ordering = order(cumulative_me_cis_total$pvalue)
ordered_pvals_list = cumulative_me_cis_total$pvalue[ordering]

cumulative_me_cis_total = cumulative_me_cis_total[ordering,]
alpha_hat = sapply(c(1:length(ordered_pvals_list)), function(x) {ordered_pvals_list[x] * all_n_tests_cis_total / x})
alpha_hat[(alpha_hat >= 1.0)] = 1.0
cumulative_me_cis_total$FDR = alpha_hat

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

cumulative_me_cis_total$FDR = alpha_hat_copy

# also do cis processing, permuted
ordering = order(cumulative_me_cis_permute_total$pvalue)
ordered_pvals_list = cumulative_me_cis_permute_total$pvalue[ordering]

cumulative_me_cis_permute_total = cumulative_me_cis_permute_total[ordering,]
alpha_hat = sapply(c(1:length(ordered_pvals_list)), function(x) {ordered_pvals_list[x] * all_n_tests_cis_permute_total / x})
alpha_hat[(alpha_hat >= 1.0)] = 1.0
cumulative_me_cis_permute_total$FDR = alpha_hat

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

cumulative_me_cis_permute_total$FDR = alpha_hat_copy

# Now we can finish the post-processing steps:
# Make necessary directories:
dir.create(file.path(write_dir, 'partial_qq_plots'), showWarnings = FALSE)
dir.create(file.path(write_dir, 'pvalue_histograms'), showWarnings = FALSE)
dir.create(file.path(write_dir, 'eqtl_list'), showWarnings = FALSE)
dir.create(file.path(write_dir, 'run_summary'), showWarnings = FALSE)

## Generate partial qq plot (for 1e-5 and above, for example)
num_to_plot = min(sum(cumulative_me_trans_permute_total$pvalue<1e-5), sum(cumulative_me_trans_total$pvalue<1e-5))
g = ggplot(data.frame(rand = cumulative_me_trans_permute_total$pvalue[1:num_to_plot], true = cumulative_me_trans_total$pvalue[1:num_to_plot], fdr_pass = (cumulative_me_trans_total$FDR[1:num_to_plot]<0.1)), aes(-log10(rand), -log10(true))) + geom_point(aes(colour = fdr_pass)) + scale_color_manual(values=c("black", "red")) + theme(legend.position="none")
#g = g + geom_abline(slope = 1, intercept = 0) + ggtitle(paste(tissue_name, ' partial qq plot', sep='')) + xlab('Permute variants') + ylab('cis-SNPs')
g = g + geom_abline(slope = 1, intercept = 0) + xlab('-log10(P-value) of association for permuted variants') + ylab('-log10(P-value) of association for trans-variants')
ggsave(g, filename = paste(write_dir, 'partial_qq_plots/', tissue_name, '_partial_qq_plot.png', sep=''))

## Generate p-value histograms
print('Generateing p-value histograms:')
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_cis_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('p-value distribution') + ggtitle(paste('Histogram of', tissue_name, 'with', all_n_tests_cis_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_cis.png', sep=''), plot=g)
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('p-value distribution') + ggtitle(paste('Histogram of', tissue_name, 'with', all_n_tests_trans_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_trans.png', sep=''), plot=g)
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_cis_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('p-value distribution') + ggtitle(paste('Histogram of', tissue_name, 'with', all_n_tests_cis_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_cis_permute.png', sep=''), plot=g)
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('p-value distribution') + ggtitle(paste('Histogram of', tissue_name, 'with', all_n_tests_trans_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_trans_permute.png', sep=''), plot=g)
g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_cis_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_cis_total, 'tests', sep=' '))
ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_cis.png', sep=''), plot=g)
g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_trans.png', sep=''), plot=g)
g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_cis_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_cis_total, 'tests', sep=' '))
ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_cis_permute.png', sep=''), plot=g)
g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_hist_trans_permute.png', sep=''), plot=g)

# Save necessary files and write out summary
# As for actual eQTL list, write out only for pvalue<save_threshold, as requested by Yuan

save_threshold = 1e-5
null_df = data.frame(snps = character(), gene = character(), statistic = numeric(), pvalue = numeric(), FDR = numeric(), beta = numeric(), snp_chr = character(), snp_pos = numeric(), gene_chr = character(), gene_start = numeric(), gene_end = numeric(), mappability = numeric(), tissue_MAF = numeric())

if (sum(cumulative_me_trans_total$pvalue<save_threshold) > 0) {
    write.table(cumulative_me_trans_total[c(1:sum(cumulative_me_trans_total$pvalue<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, "_trans_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
} else {
    write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, "_trans_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
}
if (sum(cumulative_me_cis_total$pvalue<save_threshold) > 0) {
    write.table(cumulative_me_cis_total[c(1:sum(cumulative_me_cis_total$pvalue<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, "_cis_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
} else {
    write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, "_cis_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
}
if (sum(cumulative_me_trans_permute_total$pvalue<save_threshold) > 0) {
    write.table(cumulative_me_trans_permute_total[c(1:sum(cumulative_me_trans_permute_total$pvalue<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, "_trans_permute_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
} else {
    write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, "_trans_permute_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
}
if (sum(cumulative_me_cis_permute_total$pvalue<save_threshold) > 0) {
    write.table(cumulative_me_cis_permute_total[c(1:sum(cumulative_me_cis_permute_total$pvalue<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, "_cis_permute_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
} else {
    write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, "_cis_permute_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
}

save(all_n_tests_cis_total, all_n_tests_trans_total, all_hist_vals_cis_total, all_hist_vals_trans_total, all_n_tests_cis_permute_total, all_n_tests_trans_permute_total, all_hist_vals_cis_permute_total, all_hist_vals_trans_permute_total, cumulative_me_trans_total, cumulative_me_cis_total, cumulative_me_trans_permute_total, cumulative_me_cis_permute_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_run_total_data.RData", sep=""))

num_eqtls = sum(cumulative_me_trans_total$FDR<=0.1)
if (num_eqtls > 0) {
  num_trans_genes = length(unique(cumulative_me_trans_total$gene[1:num_eqtls]))
} else {
  num_trans_genes = 0
}

summary_file = paste(write_dir, 'run_summary/', tissue_name, '_run_summary.txt', sep='')
fileConn<-file(summary_file)
line1 = paste('MatrixEQTL run summary for:', tissue_name, sep='\t')
line2 = paste('Number of genes tested:', gene_set_size, sep='\t')
line3 = paste('Number of SNPs before MAF filter:', snp_set_size_init_total, sep='\t')
line4 = paste('Number of SNPs tested:', snp_set_size_final_total, sep='\t')
line5 = paste('Number of trans (off-chromosome) tests:', all_n_tests_trans_total, sep='\t')
line6 = paste('Number of trans tests for permute genotypes:', round(mean(all_n_tests_trans_permute_total)), sep='\t')
line7 = paste('Total number of trans-eQTLs:', num_eqtls, sep='\t')
line8 = paste('Total number of unique trans-Genes before cross-mapping filter:', num_trans_genes, sep='\t')
writeLines(c(line1, line2, line3, line4, line5, line6, line7, line8), fileConn)
close(fileConn)