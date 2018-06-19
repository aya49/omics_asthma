##################################################
#  MatrixEQTL_FDR_control_subset_runs_final.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/MatrixEQTL_FDR_control_subset_runs_final.R
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
# original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/subset_runs/'
# write_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/subset_runs/'

original_dir = args[1]
# This is really fast - just iterate through the tissues
#tissue_name = args[2]
write_dir = args[2]

# library(ggplot2)
# library(dplyr)
dir.create(file.path(write_dir), showWarnings = FALSE)
dir.create(file.path(write_dir, 'eqtl_list'), showWarnings = FALSE)

tissues = list.files(path=paste(original_dir))

for (tissue_name in tissues) {
    print(tissue_name)
    list_files = list.files(path=paste(original_dir, tissue_name, '/', sep=''))
    print(length(list_files))

    # Required data to aggregate - 3 subsetting runs
    cis_subset_cumulative_me_trans_total = vector('list', length(list_files))
    cis_subset_all_hist_vals_trans_total = rep(0,100)
    cis_subset_all_hist_vals_cis_total = rep(0,100)
    cis_subset_all_n_tests_trans_total = 0
    cis_subset_snp_set_size_final_total = 0

    gwas_subset_cumulative_me_trans_total = vector('list', length(list_files))
    gwas_subset_all_hist_vals_trans_total = rep(0,100)
    gwas_subset_all_hist_vals_cis_total = rep(0,100)
    gwas_subset_all_n_tests_trans_total = 0
    gwas_subset_snp_set_size_final_total = 0

    ld_subset_cumulative_me_trans_total = vector('list', length(list_files))
    ld_subset_all_hist_vals_trans_total = rep(0,100)
    ld_subset_all_hist_vals_cis_total = rep(0,100)
    ld_subset_all_n_tests_trans_total = 0
    ld_subset_snp_set_size_final_total = 0

    i = 1
    for (item in list_files) {
        load(paste(original_dir, tissue_name, '/', item, sep=""))
        
        cis_subset_cumulative_me_trans_total[[i]] = me_cis$trans$eqtls
        cis_subset_all_hist_vals_trans_total = cis_subset_all_hist_vals_trans_total + me_cis$trans$hist.counts
        cis_subset_all_hist_vals_cis_total = cis_subset_all_hist_vals_cis_total + me_cis$cis$hist.counts
        cis_subset_all_n_tests_trans_total = cis_subset_all_n_tests_trans_total + me_cis$trans$ntests
        cis_subset_snp_set_size_final_total = cis_subset_snp_set_size_final_total + cis_subset_size

        gwas_subset_cumulative_me_trans_total[[i]] = me_gwas$trans$eqtls
        gwas_subset_all_hist_vals_trans_total = gwas_subset_all_hist_vals_trans_total + me_gwas$trans$hist.counts
        gwas_subset_all_hist_vals_cis_total = gwas_subset_all_hist_vals_cis_total + me_gwas$cis$hist.counts
        gwas_subset_all_n_tests_trans_total = gwas_subset_all_n_tests_trans_total + me_gwas$trans$ntests
        gwas_subset_snp_set_size_final_total = gwas_subset_snp_set_size_final_total + gwas_subset_size

        ld_subset_cumulative_me_trans_total[[i]] = me_ld$trans$eqtls
        ld_subset_all_hist_vals_trans_total = ld_subset_all_hist_vals_trans_total + me_ld$trans$hist.counts
        ld_subset_all_hist_vals_cis_total = ld_subset_all_hist_vals_cis_total + me_ld$cis$hist.counts
        ld_subset_all_n_tests_trans_total = ld_subset_all_n_tests_trans_total + me_ld$trans$ntests
        ld_subset_snp_set_size_final_total = ld_subset_snp_set_size_final_total + ld_subset_size

        i = i+1
    }

    cis_subset_cumulative_me_trans_total = do.call("rbind", cis_subset_cumulative_me_trans_total)
    gwas_subset_cumulative_me_trans_total = do.call("rbind", gwas_subset_cumulative_me_trans_total)
    ld_subset_cumulative_me_trans_total = do.call("rbind", ld_subset_cumulative_me_trans_total)

    # cumulative_me_trans_permute_total = do.call("rbind", cumulative_me_trans_permute_total)
    print(dim(cis_subset_cumulative_me_trans_total))
    print(dim(gwas_subset_cumulative_me_trans_total))
    print(dim(ld_subset_cumulative_me_trans_total))

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

    cis_subset_cumulative_me_trans_total = cis_subset_cumulative_me_trans_total[order(cis_subset_cumulative_me_trans_total$pvalue),]
    cis_subset_cumulative_me_trans_total$FDR = BH_FDR(cis_subset_cumulative_me_trans_total$pvalue, cis_subset_all_n_tests_trans_total)

    gwas_subset_cumulative_me_trans_total = gwas_subset_cumulative_me_trans_total[order(gwas_subset_cumulative_me_trans_total$pvalue),]
    gwas_subset_cumulative_me_trans_total$FDR = BH_FDR(gwas_subset_cumulative_me_trans_total$pvalue, gwas_subset_all_n_tests_trans_total)

    ld_subset_cumulative_me_trans_total = ld_subset_cumulative_me_trans_total[order(ld_subset_cumulative_me_trans_total$pvalue),]
    ld_subset_cumulative_me_trans_total$FDR = BH_FDR(ld_subset_cumulative_me_trans_total$pvalue, ld_subset_all_n_tests_trans_total)

    ## Generate p-value histograms
    # print('Generateing p-value histograms:')
    # g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
    # ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans.png', sep=''), plot=g)
    # g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
    # ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans_permute.png', sep=''), plot=g)

    # Save necessary files and write out summary
    # As for actual eQTL list, write out only for pvalue<save_threshold, as requested by Yuan

    save(cis_subset_cumulative_me_trans_total, cis_subset_all_hist_vals_trans_total, cis_subset_all_hist_vals_cis_total, cis_subset_all_n_tests_trans_total, cis_subset_snp_set_size_final_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_cis_subset_list_by_pvalue.RData", sep=""))
    save(gwas_subset_cumulative_me_trans_total, gwas_subset_all_hist_vals_trans_total, gwas_subset_all_hist_vals_cis_total, gwas_subset_all_n_tests_trans_total, gwas_subset_snp_set_size_final_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_gwas_subset_list_by_pvalue.RData", sep=""))
    save(ld_subset_cumulative_me_trans_total, ld_subset_all_hist_vals_trans_total, ld_subset_all_hist_vals_cis_total, ld_subset_all_n_tests_trans_total, ld_subset_snp_set_size_final_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_ld_subset_list_by_pvalue.RData", sep=""))

}