##################################################
#  MatrixEQTL_FDR_control_subset_runs_match_random.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/MatrixEQTL_FDR_control_subset_runs_match_random.R
# 
#  Post-hoc analysis where we have a set of random-matched SNPs tested
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')
library(ggplot2)

# FDR Control of p-values from the null distribution
original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/subset_runs/'
plot_dir = '/tigress/BEE/RNAseq/Analysis/Figures/manuscript_figures/subset_enrichment/'

original_dir = args[1]
# This is really fast - just iterate through the tissues
#tissue_name = args[2]
plot_dir = args[2]

# library(ggplot2)
# library(dplyr)
# dir.create(file.path(write_dir), showWarnings = FALSE)
# dir.create(file.path(write_dir, 'eqtl_list'), showWarnings = FALSE)

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

# plot_pval_hist = function(pval_hist, pval_hist_random) {
#     pval_hist_df = data.frame(ind = c(1:50)/50)
#     pval_hist_df$counts_true = sapply(c(1:50), function(x) {pval_hist[2*x -1] + pval_hist[2*x]})
#     pval_hist_df$counts_rand = sapply(c(1:50), function(x) {pval_hist_random[2*x -1] + pval_hist_random[2*x]})

#     g = ggplot(pval_hist_df, aes(ind))+
#         geom_bar(aes(y=counts_true, fill = "Subset"), col = I("black"), size = 0.1, stat="identity") +
#         geom_bar(aes(y=counts_rand, fill = "Matched Random"), col = I("black"), size = 0.1, stat="identity") +
#       # geom_histogram(aes(x = chr9$pvalue, fill = "rs13234269"), 
#       #                breaks = brx1, bins = 50, col = I("black"), size = 0.1)+
#       # geom_histogram(aes(x = chrx$pvalue, fill = "matched cis-subset"), 
#       #                breaks = brx2, bins = 50, col = I("black"), size = 0.1)+
#       theme_classic() + 
#       xlab("GTEx trans association p-values") + 
#       ylab('Frequency') + 
#       scale_x_continuous(expand = c(0.005, 0))+
#       #scale_y_continuous(label= scientific_format(), expand = c(0.005, 0)) +
#       scale_y_continuous(expand = c(0.005, 0)) +
#       scale_fill_manual(name = '', values = c("Matched Random" = "#0000FF80",
#                                               "Subset" = "#FF000080"), 
#                         guide = guide_legend(reverse = F)) +
#       guides(color = guide_legend(nrow = 2, ncol = 1))+
#       theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#             text=element_text(family="Helvetica", size = 6),
#             axis.title.y = element_text(size = 6,family = "Helvetica"),
#             # axis.title.y = element_blank(),
#             axis.title.x = element_text(size = 6,family = "Helvetica"),
#             axis.text.x = element_text(size = 4,family = "Helvetica", hjust = 0.5, vjust = 0.5),
#             axis.text.y = element_text(size = 4,family = "Helvetica"),
#             # legend.direction = 'horizontal', 
#             legend.key.size = unit(2, 'line'),
#             # legend.key = element_rect(color = white),
#             axis.ticks.y = element_line(size = 0.1),
#             axis.ticks.x=element_line(size = 0.1),
#             legend.justification=c(0.9,0.6), 
#             legend.position=c(1,1), 
#             legend.text=element_text(size=6),
#             legend.background = element_blank(),
#             legend.key.height=unit(0.5,"line"), 
#             legend.key.width = unit(0.5,"line"),
#             plot.margin=unit(c(0.5,1.7,0,0),"mm"))

#     return(g)
# }

qq_plot = function(me_test, me_random) {
    rows = min(nrow(me_test), nrow(me_random))
    plot_df = data.frame(test_pvals = -log10(me_test$pvalue)[c(1:rows)], random_pvals = -log10(me_random$pvalue)[c(1:rows)], significant = (me_test$FDR < 0.1)[c(1:rows)])

    g = ggplot(plot_df, aes(random_pvals, test_pvals, shape = factor(significant))) + geom_point() + labs(shape='FDR < 0.1') + geom_abline(slope = 1, intercept = 0)
    return(g)
}

tissues = list.files(path=paste(original_dir))

z_cis_stratified = array(0, dim = c(2,2,44))
z_gwas_stratified = array(0, dim = c(2,2,44))
z_gwas_simple = array(0, dim = c(2,2))
gwas_sig = c()
gwas_rand = c()
gwas_total = c()

ind = 1
for (tissue_name in tissues) {
    print(tissue_name)
    list_files = list.files(path=paste(original_dir, tissue_name, '/', sep=''))
    print(length(list_files))

    # Required data to aggregate - 3 subsetting runs
    cis_subset_cumulative_me_trans_total = vector('list', length(list_files))
    cis_subset_all_hist_vals_trans_total = rep(0,100)
    cis_subset_all_n_tests_trans_total = 0
    cis_subset_snp_set_size_final_total = 0

    cis_random_subset_cumulative_me_trans_total = vector('list', length(list_files))
    cis_random_subset_all_hist_vals_trans_total = rep(0,100)
    cis_random_subset_all_n_tests_trans_total = 0
    cis_random_subset_snp_set_size_final_total = 0

    gwas_subset_cumulative_me_trans_total = vector('list', length(list_files))
    gwas_subset_all_hist_vals_trans_total = rep(0,100)
    gwas_subset_all_n_tests_trans_total = 0
    gwas_subset_snp_set_size_final_total = 0

    gwas_random_subset_cumulative_me_trans_total = vector('list', length(list_files))
    gwas_random_subset_all_hist_vals_trans_total = rep(0,100)
    gwas_random_subset_all_n_tests_trans_total = 0
    gwas_random_subset_snp_set_size_final_total = 0

    # ld_subset_cumulative_me_trans_total = vector('list', length(list_files))
    # ld_subset_all_hist_vals_trans_total = rep(0,100)
    # ld_subset_all_hist_vals_cis_total = rep(0,100)
    # ld_subset_all_n_tests_trans_total = 0
    # ld_subset_snp_set_size_final_total = 0

    i = 1
    for (item in list_files) {
        load(paste(original_dir, tissue_name, '/', item, sep=""))
        
        cis_subset_cumulative_me_trans_total[[i]] = me_cis$trans$eqtls
        cis_subset_all_hist_vals_trans_total = cis_subset_all_hist_vals_trans_total + me_cis$trans$hist.counts
        cis_subset_all_n_tests_trans_total = cis_subset_all_n_tests_trans_total + me_cis$trans$ntests
        cis_subset_snp_set_size_final_total = cis_subset_snp_set_size_final_total + cis_subset_size

        cis_random_subset_cumulative_me_trans_total[[i]] = me_cis_random$trans$eqtls
        cis_random_subset_all_hist_vals_trans_total = cis_subset_all_hist_vals_trans_total + me_cis_random$trans$hist.counts
        cis_random_subset_all_n_tests_trans_total = cis_subset_all_n_tests_trans_total + me_cis_random$trans$ntests
        cis_random_subset_snp_set_size_final_total = cis_subset_snp_set_size_final_total + cis_subset_size_random

        gwas_subset_cumulative_me_trans_total[[i]] = me_gwas$trans$eqtls
        gwas_subset_all_hist_vals_trans_total = gwas_subset_all_hist_vals_trans_total + me_gwas$trans$hist.counts
        gwas_subset_all_n_tests_trans_total = gwas_subset_all_n_tests_trans_total + me_gwas$trans$ntests
        gwas_subset_snp_set_size_final_total = gwas_subset_snp_set_size_final_total + gwas_subset_size

        gwas_random_subset_cumulative_me_trans_total[[i]] = me_gwas_random$trans$eqtls
        gwas_random_subset_all_hist_vals_trans_total = gwas_subset_all_hist_vals_trans_total + me_gwas_random$trans$hist.counts
        gwas_random_subset_all_n_tests_trans_total = gwas_subset_all_n_tests_trans_total + me_gwas_random$trans$ntests
        gwas_random_subset_snp_set_size_final_total = gwas_subset_snp_set_size_final_total + gwas_subset_size_random

        # ld_subset_cumulative_me_trans_total[[i]] = me_ld$trans$eqtls
        # ld_subset_all_hist_vals_trans_total = ld_subset_all_hist_vals_trans_total + me_ld$trans$hist.counts
        # ld_subset_all_hist_vals_cis_total = ld_subset_all_hist_vals_cis_total + me_ld$cis$hist.counts
        # ld_subset_all_n_tests_trans_total = ld_subset_all_n_tests_trans_total + me_ld$trans$ntests
        # ld_subset_snp_set_size_final_total = ld_subset_snp_set_size_final_total + ld_subset_size

        i = i + 1
    }

    cis_subset_cumulative_me_trans_total = do.call("rbind", cis_subset_cumulative_me_trans_total)
    cis_random_subset_cumulative_me_trans_total = do.call("rbind", cis_random_subset_cumulative_me_trans_total)
    gwas_subset_cumulative_me_trans_total = do.call("rbind", gwas_subset_cumulative_me_trans_total)
    gwas_random_subset_cumulative_me_trans_total = do.call("rbind", gwas_random_subset_cumulative_me_trans_total)

    # ld_subset_cumulative_me_trans_total = do.call("rbind", ld_subset_cumulative_me_trans_total)

    # cumulative_me_trans_permute_total = do.call("rbind", cumulative_me_trans_permute_total)
    print(dim(cis_subset_cumulative_me_trans_total))
    print(dim(cis_random_subset_cumulative_me_trans_total))
    print(dim(gwas_subset_cumulative_me_trans_total))
    print(dim(gwas_random_subset_cumulative_me_trans_total))
    # print(dim(ld_subset_cumulative_me_trans_total))

    cis_subset_cumulative_me_trans_total = cis_subset_cumulative_me_trans_total[order(cis_subset_cumulative_me_trans_total$pvalue),]
    cis_subset_cumulative_me_trans_total$FDR = BH_FDR(cis_subset_cumulative_me_trans_total$pvalue, cis_subset_all_n_tests_trans_total)

    cis_random_subset_cumulative_me_trans_total = cis_random_subset_cumulative_me_trans_total[order(cis_random_subset_cumulative_me_trans_total$pvalue),]
    cis_random_subset_cumulative_me_trans_total$FDR = BH_FDR(cis_random_subset_cumulative_me_trans_total$pvalue, cis_random_subset_all_n_tests_trans_total)

    gwas_subset_cumulative_me_trans_total = gwas_subset_cumulative_me_trans_total[order(gwas_subset_cumulative_me_trans_total$pvalue),]
    gwas_subset_cumulative_me_trans_total$FDR = BH_FDR(gwas_subset_cumulative_me_trans_total$pvalue, gwas_subset_all_n_tests_trans_total)

    gwas_random_subset_cumulative_me_trans_total = gwas_random_subset_cumulative_me_trans_total[order(gwas_random_subset_cumulative_me_trans_total$pvalue),]
    gwas_random_subset_cumulative_me_trans_total$FDR = BH_FDR(gwas_random_subset_cumulative_me_trans_total$pvalue, gwas_random_subset_all_n_tests_trans_total)

    # ld_subset_cumulative_me_trans_total = ld_subset_cumulative_me_trans_total[order(ld_subset_cumulative_me_trans_total$pvalue),]
    # ld_subset_cumulative_me_trans_total$FDR = BH_FDR(ld_subset_cumulative_me_trans_total$pvalue, ld_subset_all_n_tests_trans_total)

    # save(cis_subset_cumulative_me_trans_total, cis_random_subset_cumulative_me_trans_total, gwas_subset_cumulative_me_trans_total, gwas_random_subset_cumulative_me_trans_total, cis_subset_all_hist_vals_trans_total, cis_random_subset_all_hist_vals_trans_total, gwas_subset_all_hist_vals_trans_total, gwas_random_subset_all_hist_vals_trans_total, file = paste0(plot_dir, 'data/', tissue_name, '_subset_enrichment.RData'))
    
    z = matrix(c(sum(cis_subset_cumulative_me_trans_total$FDR < 0.5), sum(cis_random_subset_cumulative_me_trans_total$FDR < 0.5), cis_subset_snp_set_size_final_total, cis_random_subset_snp_set_size_final_total), nrow=2, ncol=2)
    z[1,1] = length(unique(cis_subset_cumulative_me_trans_total$snps[1:z[1,1]]))
    z[2,1] = length(unique(cis_random_subset_cumulative_me_trans_total$snps[1:z[2,1]]))
    z[1,2] = z[1,2] - z[1,1]
    z[2,2] = z[2,2] - z[2,1]
    print(z)
    f_cis = fisher.test(z)

    z_cis_stratified[1,1,ind] = z[1,1]
    z_cis_stratified[1,2,ind] = z[1,2]
    z_cis_stratified[2,1,ind] = z[2,1]
    z_cis_stratified[2,2,ind] = z[2,2]

    z = matrix(c(sum(gwas_subset_cumulative_me_trans_total$FDR < 0.5), sum(gwas_random_subset_cumulative_me_trans_total$FDR < 0.5), gwas_subset_snp_set_size_final_total, gwas_random_subset_snp_set_size_final_total), nrow=2, ncol=2)
    z[1,1] = length(unique(gwas_subset_cumulative_me_trans_total$snps[1:z[1,1]]))
    z[2,1] = length(unique(gwas_random_subset_cumulative_me_trans_total$snps[1:z[2,1]]))
    z[1,2] = z[1,2] - z[1,1]
    z[2,2] = z[2,2] - z[2,1]
    print(z)
    f_gwas = fisher.test(z)

    z_gwas_stratified[1,1,ind] = z[1,1]
    z_gwas_stratified[1,2,ind] = z[1,2]
    z_gwas_stratified[2,1,ind] = z[2,1]
    z_gwas_stratified[2,2,ind] = z[2,2]

    ind = ind + 1

    gwas_sig = unique(c(gwas_sig, as.character(gwas_subset_cumulative_me_trans_total$snps[1:z[1,1]])))
    gwas_rand = unique(c(gwas_rand, as.character(gwas_random_subset_cumulative_me_trans_total$snps[1:z[2,1]])))
    gwas_total = unique(c(gwas_total, as.character(gwas_subset_cumulative_me_trans_total$snps)))

    # # Generate p-value histograms
    # print('Generateing p-value histograms:')
    # g = plot_pval_hist(cis_subset_all_hist_vals_trans_total, cis_random_subset_all_hist_vals_trans_total)
    
    # # Generate qq-plots:
    # print('Generating qq-plot for:')
    # print(tissue_name)
    # g = qq_plot(cis_subset_cumulative_me_trans_total, cis_random_subset_cumulative_me_trans_total)
    # g = g + xlab('Trans -log10(p-value) for randomly matched variants') + ylab('Trans -log10(p-value) for cis-eQTLs') + ggtitle(paste0(tissue_name, ', Fisher p-value: ', f_cis$p.value))
    # ggsave(filename = paste0(plot_dir, 'plots/', tissue_name, '_cis_qq_plot.png'), plot=g)

    # g = qq_plot(gwas_subset_cumulative_me_trans_total, gwas_random_subset_cumulative_me_trans_total)
    # g = g + xlab('Trans -log10(p-value) for randomly matched variants') + ylab('Trans -log10(p-value) for trait-associated SNPs') + ggtitle(paste0(tissue_name, ', Fisher p-value: ', f_gwas$p.value))
    # ggsave(filename = paste0(plot_dir, 'plots/', tissue_name, '_gwas_qq_plot.png'), plot=g)

    # g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
    # ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans.png', sep=''), plot=g)
    # g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
    # ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans_permute.png', sep=''), plot=g)

    # Save necessary files and write out summary
    # As for actual eQTL list, write out only for pvalue<save_threshold, as requested by Yuan

    # save(cis_subset_cumulative_me_trans_total, cis_subset_all_hist_vals_trans_total, cis_subset_all_hist_vals_cis_total, cis_subset_all_n_tests_trans_total, cis_subset_snp_set_size_final_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_cis_subset_list_by_pvalue.RData", sep=""))
    # save(gwas_subset_cumulative_me_trans_total, gwas_subset_all_hist_vals_trans_total, gwas_subset_all_hist_vals_cis_total, gwas_subset_all_n_tests_trans_total, gwas_subset_snp_set_size_final_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_gwas_subset_list_by_pvalue.RData", sep=""))
    # save(ld_subset_cumulative_me_trans_total, ld_subset_all_hist_vals_trans_total, ld_subset_all_hist_vals_cis_total, ld_subset_all_n_tests_trans_total, ld_subset_snp_set_size_final_total, file=paste(write_dir, 'eqtl_list/', tissue_name, "_ld_subset_list_by_pvalue.RData", sep=""))

}

# # In the case of FDR < 0.5:
# z_cis = matrix(c(417, 120, 202871, 209019), nrow = 2, ncol = 2)
# fisher.test(z_cis)
# z_gwas = matrix(c(193, 122, 477863, 485488), nrow = 2, ncol = 2)
# fisher.test(z_gwas)

# # In the case of FDR < 0.1:
# z_cis = matrix(c(72, 50, 202871, 209019), nrow = 2, ncol = 2)
# fisher.test(z_cis)
# z_gwas = matrix(c(57, 47, 477863, 485488), nrow = 2, ncol = 2)
# fisher.test(z_gwas)

# Revised:

# One: Cochran-Mantel-Haenszel

# mantelhaen.test(z_cis_stratified, alternative = 'greater')

mantelhaen.test(z_cis_stratified)

mantelhaen.test(z_gwas_stratified)

z_gwas = matrix(c(178, 122, 10819, 10875), nrow = 2, ncol = 2)
fisher.test(z_gwas)


cis_tissues = c('adiposesubcutaneous', 'cellstransformedfibroblasts', 'muscleskeletal')
cis_tissues = c('adiposesubcutaneous', 'muscleskeletal')
gwas_tissues = c('lung', 'muscleskeletal')

tissue_name = 'adiposesubcutaneous'
load(paste0(plot_dir, 'data/', tissue_name, '_subset_enrichment.RData'))
me_test = cis_subset_cumulative_me_trans_total
me_random = cis_random_subset_cumulative_me_trans_total
rows = min(nrow(me_test), nrow(me_random))
cis_plot_df = data.frame(test_pvals = -log10(me_test$pvalue)[c(1:rows)], random_pvals = -log10(me_random$pvalue)[c(1:rows)], significant = (me_test$FDR < 0.1)[c(1:rows)], tissue = 'Subcutaneous Adipose')

tissue_name = 'lung'
load(paste0(plot_dir, 'data/', tissue_name, '_subset_enrichment.RData'))
me_test = gwas_subset_cumulative_me_trans_total
me_random = gwas_random_subset_cumulative_me_trans_total
rows = min(nrow(me_test), nrow(me_random))
gwas_plot_df = data.frame(test_pvals = -log10(me_test$pvalue)[c(1:rows)], random_pvals = -log10(me_random$pvalue)[c(1:rows)], significant = (me_test$FDR < 0.1)[c(1:rows)], tissue = 'Lung')

tissue_name = 'muscleskeletal'
load(paste0(plot_dir, 'data/', tissue_name, '_subset_enrichment.RData'))
me_test = cis_subset_cumulative_me_trans_total
me_random = cis_random_subset_cumulative_me_trans_total
rows = min(nrow(me_test), nrow(me_random))
cis_plot_df = rbind(cis_plot_df, data.frame(test_pvals = -log10(me_test$pvalue)[c(1:rows)], random_pvals = -log10(me_random$pvalue)[c(1:rows)], significant = (me_test$FDR < 0.1)[c(1:rows)], tissue = 'Skeletal Muscle'))

me_test = gwas_subset_cumulative_me_trans_total
me_random = gwas_random_subset_cumulative_me_trans_total
rows = min(nrow(me_test), nrow(me_random))
gwas_plot_df = rbind(gwas_plot_df, data.frame(test_pvals = -log10(me_test$pvalue)[c(1:rows)], random_pvals = -log10(me_random$pvalue)[c(1:rows)], significant = (me_test$FDR < 0.1)[c(1:rows)], tissue = 'Skeletal Muscle'))

g = ggplot(cis_plot_df, aes(random_pvals, test_pvals, colour = factor(tissue))) + scale_colour_manual(values = c("red","blue")) + geom_point(aes(size = significant)) + scale_size_manual(values = c(0.2, 1)) + labs(colour='tissue', size='FDR < 0.1')
g = g + geom_abline(slope = 1, intercept = 0) + xlab('Trans -log10(p-value) for all matched random variants') + ylab('Trans -log10(p-value) for cis-eQTL subset')

g = ggplot(gwas_plot_df, aes(random_pvals, test_pvals, colour = factor(tissue))) + scale_colour_manual(values = c("red","blue")) + geom_point(aes(size = significant)) + scale_size_manual(values = c(0.2, 1)) + labs(colour='tissue', size='FDR < 0.1')
g = g + geom_abline(slope = 1, intercept = 0) + xlab('Trans -log10(p-value) for all matched random variants') + ylab('Trans -log10(p-value) for cis-eQTL subset')
