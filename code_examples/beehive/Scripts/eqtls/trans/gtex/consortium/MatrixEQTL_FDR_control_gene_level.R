##################################################
#  MatrixEQTL_FDR_control_gene_level.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/MatrixEQTL_FDR_control_gene_level.R
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

library(VGAM); library(fitdistrplus); library(qvalue); library(ggplot2)

# FDR Control of p-values from the null distribution
original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/gene_level_FDR/'
# tissue_name = 'muscleskeletal'
write_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/gene_level_FDR/'
plot_dir = '/tigress/BEE/RNAseq/Analysis/Figures/manuscript_figures/gene_level_FDR/plots/'

# original_dir = args[1]
# write_dir = args[2]
# plot_dir = args[3]]

library(ggplot2)
# library(dplyr)
# dir.create(file.path(write_dir), showWarnings = FALSE)
# dir.create(file.path(write_dir, 'eqtl_list'), showWarnings = FALSE)

tissue_table = read.table('/tigress/BEE/RNAseq/Data/Resources/gtex/tables/tissue_table.txt', sep='\t')

tissues = list.files(path=paste(original_dir))
main_trans_list = read.table('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/main_consortium_trans/trans-eQTLs_FDR-0.1.txt', header=T, stringsAsFactors=F) 
result_list = vector('list')

for (tissue_name in tissues) {
    if (tissue_name == 'genes_represented') {
        next
    }

    print(tissue_name)
    list_files = list.files(path=paste0(original_dir, tissue_name, '/'))
    print(length(list_files))

    # Read in the current extreme values by gene:
    test_stats = read.table(paste0(original_dir, 'genes_represented/', tissue_name, '_genes_present.txt'), header=T, stringsAsFactors=F)

    me_permute_joint = vector('list')
    # me_joint = vector('list')
    i = 1
    for (file in list_files) {
        load(paste0(original_dir, tissue_name, '/', file))
        me_permute_joint[[i]] = me_permute$trans$eqtls
        # me_joint[[i]] = me$trans$eqtls
        i = i + 1
    }

    me_permute_joint = do.call('rbind', me_permute_joint)
    # me_joint = do.call('rbind', me_joint)
    me_permute_joint = me_permute_joint[order(me_permute_joint$pvalue),]
    me_permute_joint = me_permute_joint[!duplicated(me_permute_joint$gene),]

    # if (nrow(me_joint) > 0) {
    #     me_permute = me_permute[order(me_permute$pvalue),]
    #     me_permute = me_permute[!duplicated(me_permute$gene),]
    # }

    gumbel_fit = fitdist(-log10(me_permute_joint$pvalue), 'gumbel', start=list(location = 0, scale = 1))
    mu = as.list(gumbel_fit$estimate)$location
    beta = as.list(gumbel_fit$estimate)$scale

    p_bonhack = sapply(test_stats$pvalue, function(x) {min(x*1e6, 1)})
    q_bonhack = qvalue(p_bonhack)
    p_cdf = sapply(-log10(test_stats$pvalue), function(x) {1 - pgumbel(x, location = mu, scale = beta)})
    q_cdf = qvalue(p_cdf)

    gumbel_emp_quantile_null = qgumbel(1-(c(1:nrow(test_stats)) / (nrow(test_stats) + 1)), location = mu, scale = beta)

    gumbel_fit = fitdist(-log10(test_stats$pvalue), 'gumbel', start=list(location = 0, scale = 1))
    mu = as.list(gumbel_fit$estimate)$location
    beta = as.list(gumbel_fit$estimate)$scale

    gumbel_emp_quantile_test = qgumbel(1-(c(1:nrow(test_stats)) / (nrow(test_stats) + 1)), location = mu, scale = beta)

    test_stats$pval_Uniform_Bonferroni = p_bonhack
    test_stats$qval_Uniform_Bonferroni = q_bonhack$qvalue
    test_stats$pval_gene_FDR_cdf = p_cdf
    test_stats$qval_gene_FDR_cdf = q_cdf$qvalue
    test_stats$gumbel_emp_quantile_null = gumbel_emp_quantile_null
    test_stats$gumbel_emp_quantile_test = gumbel_emp_quantile_test
    test_stats$null_pvals = me_permute_joint$pvalue[c(1:nrow(test_stats))]

    test_stats$FDR_gumbel = sapply(c(1:nrow(test_stats)), function(x) {x / (sum(-log10(test_stats$pvalue) >= test_stats$gumbel_emp_quantile_null[x]))})
    test_stats$FDR_pval = sapply(c(1:nrow(test_stats)), function(x) {x / (sum(-log10(test_stats$pvalue) >= -log10(test_stats$null_pvals[x])))})

    g = qplot(test_stats$gumbel_emp_quantile_null, -log10(test_stats$null_pvals)) + geom_abline() + xlab('Quantiles of Gumbel fit on null') + ylab('-log10(null p-value)')
    ggsave(g, file = paste0(plot_dir, tissue_name, '_qqplot1.png'))
    g = qplot(test_stats$gumbel_emp_quantile_null, -log10(test_stats$pvalue)) + geom_abline() + xlab('Quantiles of Gumbel fit on null') + ylab('-log10(test p-value)')
    ggsave(g, file = paste0(plot_dir, tissue_name, '_qqplot1.png'))
    g = qplot(test_stats$gumbel_emp_quantile_test, -log10(test_stats$pvalue)) + geom_abline() + xlab('Quantiles of Gumbel fit on test') + ylab('-log10(test p-value)')
    ggsave(g, file = paste0(plot_dir, tissue_name, '_qqplot1.png'))

    write.table(test_stats, file = paste0(write_dir, '/tissues/', tissue_name, '.txt'), col.names=T, row.names=F, quote=F, sep='\t')

    ind = which(rownames(tissue_table) == tissue_name)
    result_array = c(tissue_table$num_samples_with_geno[ind], length(unique(main_trans_list[(main_trans_list$tissue == tissue_name),]$gene)), sum(test_stats$qval_Uniform_Bonferroni < 0.1), sum(test_stats$qval_gene_FDR_cdf < 0.1), sum(test_stats$FDR_gumbel < 0.1), sum(test_stats$FDR_pval < 0.1))
    print(result_array)
    result_list[[as.character(tissue_table$SMTSD[ind])]] = result_array

    result_df = do.call('rbind', result_list)
    colnames(result_df) = c('nSamples', 'eGenes', 'eGenes_BonHack_qval', 'eGenes_Gumbel_qval', 'FDR_gumbel', 'FDR_pval')
    write.table(result_df, file = paste0(write_dir, 'num_eGenes.txt'), row.names=T, quote=F, sep='\t')

}

result_df = do.call('rbind', result_list)
colnames(result_df) = c('eGenes', 'eGenes_BonHack_qval', 'eGenes_Gumbel_qval', 'FDR_gumbel', 'FDR_pval')
write.table(result_df, file = paste0(write_dir, 'num_eGenes.txt'), row.names=T, quote=F, sep='\t')


#     gumbel_fit = fitdist(abs(me_permute_joint$statistic), 'gumbel', start=list(location = 0, scale = 1))

#     mu = as.list(gumbel_fit$estimate)$location
#     beta = as.list(gumbel_fit$estimate)$scale

#     p_bonhack = sapply(test_stats$pvalue, function(x) {min(x*1e6, 1)})
#     p_ks = sapply(abs(test_stats$statistic), function(x) {ks.test(x, "pgumbel", mu, beta)$p.value})
#     p_cdf = sapply(abs(test_stats$statistic), function(x) {1 - pgumbel(x, location = mu, scale = beta)})
#     q_bonhack = qvalue(p_bonhack)
#     q_ks = qvalue(p_ks)
#     q_cdf = qvalue(p_cdf)

#     test_stats$pval_Uniform_Bonferroni = p_bonhack
#     test_stats$pval_gene_FDR_ks = p_ks
#     test_stats$pval_gene_FDR_cdf = p_cdf
#     test_stats$qval_Uniform_Bonferroni = q_bonhack$qvalue
#     test_stats$qval_gene_FDR_ks = q_ks$qvalue
#     test_stats$qval_gene_FDR_cdf = q_cdf$

#     # ggplot scripts

#     df = rbind(data.frame(x = abs(test_stats$statistic), type = 'test'), data.frame(x = abs(me_permute_joint$statistic), type = 'null'))

#     g = ggplot(df) + geom_point(aes(x, exp(-exp(-(x - mu)/beta)), colour = type)) + 
#         stat_function(fun = pgumbel, args = list(location = mu, scale = beta), colour = "black") + 
#         theme_classic() + xlab('-log10(p-value)') + ylab('CDF') + guides(colour=guide_legend(title="Gene type"))

#     ggsave(paste0(plot_dir, tissue_name, "_eGene_calibration_CDF.pdf"), plot = g, width = 5.5, height = 5)
#     ggsave(paste0(plot_dir, tissue_name, "_eGene_calibration_CDF.png"), plot = g, width = 5.5, height = 5)

#     # p_bonhack = sapply(me_permute_joint$pvalue, function(x) {min(x*1e6, 1)})
#     # p_ks = sapply(abs(test_stats$statistic), function(x) {ks.test(x, "pgumbel", mu, beta)$p.value})

#     # hist(qnorm(1 - p_bonhack))
#     # hist(qnorm(1 - p_ks))

#     # result_list[[tissue_name]] = c(sum(main_trans_list$tissue == tissue_name), sum(test_stats$pvalue < 5e-11), sum(q$qvalue < 0.10))
#     # result_array = c(length(unique(main_trans_list[which(main_trans_list$tissue == tissue_name),]$gene)), sum(test_stats$pvalue < 5e-11), sum(q$qvalue < 0.10))
#     # print(result_array)

# # }

# # result_df = data.frame(tissue_name = result_array)
# # rownames(result_df) = c('BH on SNP-gene Pairs', 'Bon Hack', 'Gumbel fit')
# write.table(test_stats, file = paste0(write_dir, tissue_name, '.txt'), col.names=T, row.names=F, quote=F, sep='\t')
