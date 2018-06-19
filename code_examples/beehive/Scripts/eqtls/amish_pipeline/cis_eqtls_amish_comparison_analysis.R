##################################################
#  cis_eqtls_amish_comparison_analysis.R
#
#  $proj/Scripts/eqtls/amish_pipeline/cis_eqtls_amish_comparison_analysis.R
#
#  Script for generating comparison stats for Amish and GTEx cis-eQTLs
#
#  Author: Brian Jo
#
##################################################

save_dir = '/tigress/BEE/RNAseq/Output/cis-mapping/amish/cellsebvtransformedlymphocytes/'
fig_dir = '/tigress/BEE/RNAseq/Analysis/Figures/amish_comparison/'

if (!file.exists(paste0(fig_dir, 'comp_table_all.RData'))) {

	comp_table_all = read.table(paste0(save_dir, 'gtex_amish_comparison_ciseqtls_chr1.txt'), stringsAsFactors=F, sep='\t', header=T)
	for (i in c(2:22)) {
		print(i)
		comp_table = read.table(paste0(save_dir, 'gtex_amish_comparison_ciseqtls_chr', i, '.txt'), stringsAsFactors=F, sep='\t', header=T)
		comp_table_all = rbind(comp_table_all, comp_table)
	}
	comp_table_all = comp_table_all[!is.na(comp_table_all$p_wald),]
	# dim(comp_table_all)
	# [1] 10622728       26

	# For amish data, the categories are ~5%, ~10%, and rest
	comp_table_all$amish_af_cat = 'amish_common'
	comp_table_all$amish_af_cat[comp_table_all$amish_af <= 0.05] = 'amish_rare'

	# For GTEx data, the categories are ~
	comp_table_all$gtex_af_cat = 'gtex_common'
	comp_table_all$gtex_af_cat[comp_table_all$MAF <= 0.05] = 'gtex_rare'
	comp_table_all$gtex_af_cat[comp_table_all$MAF <= 0.01] = 'gtex_very_rare'

	# About 18M entries total
	# save(comp_table_all, file = paste0(fig_dir, 'comp_table_all.RData'))

	# Calculate FDR approximations for both datasets
	comp_table_all$FDR = p.adjust(comp_table_all$pvalue)
	comp_table_all$amish_FDR = p.adjust(comp_table_all$p_wald)

	vec1 = comp_table_all$FDR <= 0.05
	vec2 = comp_table_all$amish_FDR <= 0.05
	vec = vec1 + vec2
	subset = comp_table_all[(vec>0),]

	write.table(subset, quote = F, sep = '\t', row.names = F, paste0(save_dir, 'gtex_amish_sig_ciseqtls.txt'))
	save(comp_table_all, subset, file = paste0(fig_dir, 'comp_table_all.RData'))

} else {
	load(paste0(fig_dir, 'comp_table_all.RData'))
}

# for lymphocytes
# > length(unique(comp_table_all$gene))
# [1] 12108
# > length(unique(comp_table_all[vec1,]$gene))
# [1] 392
# > length(unique(comp_table_all[vec2,]$gene))
# [1] 74
# > length(unique(comp_table_all[vec==2,]$gene))
# [1] 20

# # Make three qq plots
library(ggplot2)
library(cowplot)
qq_df = data.frame(GTEx_pval = comp_table_all$pvalue[order(comp_table_all$pvalue)], amish_pval = comp_table_all$p_wald[order(comp_table_all$p_wald)], GTEx_permuted_pval = comp_table_all$pvalue_permuted[order(comp_table_all$pvalue_permuted)])
qq_df = qq_df[1:50000,]

g1 = ggplot(qq_df, aes(-log10(amish_pval), -log10(GTEx_pval))) + geom_point() + geom_abline(slope=1) + xlab('Amish log p-vals') + ylab('GTEx log p-vals')
g2 = ggplot(qq_df, aes(-log10(GTEx_permuted_pval), -log10(amish_pval))) + geom_point() + geom_abline(slope=1) + xlab('Permuted log p-vals') + ylab('Amish log p-vals')
g3 = ggplot(qq_df, aes(-log10(GTEx_permuted_pval), -log10(GTEx_pval))) + geom_point() + geom_abline(slope=1) + xlab('Permuted log p-vals') + ylab('GTEx log p-vals')
png(paste0(fig_dir, 'qq_plots_', 'cellsebvtransformedlymphocytes', '.png'), width = 1500, height = 500)
plot_grid(g1, g2, g3, labels = c('A', 'B', 'C'), ncol = 3)
dev.off()

# # Make gene-specfic data frame
# gene_list = unique(comp_table_all$gene)

# min_pvals = sapply(gene_list, function(x) {min(comp_table_all[which(comp_table_all$gene == x),]$pvalue)})
# min_amish_pvals = sapply(gene_list, function(x) {min(comp_table_all[which(comp_table_all$gene == x),]$p_wald)})
# min_perm_pvals = sapply(gene_list, function(x) {min(comp_table_all[which(comp_table_all$gene == x),]$pvalue_permuted)})

# qq_df_gene = data.frame(GTEx_pval = min_pvals[order(min_pvals)], amish_pval = min_amish_pvals[order(min_amish_pvals)], GTEx_permuted_pval = min_perm_pvals[order(min_perm_pvals)])
# g1 = ggplot(qq_df_gene, aes(-log10(amish_pval), -log10(GTEx_pval))) + geom_point() + geom_abline(slope=1) + xlab('Amish log p-vals') + ylab('GTEx log p-vals')
# g2 = ggplot(qq_df_gene, aes(-log10(GTEx_permuted_pval), -log10(amish_pval))) + geom_point() + geom_abline(slope=1) + xlab('Permuted log p-vals') + ylab('Amish log p-vals')
# g3 = ggplot(qq_df_gene, aes(-log10(GTEx_permuted_pval), -log10(GTEx_pval))) + geom_point() + geom_abline(slope=1) + xlab('Permuted log p-vals') + ylab('GTEx log p-vals')
# png(paste0(fig_dir, 'qq_plots_gene_', 'cellsebvtransformedlymphocytes', '.png'), width = 1500, height = 500)
# plot_grid(g1, g2, g3, labels = c('A', 'B', 'C'), ncol = 3)
# dev.off()

# Compare effect sizes - make sure to appropriately flip the effect size signs
# z = as.character(sapply(comp_table_all$snps, function(x) {strsplit(x, '_')[[1]][3]}))
# inds = (z == comp_table_all$amish_allele1)
# comp_table_all$amish_beta_cor = comp_table_all$amish_beta
# comp_table_all$amish_beta_cor[inds] = -comp_table_all$amish_beta_cor[inds]

# g = ggplot(comp_table_all, aes(statistic, (amish_beta_cor/amish_se))) + geom_point() + facet_grid(gtex_af_cat ~ amish_af_cat) + xlab('GTEx') + ylab('amish') + geom_abline(slope = 1)
# ggsave(g, filename = paste0(fig_dir, 'effect_size_comparison_', 'cellsebvtransformedlymphocytes', '.png'), width=8, height=8)

g = ggplot(subset, aes(statistic, amish_beta_cor/amish_se)) + geom_point() + facet_grid(gtex_af_cat ~ amish_af_cat) + xlab('GTEx') + ylab('amish') + geom_abline(slope = 1)
ggsave(g, filename = paste0(fig_dir, 'effect_size_comparison_sig_', 'cellsebvtransformedlymphocytes', '.png'), width=8, height=8)

g = ggplot(subset, aes(statistic, amish_beta_cor/amish_se)) + geom_point() + ylab('amish') + geom_abline(slope = 1)
ggsave(g, filename = paste0(fig_dir, 'effect_size_comparison_sig_', 'cellsebvtransformedlymphocytes', '.png'), width=8, height=8)

# # the chi-squared test of independence for FDR 0.05
# chisq.test(as.matrix(rbind(c(sum(vec==2), sum(vec1)-sum(vec==2)), c(sum(vec2)-sum(vec==2), sum(vec==0)))))

# #         [,1] [,2]
# # [1,]   11019  613
# # [2,] 1822512 1357

# # 	Pearson's Chi-squared test with Yates' continuity
# # 	correction

# # data:  as.matrix(rbind(c(sum(vec1) - sum(vec == 2), sum(vec == 2)),     c(sum(vec == 0), sum(vec2) - sum(vec == 2))))
# # X-squared = 29053, df = 1, p-value < 2.2e-16

