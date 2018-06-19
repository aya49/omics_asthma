##################################################
#  GGR_validation_collect_stats_redux.R
#
#  $proj/Scripts/causality/GGR/gtex/GGR_validation_collect_stats_redux.R
# 
#  GGR validation pipeline with trans-eQTLs and MR implementations.
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# FDR Control of p-values from the null distribution
condition = '0mean-unnormalized_g-null_l-fdr'
original_dir = paste0('/tigress/BEE/RNAseq/Output/causality/GGR/output/all_by_all/', condition, '/')
regression = 'enet-1'
write_dir = paste0('/tigress/BEE/RNAseq/Output/causality/GGR/final_list/all_by_all/', condition, '/')
graph_dir = paste0('/tigress/BEE/RNAseq/Analysis/Figures/GGR_validation_plots/all_by_all/', condition, '/')

condition = args[4]
regression = args[5]
original_dir = paste0(args[1], condition, '/')
write_dir = paste0(args[2], condition, '/')
graph_dir = paste0(args[3], condition, '/')
dir.create(write_dir, showWarnings = FALSE)
dir.create(graph_dir, showWarnings = FALSE)

# name of files in the directory
file_identifier = paste0('prot2TPM-er-reps_', condition, '-0.05_', regression, '-union-network_part_')

library(ggplot2)
library(dplyr)
library(qvalue)

# Function for Benjamini-Hochberg FDR correction
BH_FDR = function(pvalues, ntests) {
  alpha_hat = sapply(c(1:length(pvalues)), function(x) {pvalues[x] * ntests / x})
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
  }
  return(alpha_hat_copy)
}

# Plot p-value histogram
plot_pval_hist = function(all_hist_vals_trans_total, all_hist_vals_trans_permute_total) {
	pval_hist_df = data.frame(ind = c(1:50)/50)
	pval_hist_df$counts_true = sapply(c(1:50), function(x) {all_hist_vals_trans_total[2*x -1] + all_hist_vals_trans_total[2*x]})
	pval_hist_df$counts_rand = sapply(c(1:50), function(x) {all_hist_vals_trans_permute_total[2*x -1] + all_hist_vals_trans_permute_total[2*x]})

	g = ggplot(pval_hist_df, aes(ind)) +
		geom_bar(aes(y=counts_true, fill = "Proposed gene pairs"), col = I("black"), size = 0.1, stat="identity") +
		geom_bar(aes(y=counts_rand, fill = "Permuted tests"), col = I("black"), size = 0.1, stat="identity") +
		theme_classic() + 
		xlab("GTEx trans association p-values") + 
		ylab('Frequency') + 
		scale_x_continuous(expand = c(0.005, 0)) +
		#scale_y_continuous(label= scientific_format(), expand = c(0.005, 0)) +
		scale_y_continuous(expand = c(0.005, 0)) +
		scale_fill_manual(name = '', values = c("Proposed gene pairs" = "#0000FF80", "Permuted tests" = "#FF000080"), guide = guide_legend(reverse = F)) + 
		guides(color = guide_legend(nrow = 2, ncol = 1)) + 
		theme(legend.justification=c(0.9,0.6), legend.position=c(0.97,0.97))

	return(g)
}

# Plot qq plots - both for a single tissue and for multiple tissues
plot_qq_plot = function(pval_df) {

	if (length(unique(pval_df$tissue)) == 1) {
		g = ggplot(pval_df, aes(-log10(null), -log10(true), colour = 'black')) + 
			geom_point() + 
			geom_abline(slope = 1, intercept = 0) + 
			theme_classic() + 
			ylab('-log(p-value) for proposed gene pairs') + 
			xlab('-log(p-value) for permuted tests') + 
			guides(colour = FALSE)
	} else {
		pval_df$tissue[pval_df$tissue == 'lung'] = 'Lung'

		g = ggplot(pval_df, aes(-log10(null), -log10(true), colour = factor(tissue))) + 
			geom_point() + 
			geom_abline(slope = 1, intercept = 0) + 
			theme_classic() + 
			labs(tissue="Tissue") +
			ylab('-log(p-value) for proposed gene pairs') + 
			xlab('-log(p-value) for permuted tests') + 
			guides(colour = guide_legend(title = "Tissue")) +
			theme(legend.justification=c(0.9,0.1), legend.position=c(0.9,0.1))
	}
	return(g)
}

setwd(original_dir)
list_files = list.files(path=original_dir, pattern = regression)
print(length(list_files))

# We test on five tissues
tissues = c('lung', 'adiposesubcutaneous', 'cellstransformedfibroblasts', 'arterytibial', 'thyroid')
num_parts = max(as.numeric(unique(sapply(list_files, function(x) {strsplit(x, '_')[[1]][7]}))))

paste0('condition: ', condition, ', regression: ', regression)
# Lists that will collect various statistics for analysis
trans_pi0 = list()
gene_MR_pi0 = list()
n_below_0.1 = list()
n_below_0.2 = list()
edges_below_0.1 = list()
edges_below_0.2 = list()
n_edges_below_0.1 = list()
n_edges_below_0.2 = list()
ntests = list()
pi_zeros = list()
qq_plot_collect = list()

# Iterate through tissues
for (tissue in tissues) {
	print(tissue)

	all_hist_vals_trans_total = rep(0,100)
	all_hist_vals_trans_permute_total = rep(0,100)
	all_n_tests_trans_total = 0
	all_n_tests_trans_permute_total = 0

	full_master_table = vector('list', num_parts)
	full_master_table_perm = vector('list', num_parts)

	# Iterate through files and collect statistics
	for (i in c(1:num_parts)) {
		file = paste0(file_identifier, i, '_', tissue, '.RData')
		if (file.exists(file)) {
			load(file)
			all_hist_vals_trans_total = all_hist_vals_trans_total + cumulative_me_total$trans$hist.counts
			all_hist_vals_trans_permute_total = all_hist_vals_trans_permute_total + cumulative_me_permute_total$trans$hist.counts

			all_n_tests_trans_total = all_n_tests_trans_total + cumulative_me_total$trans$ntests
			all_n_tests_trans_permute_total = all_n_tests_trans_permute_total + cumulative_me_permute_total$trans$ntests
			full_master_table[[i]] = cumulative_me_total$trans$eqtls
			full_master_table_perm[[i]] = cumulative_me_permute_total$trans$eqtls
		} else {
			print(paste('File missing:', condition, regression, tissue, i, sep = ' '))
		}
	}

	# Uncomment if pval histogram is wanted:
	# g = plot_pval_hist(all_hist_vals_trans_total, all_hist_vals_trans_permute_total)
	# ggsave(filename = paste0(graph_dir, regression, '_', tissue, '_pvalue_hist.pdf'), plot = g)

	full_master_table = do.call("rbind", full_master_table)
	full_master_table_perm = do.call("rbind", full_master_table_perm)
	full_master_table = full_master_table[order(full_master_table$pvalue),]
	full_master_table_perm = full_master_table_perm[order(full_master_table_perm$pvalue),]

	# Uncomment if qq plot is wanted
	# temp_sum = sum(full_master_table$pvalue < 1e-4)
	# qq_plot_df = data.frame(true = full_master_table$pvalue[c(1:temp_sum,seq(from=temp_sum+1, to=nrow(full_master_table), by = 30))], null = full_master_table_perm$pvalue[c(1:temp_sum,seq(from=temp_sum+1, to=nrow(full_master_table), by = 30))], tissue = tissue, stringsAsFactors=FALSE)
	# qq_plot_collect[[tissue]] = qq_plot_df
	# g = plot_qq_plot(qq_plot_df)
	# ggsave(filename = paste0(graph_dir, regression, '_', tissue, '_qq_plot.pdf'), plot = g)

	qobj = qvalue(full_master_table$pvalue)
	# print(qobj$pi0)
	full_master_table$FDR = qobj$qvalues

	qobj_perm = qvalue(full_master_table_perm$pvalue)
	# print(qobj_perm$pi0)
	full_master_table_perm$FDR = qobj$qvalues

	n_below_0.1[tissue] = sum(full_master_table$FDR<0.1)
	n_below_0.2[tissue] = sum(full_master_table$FDR<0.2)
	edges_below_0.1[[tissue]] = paste(full_master_table$cis_gene, full_master_table$gene, sep="-")[full_master_table$FDR<0.1]
	edges_below_0.2[[tissue]] = paste(full_master_table$cis_gene, full_master_table$gene, sep="-")[full_master_table$FDR<0.2]
	n_edges_below_0.1[tissue] = length(unique(as.character(edges_below_0.1[[tissue]])))
	n_edges_below_0.2[tissue] = length(unique(as.character(edges_below_0.2[[tissue]])))

	ntests[tissue] = all_n_tests_trans_total
	pi_zeros[tissue] = qobj$pi0

	# save the final table
	# write.table(full_master_table[full_master_table$pvalue < 0.01,], file=paste0(write_dir, tissue, "_", regression, "_trans_eQTL_table.txt"), sep='\t', quote=F, row.names=F)
}

# Uncomment if qq plot is wanted
# qq_plot_collected = do.call('rbind', qq_plot_collect)
# qq_plot_collected$tissue[qq_plot_collected$tissue == 'lung'] = 'Lung'
# qq_plot_collected$tissue[qq_plot_collected$tissue == 'adiposesubcutaneous'] = 'Adipose Subcutaneous'
# qq_plot_collected$tissue[qq_plot_collected$tissue == 'cellstransformedfibroblasts'] = 'Transformed Fibroblasts'
# qq_plot_collected$tissue[qq_plot_collected$tissue == 'arterytibial'] = 'Artery Tibial'
# qq_plot_collected$tissue[qq_plot_collected$tissue == 'thyroid'] = 'Thyroid'

# g = plot_qq_plot(qq_plot_collected)
# ggsave(filename = paste0(graph_dir, regression, '_all_tissues_qq_plot.pdf'), plot = g)

all_tissues_list_0.1 = as.character(unlist(edges_below_0.1))
all_tissues_list_0.2 = as.character(unlist(edges_below_0.2))
other_tissues_list_0.1 = all_tissues_list_0.1[c((n_below_0.1[['lung']]+1):length(all_tissues_list_0.1))]
other_tissues_list_0.2 = all_tissues_list_0.2[c((n_below_0.2[['lung']]+1):length(all_tissues_list_0.2))]

# n_specific_edges_0.1 = sum(sapply(unique(edges_below_0.1[['lung']]), function(x) {!(x %in% other_tissues_list_0.1)}))
# n_specific_edges_0.2 = sum(sapply(unique(edges_below_0.2[['lung']]), function(x) {!(x %in% other_tissues_list_0.2)}))
# print(n_specific_edges_0.1)
# print(n_specific_edges_0.2)

summary_table = cbind(do.call('rbind', n_below_0.1), do.call('rbind', n_below_0.2))
summary_table = cbind(summary_table, do.call('rbind', n_edges_below_0.1))
summary_table = cbind(summary_table, do.call('rbind', n_edges_below_0.2))
summary_table = cbind(summary_table, do.call('rbind', ntests))
summary_table = cbind(summary_table, do.call('rbind', pi_zeros))
colnames(summary_table) = c('n_below_0.1', 'n_below_0.2', 'n_edges_below_0.1', 'n_edges_below_0.2', 'ntests', 'pi_zeros')
print(summary_table)

total_edges = length(unique(paste0(full_master_table$cis_gene, full_master_table$gene)))
print(total_edges)

a = summary_table['lung','n_edges_below_0.2']
b = as.numeric(colSums(summary_table)['n_edges_below_0.2'] - a)
c = total_edges
d = total_edges*4

print(paste0("Fisher's exact test p-value for lung enrichment: ", fisher.test(matrix( c(a,c,b,d), nrow=2 ))$p.value))
