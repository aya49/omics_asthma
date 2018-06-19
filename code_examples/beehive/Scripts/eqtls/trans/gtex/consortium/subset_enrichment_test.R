# MuTHER section: test basal pi0 for effect-size matched SNPs for rs13234269

# effect size = 






# test subset enrichment

tissues = rownames(read.table('/tigress/BEE/RNAseq/Data/Resources/gtex/tables/tissue_table.txt', sep='\t', stringsAsFactors=F, header=T))

all_by_all_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/eqtl_list_pval_thresh_1e-5/'
subset_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/subset_runs/eqtl_list/'

setwd('/tigress/BEE/RNAseq/Analysis/Figures/manuscript_figures/subset_enrichment/')

trans_lists = list.files(all_by_all_dir, pattern='_trans_list_by_pvalue.RData')

# Code for Wilcox test:

total_cis = data.frame(ind = numeric(), pvalue = numeric())
total_gwas = data.frame(ind = numeric(), pvalue = numeric())

sig_list = list()
for (tissue in tissues) {
	print(tissue)
	load(paste0(all_by_all_dir, trans_lists[which(sapply(trans_lists, function(x) {grepl(tissue, x)}))]))
	load(paste0(subset_dir, tissue, '_cis_subset_list_by_pvalue.RData'))
	load(paste0(subset_dir, tissue, '_gwas_subset_list_by_pvalue.RData'))
	cumulative_me_trans_total = cumulative_me_trans_total[cumulative_me_trans_total$pvalue < 1e-7,]
	cis_subset_cumulative_me_trans_total = cis_subset_cumulative_me_trans_total[cis_subset_cumulative_me_trans_total$pvalue < 1e-7,]
	gwas_subset_cumulative_me_trans_total = gwas_subset_cumulative_me_trans_total[gwas_subset_cumulative_me_trans_total$pvalue < 1e-7,]

	z = data.frame(ind = 1, pvalue = cis_subset_cumulative_me_trans_total$pvalue)
	z = rbind(z, data.frame(ind = 2, pvalue = cumulative_me_trans_total$pvalue))
	z = z[order(z$pvalue),]

	total_cis = rbind(total_cis, z)

	q = data.frame(ind = 1, pvalue = gwas_subset_cumulative_me_trans_total$pvalue)
	q = rbind(q, data.frame(ind = 2, pvalue = cumulative_me_trans_total$pvalue))
	q = q[order(q$pvalue),]

	total_gwas = rbind(total_gwas, q)

	sig_list[[tissue]] = c(wilcox.test(z$pvalue~factor(z$ind))$p.value, wilcox.test(q$pvalue~factor(q$ind))$p.value)
	# z = z[order(z$pvalue),]
	# print(wilcox.test(z$pval, z$ind))
	# print(sum(which(z$ind ==1) < nrow(z)/2))
	# print(sum(which(z$ind ==1) > nrow(z)/2))
}
sig_df = do.call('rbind', sig_list)

total_cis = total_cis[order(total_cis$pvalue)]
total_gwas = total_gwas[order(total_gwas$pvalue)]
print(wilcox.test(total_cis$pval ~ total_cis$ind))
print(wilcox.test(total_gwas$pval ~ total_gwas$ind))


# Code for chi-sq test:

count_list = list()
for (tissue in tissues) {
	print(tissue)
	load(paste0(all_by_all_dir, trans_lists[which(sapply(trans_lists, function(x) {grepl(tissue, x)}))]))
	load(paste0(subset_dir, tissue, '_cis_subset_list_by_pvalue.RData'))
	load(paste0(subset_dir, tissue, '_gwas_subset_list_by_pvalue.RData'))
	cumulative_me_trans_total = cumulative_me_trans_total[cumulative_me_trans_total$pvalue < 1e-7,]
	cis_subset_cumulative_me_trans_total = cis_subset_cumulative_me_trans_total[cis_subset_cumulative_me_trans_total$pvalue < 1e-7,]
	gwas_subset_cumulative_me_trans_total = gwas_subset_cumulative_me_trans_total[gwas_subset_cumulative_me_trans_total$pvalue < 1e-7,]

	sum_all = read.table(paste0('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/run_summary/', tissue, '_run_summary.txt'), header=F, stringsAsFactors=F, sep='\t')
	sum_cis = read.table(paste0('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/subset_runs/run_summary/', tissue, '_cis_subset_run_summary.txt'), header=F, stringsAsFactors=F, sep='\t')
	sum_gwas = read.table(paste0('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/subset_runs/run_summary/', tissue, '_gwas_subset_run_summary.txt'), header=F, stringsAsFactors=F, sep='\t')

	count_list[[tissue]] = c(as.numeric(sum_all$V2[4]), as.numeric(sum_cis$V2[4]), as.numeric(sum_gwas$V2[4]), nrow(cumulative_me_trans_total), nrow(cis_subset_cumulative_me_trans_total), nrow(gwas_subset_cumulative_me_trans_total))
}

count_df7 = data.frame(do.call('rbind', count_list))
count_df7$cis_chisq = 0
count_df7$gwas_chisq = 0

for (tissue in tissues) {
	print(tissue)
	row = count_df7[tissue,]
	count_df7[tissue, 'cis_chisq'] = chisq.test(matrix(c(row$X5, row$X4, row$X2, row$X1), nrow=2, ncol=2))$p.value
	if ((row$X5 / row$X2) < row$X4 / row$X1) {
		count_df7[tissue, 'cis_chisq'] = 1 - count_df7[tissue, 'cis_chisq']
	}
	count_df7[tissue, 'gwas_chisq'] = chisq.test(matrix(c(row$X6, row$X4, row$X3, row$X1), nrow=2, ncol=2))$p.value
	if ((row$X6 / row$X3) < row$X4 / row$X1) {
		count_df7[tissue, 'gwas_chisq'] = 1 - count_df7[tissue, 'gwas_chisq']
	}
}

chisq.test(matrix(c(sum(count_df7$X5), sum(count_df7$X4), sum(count_df7$X2), sum(count_df7$X1)), nrow=2, ncol=2))
chisq.test(matrix(c(sum(count_df7$X6), sum(count_df7$X4), sum(count_df7$X3), sum(count_df7$X1)), nrow=2, ncol=2))
