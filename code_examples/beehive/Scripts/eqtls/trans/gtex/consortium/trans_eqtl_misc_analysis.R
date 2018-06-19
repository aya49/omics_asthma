##################################################
#  trans_eqtl_misc_analysis.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/trans_eqtl_misc_analysis.R
# 
#  Generating misc stats for the manuscript
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# How many trans-eQTLs are also cis-eQTLs?
input_trans_list = read.csv(file = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/trans-eQTLs_FDR-0.1.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)

tissue_list = unique(input_trans_list$tissue)

fastqtl_file = '/tigress/BEE/RNAseq/Data/Resources/gtex/dbGaP/GTEx_phs000424/v6p.egenes.auto.mlinc.bh.RData'
if (file.exists(fastqtl_file)) {
	load(fastqtl_file)
} else {
	fastqtl_list = read.table(file = '/tigress/BEE/RNAseq/Data/Resources/gtex/dbGaP/GTEx_phs000424/v6p.egenes.auto.mlinc.bh.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	fastqtl_list$tissue_name = as.character(sapply(fastqtl_list$tissue, function(x) {tolower(gsub(" ", "", gsub("[[:punct:]]", "", x)))}))
	save(fastqtl_list, file = '/tigress/BEE/RNAseq/Data/Resources/gtex/dbGaP/GTEx_phs000424/v6p.egenes.auto.mlinc.bh.RData')
}

total_in_cis = 0
for (tissue in tissue_list) {
	print(tissue)
	input_trans_list_subset = input_trans_list[input_trans_list$tissue == tissue,]
	fastqtl_list_subset = fastqtl_list[fastqtl_list$tissue_name == tissue,]
	fastqtl_list_subset = fastqtl_list_subset[fastqtl_list_subset$fdr < 0.05,]
	print(nrow(fastqtl_list_subset))
	print(sum(sapply(input_trans_list_subset$snps, function(x) {x %in% fastqtl_list_subset$rs_id_dbSNP142_GRCh37p13})))
}


total_in_cis = 0
for (tissue in tissue_list) {
	print(tissue)
	input_trans_list_subset = input_trans_list[input_trans_list$tissue == tissue,]
	fastqtl_list_subset = fastqtl_list[fastqtl_list$tissue_name == tissue,]
	fastqtl_list_subset = fastqtl_list_subset[fastqtl_list_subset$fdr < 0.05,]
	print(nrow(fastqtl_list_subset))
	print(sum(sapply(input_trans_list_subset$snps, function(x) {x %in% fastqtl_list_subset$rs_id_dbSNP142_GRCh37p13})))
}

# Compile gene level FDR results:
gene_fdr_directory = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/gene_level_FDR/tissues/'
files = list.files(gene_fdr_directory)
trans_table = read.table('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/trans-eQTLs_FDR-0.1.txt', header=T, sep='\t')
tissue_table = read.table('/tigress/BEE/RNAseq/Data/Resources/gtex/tables/tissue_table.txt', header=T, sep='\t')

count_list = list('vector', length(files))
gene_list_1 = list('vector', length(files))
gene_list_2 = list('vector', length(files))
gene_list_3 = list('vector', length(files))
tissue_list = list('vector', length(files))

i = 1
for (tissue in files) {
	tissue_name = strsplit(tissue, '\\.')[[1]][1]
	ind = which(rownames(tissue_table) == tissue_name)
	table = read.table(paste0(gene_fdr_directory, tissue), header=T, sep='\t')
	count_list[[i]] = c(tissue_table$num_samples_with_geno[ind], length(unique(trans_table[trans_table$tissue == tissue_name,]$gene)), sum(table$qval_Uniform_Bonferroni < 0.1), sum(table$qval_gene_FDR_cdf < 0.1))
	gene_list_1[[i]] = as.character(unique(trans_table[trans_table$tissue == tissue_name,]$gene))
	gene_list_2[[i]] = as.character(table[table$qval_Uniform_Bonferroni < 0.1,]$gene)
	gene_list_3[[i]] = as.character(table[table$qval_gene_FDR_cdf < 0.1,]$gene)
	tissue_list[[i]] = as.character(tissue_table$SMTSD[ind])

	i = i + 1
}
count_table = data.frame(do.call('rbind', count_list))
colnames(count_table) = c('num_samples', 'eGenes', 'eGenes_BonHack', 'eGenes_Gumbel')

tissue_list = unlist(tissue_list)
tissue_list = tissue_list[order(-count_table$num_samples)]
count_table = count_table[order(-count_table$num_samples),]

total = c(sum(count_table$num_samples), length(unique(unlist(gene_list_1))), length(unique(unlist(gene_list_2))), length(unique(unlist(gene_list_3))))
count_table = rbind(count_table, total)
rownames(count_table) = c(tissue_list, 'Total')
write.table(count_table, file = '/tigress/BEE/RNAseq/Analysis/Figures/manuscript_figures/supplement_bundle/eGene_table.csv', sep = ',', quote=F)

# Compile intrachromosomal results:
intrachrom_results = read.table('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/intrachromosomal/trans-eQTLs_FDR-0.1.txt', header=T, sep='\t')

count_list = list('vector', length(files))
gene_list = list('vector', length(files))
snps_list = list('vector', length(files))
tissue_list = list('vector', length(files))

i = 1
for (tissue in files) {
	tissue_name = strsplit(tissue, '\\.')[[1]][1]
	ind = which(rownames(tissue_table) == tissue_name)
	count_list[[i]] = c(tissue_table$num_samples_with_geno[ind], length(unique(intrachrom_results[intrachrom_results$tissue == tissue_name,]$gene)), length(unique(intrachrom_results[intrachrom_results$tissue == tissue_name,]$snps)))
	gene_list[[i]] = as.character(intrachrom_results[intrachrom_results$tissue == tissue_name,]$gene)
	snps_list[[i]] = as.character(intrachrom_results[intrachrom_results$tissue == tissue_name,]$snps)
	tissue_list[[i]] = as.character(tissue_table$SMTSD[ind])

	i = i + 1
}
count_table = data.frame(do.call('rbind', count_list))
colnames(count_table) = c('num_samples', 'eGenes_5MB', 'eVariants_5MB')

tissue_list = unlist(tissue_list)
tissue_list = tissue_list[order(-count_table$num_samples)]
count_table = count_table[order(-count_table$num_samples),]

total = c(sum(count_table$num_samples), length(unique(unlist(gene_list))), length(unique(unlist(snps_list))))
count_table = rbind(count_table, total)
rownames(count_table) = c(tissue_list, 'Total')
write.table(count_table, file = '/tigress/BEE/RNAseq/Analysis/Figures/manuscript_figures/supplement_bundle/intrachrom_table.csv', sep = ',', quote=F)




