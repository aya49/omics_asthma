##################################################
#  Final_trans_list_boxplots.R
#
#  $proj/Scripts/eqtls/trans/gtex/Final_trans_list_boxplots.R
# 
#  Finalizes the list of trans-eQTLs and draws eQTL boxplots
#
#  Author: Brian Jo
#
##################################################

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# in_dir = paste0(proj_dir, '/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/')
# FDR_thresh = 0.10
# dist_thresh = 1e5
# method = 'all-by-all'
in_dir = paste0(proj_dir, args[1])
FDR_thresh = as.numeric(args[2])
dist_thresh = as.numeric(args[3])
method = args[4]

library(ggplot2)
library(pracma)
library(dplyr)

if (file.exists(paste0(in_dir, 'Final_trans_eQTL_list_', FDR_thresh, '.txt'))) {
	total_trans_df_final = read.table(paste0(in_dir, 'Final_trans_eQTL_list_', FDR_thresh, '.txt'), stringsAsFactors=FALSE, header=TRUE, sep='\t')
} else {

	pairwise_conflict_file = paste0(proj_dir, '/Data/Resources/annotations/pairwise_conflict.txt')
	conflict_df = read.table(pairwise_conflict_file, stringsAsFactors=FALSE, header=FALSE, sep='\t')

	list_files = list.files(path=paste0(in_dir, 'eqtl_list/'), pattern='_trans_list_by_pvalue.txt')
	total_trans_df = vector('list', length(list_files))

	i = 1
	for (file in list_files) {
		print(i)
		tissue_name = strsplit(file, '_trans_list_by_pvalue.txt')[[1]][1]
		temp_df = read.table(paste0(in_dir, 'eqtl_list/', file), stringsAsFactors=FALSE, header=TRUE, sep='\t')
		if (sum(temp_df$FDR<=FDR_thresh) > 0) {
			total_trans_df[[i]] = temp_df[(temp_df$FDR<=FDR_thresh),]
			total_trans_df[[i]]$tissue = tissue_name
		}
		i = i+1
	}

	total_trans_df = do.call('rbind', total_trans_df)

	# read in gene locations - processing on the gene annotations:
	gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt')
	gene_position_df = read.table(gene_position_file, stringsAsFactors=FALSE, header=TRUE, sep='\t')
	gene_position_df = gene_position_df[,c('gene_id', 'chr', 'start', 'end')]
	gene_position_df = mutate(gene_position_df, chr = paste0("chr", chr))
	gene_position_df = select(gene_position_df, c(gene_id, chr, start, end))
	gene_position_df = gene_position_df[!duplicated(gene_position_df$gene_id),]
	gene_position_df$start = as.numeric(gene_position_df$start)
	gene_position_df$end = as.numeric(gene_position_df$end)
	rownames(gene_position_df) = gene_position_df$gene_id

	# for each row, determine whether to include to the final list or not:
	trans_bool = rep(TRUE, dim(total_trans_df)[1])
	conflict_list = list()

	for (i in c(1:dim(total_trans_df)[1])) {
		print(i)
		trans_gene = total_trans_df$gene[i]
		if (length(conflict_list[[trans_gene]]) == 0) {
			target_gene_list = unique(c(conflict_df$V2[conflict_df$V1 == trans_gene], conflict_df$V1[conflict_df$V2 == trans_gene]))
			target_gene_list = target_gene_list[target_gene_list %in% rownames(gene_position_df)]
			conflict_list[[trans_gene]] = target_gene_list
		}
		test_gene_set = conflict_list[[trans_gene]][gene_position_df[conflict_list[[trans_gene]], 'chr'] == total_trans_df$snp_chr[i]]
		if (length(test_gene_set)>0) {
			bool1 = sapply(gene_position_df[test_gene_set,]$start, function(x) {(total_trans_df$snp_pos[i] - x + dist_thresh)>=0})
			bool2 = sapply(gene_position_df[test_gene_set,]$end, function(x) {(total_trans_df$snp_pos[i] - x - dist_thresh)<=0})
			if (sum(rowSums(cbind(bool1,bool2)) == 2) > 0) {
				trans_bool[i] = FALSE
			}
		}
	}

	# Write out the final table of trans-eQTLs and also the summary of cross-mapping filter
	tissue_list = as.character(sapply(list_files, function(x) {strsplit(x, '_trans_list_by_pvalue.txt')[[1]][1]}))
	num_trans = as.numeric(sapply(tissue_list, function(x) {sum(total_trans_df$tissue == x)}))
	num_trans = c(num_trans, dim(total_trans_df)[1])

	total_trans_df_final = total_trans_df[trans_bool,]
	num_trans_filter = as.numeric(sapply(tissue_list, function(x) {sum(total_trans_df_final$tissue == x)}))
	num_trans_filter = c(num_trans_filter, dim(total_trans_df_final)[1])

	cross_mapping_df = data.frame(num_trans = num_trans, num_trans_filter = num_trans_filter)
	rownames(cross_mapping_df) = c(tissue_list, 'total')
	# proportion - will naturally have a NaN if there were 0
	cross_mapping_df$proportion = cross_mapping_df$num_trans_filter / cross_mapping_df$num_trans

	tissue_table = read.table(paste0(proj_dir, '/Data/Resources/gtex/tables/tissue_table.txt'), sep = '\t', stringsAsFactors = FALSE)
	cross_mapping_df$num_samples = tissue_table[rownames(cross_mapping_df), 'num_samples_with_geno']
	cross_mapping_df$num_samples[dim(cross_mapping_df)[1]] = sum(cross_mapping_df$num_samples[c(1:(dim(cross_mapping_df)[1]-1))])

	write.table(total_trans_df_final, file=paste0(in_dir, 'Final_trans_eQTL_list_', FDR_thresh, '.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
	write.table(cross_mapping_df, file=paste0(in_dir, 'Cross_mapping_filter_summary_', FDR_thresh, '.txt'), quote = FALSE, sep = "\t", row.names = TRUE)
}
