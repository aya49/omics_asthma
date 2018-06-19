##################################################
#  Final_trans_list_PEER.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/Final_trans_list_PEER.R
# 
#  Finalizes the list of trans-eQTLs and draws eQTL boxplots
#
#  Author: Brian Jo
#
##################################################

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# in_dir = paste0(proj_dir, '/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/')
# dist_thresh = 1e5
# tissue_name = 'smallintestineterminalileum'
in_dir = args[1]
dist_thresh = as.numeric(args[2])
tissue_name = args[3]

#library(pracma)
library(dplyr)

crossmap_data = paste0(proj_dir, '/Data/Resources/annotations/pairwise_conflict.RData')
if (file.exists(crossmap_data)) {
	load(crossmap_data)
} else {
	pairwise_conflict_file = paste0(proj_dir, '/Data/Resources/annotations/pairwise_conflict.txt')
	conflict_df = read.table(pairwise_conflict_file, stringsAsFactors=FALSE, header=FALSE, sep='\t')
	colnames(conflict_df) = c('src', 'dest')
	save(conflict_df, file = crossmap_data)
}

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

# Define a function that outputs a vector of booleans when given a dataframe of associations
crossmap_filter = function (total_trans_df, conflict_df, gene_position_df) {
	trans_bool = rep(TRUE, nrow(total_trans_df))
	conflict_list = list()
	for (i in c(1:nrow(total_trans_df))) {
		print(i)
		trans_gene = as.character(total_trans_df$gene[i])
		if (length(conflict_list[[trans_gene]]) == 0) {
			target_gene_list = unique(c(conflict_df$dest[conflict_df$src == trans_gene], conflict_df$src[conflict_df$dest == trans_gene]))
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
	return(trans_bool)
}

# Process each tissue separately
list_files = list.files(path=paste0(in_dir, 'eqtl_list/'), pattern='_trans_list_by_pvalue.RData')
tissue_specific = list_files[sapply(list_files, function(x) {grepl(tissue_name, x)})]

dir.create(file.path(in_dir, 'eqtl_list_post_crossmapping'), showWarnings = FALSE)

conflict_list = list()

PEER_list = as.character(sapply(tissue_specific, function(x) {strsplit(x, '_')[[1]][2]}))
null_df = data.frame(snps = character(), gene = character(), statistic = numeric(), pvalue = numeric(), FDR = numeric(), beta = numeric(), snp_chr = character(), snp_pos = numeric(), gene_chr = character(), gene_start = numeric(), gene_end = numeric())

for (peer in PEER_list) {
	print(peer)
	file = paste0(in_dir, 'eqtl_list/', paste(tissue_name, peer, 'trans_list_by_pvalue.RData', sep='_'))
	# total_trans_df = read.table(file, stringsAsFactors=FALSE, header=TRUE, sep='\t')
	load(file)

	# Apply cross-mapping filter - protein_coding
	if (sum(cumulative_me_trans_prot_only$FDR <= 0.5) > 0) {
		total_trans_df = cumulative_me_trans_prot_only[cumulative_me_trans_prot_only$FDR <= 0.5,]
		total_trans_df = total_trans_df[crossmap_filter(total_trans_df, conflict_df, gene_position_df),]
		write.table(total_trans_df, file=paste0(in_dir, 'eqtl_list_post_crossmapping/', paste(tissue_name, peer, 'trans_prot_list_FDR_0.5.txt', sep='_')), quote = FALSE, sep = "\t", row.names = FALSE)
	} else {
		write.table(null_df, file=paste0(in_dir, 'eqtl_list_post_crossmapping/', paste(tissue_name, peer, 'trans_prot_list_FDR_0.5.txt', sep='_')), quote = FALSE, sep = "\t", row.names = FALSE)
	}
	# Apply cross-mapping filter - lincRNA
	if (sum(cumulative_me_trans_linc_only$FDR <= 0.5) > 0) {
		total_trans_df = cumulative_me_trans_linc_only[cumulative_me_trans_linc_only$FDR <= 0.5,]
		total_trans_df = total_trans_df[crossmap_filter(total_trans_df, conflict_df, gene_position_df),]
		write.table(total_trans_df, file=paste0(in_dir, 'eqtl_list_post_crossmapping/', paste(tissue_name, peer, 'trans_linc_list_FDR_0.5.txt', sep='_')), quote = FALSE, sep = "\t", row.names = FALSE)
	} else {
		write.table(null_df, file=paste0(in_dir, 'eqtl_list_post_crossmapping/', paste(tissue_name, peer, 'trans_linc_list_FDR_0.5.txt', sep='_')), quote = FALSE, sep = "\t", row.names = FALSE)
	}
	# Apply cross-mapping filter - joint
	if (sum(cumulative_me_trans_total$FDR <= 0.5) > 0) {
		total_trans_df = cumulative_me_trans_total[cumulative_me_trans_total$FDR <= 0.5,]
		total_trans_df = total_trans_df[crossmap_filter(total_trans_df, conflict_df, gene_position_df),]
		write.table(total_trans_df, file=paste0(in_dir, 'eqtl_list_post_crossmapping/', paste(tissue_name, peer, 'trans_list_FDR_0.5.txt', sep='_')), quote = FALSE, sep = "\t", row.names = FALSE)
	} else {
		write.table(null_df, file=paste0(in_dir, 'eqtl_list_post_crossmapping/', paste(tissue_name, peer, 'trans_list_FDR_0.5.txt', sep='_')), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}

# Only for the final PEER factor, and for the joint analysis, produce the cross-mapped list for p-value threshold of 1e-5 (for joint-tissue analysis)
# dir.create(file.path(in_dir, 'files_for_ben'), showWarnings = FALSE)
# total_trans_df = cumulative_me_trans_total
# total_trans_df = total_trans_df[crossmap_filter(total_trans_df, conflict_df, gene_position_df),]
# save(total_trans_df, file=paste0(in_dir, 'files_for_ben/', paste(tissue_name, peer, 'trans_list_pval_1e-5.RData', sep='_')))

# gene set size - not necessary if the FDR step passes in gene set size
expression_matrix = read.csv(file = paste0(proj_dir, '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/', tissue_name, '_v6p_consortium_autosomes_normalized.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
mappability_cutoff = 0.8
mappability_list = read.table(paste0(proj_dir, '/Data/Resources/annotations/avg_mappability_Exon_UTR.txt'), stringsAsFactors=F)
rownames(mappability_list) = mappability_list$V1
# Arbitrary 0.8 cutoff - can be modified for a different threshold
mappability_list = mappability_list[(mappability_list$V2>mappability_cutoff),]
# Filter out genes with low mappability
expression_matrix = expression_matrix[rownames(expression_matrix) %in% rownames(mappability_list),]

dir.create(file.path(in_dir, 'run_summary'), showWarnings = FALSE)
summary_file = paste(in_dir, 'run_summary/', tissue_name, '_run_summary.txt', sep='')
fileConn<-file(summary_file)
line1 = paste('MatrixEQTL run summary for:', tissue_name, sep='\t')
line2 = paste('Number of genes tested:', nrow(expression_matrix), sep='\t')
line3 = paste('Number of SNPs tested:', snp_set_size_final_total, sep='\t')
line4 = paste('Number of trans (off-chromosome) tests:', all_n_tests_trans_total, sep='\t')
line5 = paste('Total number of trans-eQTLs before cross-mapping filter (FDR 0.1):', sum(cumulative_me_trans_total$FDR <= 0.1), sep='\t')
# total_trans_df only exists if sum(cumulative_me_trans_total$FDR <= 0.5) > 0, check above.
if (sum(cumulative_me_trans_total$FDR <= 0.5) > 0) {
	line6 = paste('Total number of trans-eQTLs after cross-mapping filter (FDR 0.1):', sum(total_trans_df$FDR <= 0.1), sep='\t')
	line7 = paste('Total number of trans-eGenes after cross-mapping filter (FDR 0.1):', length(unique(total_trans_df[c(0:sum(total_trans_df$FDR <= 0.1)),]$gene)), sep='\t')
} else {
	line6 = paste('Total number of trans-eQTLs after cross-mapping filter (FDR 0.1):', 0, sep='\t')
	line7 = paste('Total number of trans-eGenes after cross-mapping filter (FDR 0.1):', 0, sep='\t')
}
writeLines(c(line1, line2, line3, line4, line5, line6, line7), fileConn)
close(fileConn)
