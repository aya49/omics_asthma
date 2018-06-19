##################################################
#  MatrixEQTL_FDR_control_PEER_increments_final.R
#
#  $proj/Scripts/eqtls/trans/gtex/MatrixEQTL_FDR_control_PEER_increments_final.R
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
# original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/'
# tissue_name = 'thyroid'
# write_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/'
# nPEER = '10'

original_dir = args[1]
tissue_name = args[2]
write_dir = args[3]
nPEER = args[4]

library(ggplot2)
library(dplyr)

list_files = list.files(path=paste(original_dir, tissue_name, '/', sep=''),pattern=paste0('PEER', nPEER))
print(length(list_files))

# Required data to aggregate
cumulative_me_trans_total = vector('list', length(list_files))
cumulative_me_trans_prot_only = vector('list', length(list_files))
cumulative_me_trans_linc_only = vector('list', length(list_files))

# cumulative_me_trans_permute_total = vector('list', length(list_files))
all_hist_vals_trans_total = rep(0,100)
# all_hist_vals_trans_permute_total = rep(0,100)
all_n_tests_trans_total = 0
all_n_tests_prot_total = 0
all_n_tests_linc_total = 0
# all_n_tests_trans_permute_total = 0

# Load in gene annotation and the tissue expression file
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

gene_positions = read.table(file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gene_positions = gene_positions[,c('gene_id', 'chr', 'start', 'end')]
gene_positions = mutate(gene_positions, chr = paste("chr", chr, sep = ""))
gene_positions = select(gene_positions, c(gene_id, chr, start, end))
expression_matrix = expression_matrix[rownames(expression_matrix) %in% gene_positions$gene_id,]
gene_positions = gene_positions[gene_positions$gene_id %in% rownames(expression_matrix),]
rownames(gene_positions) = gene_positions$gene_id
gene_positions = gene_positions[rownames(expression_matrix),]

v6p_consortium_gene_table = read.table(paste0(proj_dir, '/Data/Resources/annotations/gtex.v6p.genes.with.without.cis.eqtls.tested.txt'), header=T, sep='\t')
rownames(v6p_consortium_gene_table) = v6p_consortium_gene_table$GENE
gene_positions$TYPE = v6p_consortium_gene_table[rownames(gene_positions),'TYPE']

i = 1
for (item in list_files) {
    load(paste(original_dir, tissue_name, '/', item, sep=""))
    
    cumulative_me_trans_total[[i]] = me$trans$eqtls

    chr = me$trans$eqtls$snp_chr[1]
    trans_gene_positions = gene_positions[gene_positions$chr != chr,]
    sum_prot = sum(trans_gene_positions$TYPE == 'protein_coding')
    sum_linc = sum(trans_gene_positions$TYPE == 'lincRNA')

    cumulative_me_trans_prot_only[[i]] = me$trans$eqtls[(gene_positions[as.character(me$trans$eqtls$gene),'TYPE'] == 'protein_coding'),]
    cumulative_me_trans_linc_only[[i]] = me$trans$eqtls[(gene_positions[as.character(me$trans$eqtls$gene),'TYPE'] == 'lincRNA'),]

    # cumulative_me_trans_permute_total[[i]] = me_permute$trans$eqtls
    all_hist_vals_trans_total = all_hist_vals_trans_total+me$trans$hist.counts
    # all_hist_vals_trans_permute_total = all_hist_vals_trans_permute_total+me_permute$trans$hist.counts
    all_n_tests_trans_total = all_n_tests_trans_total+me$trans$ntests
    all_n_tests_prot_total = all_n_tests_prot_total+(me$trans$ntests * (sum_prot / (sum_prot + sum_linc)))
    all_n_tests_linc_total = all_n_tests_linc_total+(me$trans$ntests * (sum_linc / (sum_prot + sum_linc)))
    # print(all_n_tests_prot_total)
    # print(all_n_tests_linc_total)
    # all_n_tests_trans_permute_total = all_n_tests_trans_permute_total+me_permute$trans$ntests
    i = i+1
}

cumulative_me_trans_total = do.call("rbind", cumulative_me_trans_total)
cumulative_me_trans_prot_only = do.call("rbind", cumulative_me_trans_prot_only)
cumulative_me_trans_linc_only = do.call("rbind", cumulative_me_trans_linc_only)

# cumulative_me_trans_permute_total = do.call("rbind", cumulative_me_trans_permute_total)
print(dim(cumulative_me_trans_total))
# print(dim(cumulative_me_trans_permute_total))

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

cumulative_me_trans_total = cumulative_me_trans_total[order(cumulative_me_trans_total$pvalue),]
cumulative_me_trans_total$FDR = BH_FDR(cumulative_me_trans_total$pvalue, all_n_tests_trans_total)

cumulative_me_trans_prot_only = cumulative_me_trans_prot_only[order(cumulative_me_trans_prot_only$pvalue),]
cumulative_me_trans_prot_only$FDR = BH_FDR(cumulative_me_trans_prot_only$pvalue, all_n_tests_prot_total)

cumulative_me_trans_linc_only = cumulative_me_trans_linc_only[order(cumulative_me_trans_linc_only$pvalue),]
cumulative_me_trans_linc_only$FDR = BH_FDR(cumulative_me_trans_linc_only$pvalue, all_n_tests_linc_total)

# Now we can finish the post-processing steps:
# Make necessary directories:
dir.create(file.path(write_dir, 'pvalue_histograms'), showWarnings = FALSE)
dir.create(file.path(write_dir, 'eqtl_list'), showWarnings = FALSE)

## Generate p-value histograms
# print('Generateing p-value histograms:')
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans.png', sep=''), plot=g)
# g = ggplot(data.frame(ind = (c(1:100)/100), counts = all_hist_vals_trans_permute_total), aes(ind, counts)) + geom_bar(stat="identity") + xlab('P-value distribution') + ggtitle(paste('Histogram of of P-values for', all_n_tests_trans_total, 'tests', sep=' '))
# ggsave(filename = paste(write_dir, 'pvalue_histograms/', tissue_name, '_', nPEER, '_hist_trans_permute.png', sep=''), plot=g)

# Save necessary files and write out summary
# As for actual eQTL list, write out only for pvalue<save_threshold, as requested by Yuan

save(cumulative_me_trans_total, cumulative_me_trans_prot_only, cumulative_me_trans_linc_only, file=paste(write_dir, 'eqtl_list/', tissue_name, '_', nPEER, "_trans_list_by_pvalue.RData", sep=""))

# save_threshold = 0.5
# null_df = data.frame(snps = character(), gene = character(), statistic = numeric(), pvalue = numeric(), FDR = numeric(), beta = numeric(), snp_chr = character(), snp_pos = numeric(), gene_chr = character(), gene_start = numeric(), gene_end = numeric())

# if (sum(cumulative_me_trans_total$FDR<save_threshold) > 0) {
#     write.table(cumulative_me_trans_total[c(1:sum(cumulative_me_trans_total$FDR<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, '_', nPEER, "_trans_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
# } else {
#     write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, '_', nPEER, "_trans_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
# }

# if (sum(cumulative_me_trans_prot_only$FDR<save_threshold) > 0) {
#     write.table(cumulative_me_trans_prot_only[c(1:sum(cumulative_me_trans_prot_only$FDR<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, '_', nPEER, "_trans_prot_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
# } else {
#     write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, '_', nPEER, "_trans_prot_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
# }

# if (sum(cumulative_me_trans_linc_only$FDR<save_threshold) > 0) {
#     write.table(cumulative_me_trans_linc_only[c(1:sum(cumulative_me_trans_linc_only$FDR<save_threshold)),], file=paste(write_dir, 'eqtl_list/', tissue_name, '_', nPEER, "_trans_linc_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
# } else {
#     write.table(null_df, file=paste(write_dir, 'eqtl_list/', tissue_name, '_', nPEER, "_trans_linc_list_by_pvalue.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)    
# }
