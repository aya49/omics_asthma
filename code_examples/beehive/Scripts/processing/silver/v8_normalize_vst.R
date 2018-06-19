##################################################
#  v8_normalize_vst.R
#
#  $proj/Scripts/processing/silver/v8_normalize_vst.R
#
#  Matrix normalization script using vst in DESeq
#
#  Authors: Brian Jo
#
##################################################
library("DESeq")
library("vsn")
library("matrixStats")
library("ggplot2")

args <-commandArgs(TRUE)

# tissue = 'Whole_Blood'
tissue = args[1]

count_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/v8/raw/'
gtex_v8_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
out_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/v8/vst/'
fig_dir = '/tigress/BEE/RNAseq/Analysis/Figures/exploratory/v8_expression/'

gene_counts = read.table(paste0(count_dir, 'v8_RSEMv1.3.0_gene_count_', tissue, '.txt'), header=T, stringsAsFactors=F, sep='\t')
genes = gene_counts$gene
gene_counts = gene_counts[,c(2:ncol(gene_counts))]
rownames(gene_counts) = genes
colnames(gene_counts) = sapply(colnames(gene_counts), function(x) {paste0(strsplit(x, '\\.')[[1]][1], '-', strsplit(x, '\\.')[[1]][2])})

header = readLines(gzfile(paste0(gtex_v8_dir, tissue, '.v8.normalized_expression.bed.gz')), n = 1)
header = strsplit(header, '\t')[[1]]
normalized_expression = read.csv(gzfile(paste0(gtex_v8_dir, tissue, '.v8.normalized_expression.bed.gz'), 'r'), sep = '\t', stringsAsFactors = FALSE)
colnames(normalized_expression) = header
rownames(normalized_expression) = normalized_expression$gene_id
normalized_expression = normalized_expression[,c(5:ncol(normalized_expression))]

gene_count_comparison = gene_counts[rownames(normalized_expression),]

cds = newCountDataSet(round(gene_count_comparison), data.frame(row.names = colnames(gene_count_comparison), condition = rep('treated', ncol(gene_count_comparison))))
cds = estimateSizeFactors(cds)
cdsBlind = estimateDispersions(cds, method="blind")
vsd = varianceStabilizingTransformation(cdsBlind)

png(filename = paste0(fig_dir, tissue, '_gene_dispersion.png'), width = 500, height = 500)
plotDispEsts(cdsBlind)
dev.off()
png(filename = paste0(fig_dir, tissue, '_gene_logcount_meanSd.png'), width = 800, height = 500)
meanSdPlot(log2(counts(cds) + 1))
dev.off()
png(filename = paste0(fig_dir, tissue, '_gene_vst_meanSd.png'), width = 800, height = 500)
meanSdPlot(vsd)
dev.off()
png(filename = paste0(fig_dir, tissue, '_gene_Sds_histogram.png'), width = 800, height = 500)
hist(rowSds(exprs(vsd)), 50)
dev.off()

inds = sample(c(1:nrow(gene_count_comparison)), 200)
temp_mat = exprs(vsd[inds,colnames(normalized_expression)])
png(filename = paste0(fig_dir, tissue, '_vst_count_qn_comparison.png'), width = 600, height = 600)
qplot(as.numeric(unlist(temp_mat - rowMeans(temp_mat))), as.numeric(unlist(normalized_expression[inds,]))) + xlab('Mean-centered, variance-stabilized counts') + ylab('GTEx v8 normalized')
dev.off()

# Save the normalized count matrices
write.table(round(exprs(vsd), digits=4), quote=F, sep='\t', file=paste0(out_dir, 'v8_RSEMv1.3.0_gene_count_vst_', tissue, '.txt'))

# Read in the transcript counts
transcript_counts = read.table(paste0(count_dir, 'v8_RSEMv1.3.0_transcript_count_', tissue, '.txt'), header=T, stringsAsFactors=F, sep='\t')
transcripts = transcript_counts$transcript
genes = transcript_counts$gene
transcript_counts = transcript_counts[,c(3:ncol(transcript_counts))]
rownames(transcript_counts) = transcripts
colnames(transcript_counts) = sapply(colnames(transcript_counts), function(x) {paste0(strsplit(x, '\\.')[[1]][1], '-', strsplit(x, '\\.')[[1]][2])})

inds = as.numeric(unlist(sapply(rownames(normalized_expression), function(x) {which(genes == x)})))
transcript_counts_filtered = transcript_counts[inds,]
genes_filtered = genes[inds]
transcripts_filtered = transcripts[inds]
print(dim(transcript_counts_filtered))

# thresh = ncol(transcript_counts_filtered)/10
# inds2 = sapply(c(1:nrow(transcript_counts_filtered)), function(x) {sum(transcript_counts_filtered[x,] > 0) >= thresh})
# transcript_counts_filtered = transcript_counts_filtered[inds2,]
# genes_filtered = genes_filtered[inds2]
# transcripts_filtered = transcripts_filtered[inds2]
# print(dim(transcript_counts_filtered))

cds = newCountDataSet(round(transcript_counts_filtered), data.frame(row.names = colnames(transcript_counts_filtered), condition = rep('treated', ncol(transcript_counts_filtered))))
cds = estimateSizeFactors(cds)
cdsBlind = estimateDispersions(cds, method="blind")
vsd = varianceStabilizingTransformation(cdsBlind)

png(filename = paste0(fig_dir, tissue, '_transcript_dispersion.png'), width = 500, height = 500)
plotDispEsts(cdsBlind)
dev.off()
png(filename = paste0(fig_dir, tissue, '_transcript_logcount_meanSd.png'), width = 800, height = 500)
meanSdPlot(log2(counts(cds) + 1))
dev.off()
png(filename = paste0(fig_dir, tissue, '_transcript_vst_meanSd.png'), width = 800, height = 500)
meanSdPlot(vsd)
dev.off()
png(filename = paste0(fig_dir, tissue, '_transcript_Sds_histogram.png'), width = 800, height = 500)
hist(rowSds(exprs(vsd)), 50)
dev.off()

# Save the normalized count matrices
out_df = cbind(genes_filtered, transcripts_filtered)
out_df = cbind(out_df, round(exprs(vsd), digits=4))
write.table(out_df, row.names=F, quote=F, sep='\t', file=paste0(out_dir, 'v8_RSEMv1.3.0_transcript_count_vst_', tissue, '.txt'))
