##################################################
#  v9_trans_matrix_eqtl.R
#
#  $proj/Scripts/eqtls/trans/gtex/v8_consortium/v9_trans_matrix_eqtl.R
# 
#  This version performs the all-by-all trans-mapping analysis for v9 consortium.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# Example
# args = c(1:8)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Adipose_Subcutaneous.v8.normalized_expression.bed.gz'
# args[2] = '5'
# args[3] = '1'
# args[4] = '25000'
# args[5] = '2.5e8'
# args[6] = 'Adipose_Subcutaneous'
# args[7] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
# args[8] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/v9/all-by-all/Adipose_Subcutaneous/'

expression_file_location = args[1]
part_number = as.numeric(args[2])
chr_number = args[3]
partition_size = as.numeric(args[4])
cis_dist = as.numeric(args[5])
tissue_name = args[6]
cov_dir = args[7]
out_dir = args[8]

# Read in the expression files and the gene positions
header = readLines(gzfile(expression_file_location), n = 1)
header = strsplit(header, '\t')[[1]]
expression_matrix = read.csv(gzfile(expression_file_location, 'r'), sep = '\t', stringsAsFactors = FALSE)

colnames(expression_matrix) = header
rownames(expression_matrix) = expression_matrix$gene_id

# Gene filter #1 - Filter by mappability of 0.8
mappability_list = read.table('/tigress/BEE/RNAseq/Output/processing/mappability/annotation/hg38_gene_mappability.txt', col.names = c('gene', 'mappability'), stringsAsFactors = FALSE)
rownames(mappability_list) = mappability_list$gene
expression_matrix = expression_matrix[sapply(rownames(expression_matrix), function(x) {(x %in% mappability_list$gene) && (mappability_list[x,2] >= 0.8)}),]

# Get the gene positions
gene_positions = expression_matrix[,c(4,1,2,3)]
colnames(gene_positions) = c('gene_id', 'chr', 'start', 'end')
expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]

# Gene positions in the gtex expression file are wrong - let's fix that here
gencode = read.table('/tigress/BEE/RNAseq/Data/Resources/annotations/silver/gencode.v26.annotation.gtf', stringsAsFactors = F, sep = '\t', header = F)
gencode = gencode[gencode$V3 == 'gene',]
gene_ids = sapply(c(1:nrow(gencode)), function(x) {strsplit(strsplit(gencode[x,]$V9, ';')[[1]][1], ' ')[[1]][2]})
rownames(gencode) = gene_ids
gene_positions$start = gencode[rownames(gene_positions), 'V4']
gene_positions$end = gencode[rownames(gene_positions), 'V5']

# Gene filter #1.1 - Filter out genes with no annotated position
rows_to_keep = which(!is.na(gene_positions$gene_id))
expression_matrix = expression_matrix[rows_to_keep,]
gene_positions = gene_positions[rows_to_keep,]

# Save the list of genes tested for this tissue
if (!file.exists(paste0(out_dir, 'gene_positions.txt'))) {write.table(gene_positions, file = paste0(out_dir, 'gene_positions.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)}

# Read in the genotype positions - also generate permutations here with the same seed for consistency across runs
load(paste0(proj_dir, '/Data/Genotype/gtex/v9/allelic_dosage/GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_chr', chr_number, '_dosage_MAF_05_not_in_repeat.RData'))
print('Done loading expression and genotypes')
# This loads in the data frame genotype_mat_not_in_repeat
# set.seed(42)
# order = sample(ncol(genotype_mat_not_in_repeat))
# genotype_mat_not_in_repeat_perm = genotype_mat_not_in_repeat[,order]
# colnames(genotype_mat_not_in_repeat_perm) = colnames(genotype_mat_not_in_repeat)

# Get the appropriate partition
num_parts = ceiling(nrow(genotype_mat_not_in_repeat) / partition_size)
num_inds = ceiling(nrow(genotype_mat_not_in_repeat) / num_parts)
if (part_number == num_parts) {
  genotype_matrix = genotype_mat_not_in_repeat[c((((part_number-1)*partition_size)+1):nrow(genotype_mat_not_in_repeat)),]
  # genotype_matrix_perm = genotype_mat_not_in_repeat_perm[c((((part_number-1)*partition_size)+1):nrow(genotype_mat_not_in_repeat_perm)),]
} else {
  genotype_matrix = genotype_mat_not_in_repeat[c((((part_number-1)*partition_size)+1):(part_number*partition_size)),]
  # genotype_matrix_perm = genotype_mat_not_in_repeat_perm[c((((part_number-1)*partition_size)+1):(part_number*partition_size)),]
}

# Fix column name
colnames(genotype_matrix) = as.character(sapply(colnames(genotype_matrix), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))
# colnames(genotype_matrix_perm) = as.character(sapply(colnames(genotype_matrix_perm), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))

# Make sure the columns are the same
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]
# genotype_matrix_perm = genotype_matrix_perm[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix_perm))})]

# Fix the data type
genotype_matrix_temp = data.frame(lapply(genotype_matrix,as.numeric))
rownames(genotype_matrix_temp) = rownames(genotype_matrix)
genotype_matrix = genotype_matrix_temp

# genotype_matrix_temp = data.frame(lapply(genotype_matrix_perm,as.numeric))
# rownames(genotype_matrix_temp) = rownames(genotype_matrix_perm)
# genotype_matrix_perm = genotype_matrix_temp

# SNP filter #1 - For now, remove rows that have NA in them:
genotype_matrix = genotype_matrix[as.numeric(which(!is.na(rowSums(genotype_matrix)))),]
# SNP filter #2 - remove SNPs with tissue-specific MAF < 0.05
MAF_list = rowSums(genotype_matrix / (ncol(genotype_matrix)*2))
MAF_list = sapply(MAF_list, function(x) {min(x, 1-x)})
genotype_matrix = genotype_matrix[MAF_list >= 0.05,]

# Get the SNP positions
snp_positions = data.frame(ID = rownames(genotype_matrix), chr = sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][1]}), pos = as.numeric(sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][2]})))
print('Number of SNPs:')
print(nrow(snp_positions))

# Load in the covariates
suffix = '.v8.covariates.txt'
covars = read.csv(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
colnames(covars) = as.character(sapply(colnames(covars), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

# remove covariates that have only one value
covars = covars[sapply(rownames(covars), function(x) {length(unique(as.numeric(covars[x,]))) > 1}),]

# Read in cross-mapping conflict pairs
cross_mapping_df = read.table('/tigress/BEE/RNAseq/Output/processing/mappability/annotation/hg38_mappability_graph.txt', stringsAsFactors = FALSE, col.names = c('src', 'dest'))
rownames(cross_mapping_df) = cross_mapping_df$src

cis_gene_threshold = 100000
# Create intervals of cis-gene equivalence classes
chr_genes = gene_positions[gene_positions$chr == paste0('chr', chr_number),]
interval_df = data.frame(start = 1, end = 250000000, gene_list = '', stringsAsFactors = F)

for (gene in (rownames(chr_genes))) {
  # print(gene)
  # start creating intervals
  gene_int_start = gene_positions[gene, 'start'] - cis_gene_threshold
  # we're doing TSS distance - change this to end if it changes to gene distance
  gene_int_end = gene_positions[gene, 'start'] + cis_gene_threshold

  bool_1 = 1 * (interval_df$end > gene_int_start)
  bool_2 = 1 * (interval_df$start <= gene_int_end)
  ints_to_change = which(bool_1 + bool_2 == 2)

  if (length(ints_to_change) == 1) {
    # Add two intervals
    gene_list = as.character(interval_df[ints_to_change,]$gene_list)
    new_gene_list = paste0(gene_list, ',', gene)
    interval_df = rbind(interval_df, data.frame(start = gene_int_start, end = gene_int_end, gene_list = new_gene_list))
    interval_df = rbind(interval_df, data.frame(start = gene_int_end, end = interval_df[ints_to_change,]$end, gene_list = gene_list))
    interval_df[ints_to_change,]$end = gene_int_start
  } else {
    bool_3 = 1 * (interval_df$start <= gene_int_start)
    bool_4 = 1 * (interval_df$end > gene_int_end)
    left = which(bool_1 + bool_3 == 2)
    right = which(bool_2 + bool_4 == 2)
    interval_df = rbind(interval_df, data.frame(start = interval_df[right,]$start, end = gene_int_end, gene_list = paste0(interval_df[right,]$gene_list, ',', gene)))
    interval_df[right,]$start = gene_int_end
    interval_df = rbind(interval_df, data.frame(start = gene_int_start, end = interval_df[left,]$end, gene_list = paste0(interval_df[left,]$gene_list, ',', gene)))
    interval_df[left,]$end = gene_int_start
    rest = setdiff(ints_to_change, c(left, right))
    interval_df[rest,]$gene_list = sapply(interval_df[rest,]$gene_list, function(x) {paste0(x, ',', gene)})
  }
}
interval_df = interval_df[order(interval_df$start),]

# which intervals will we consider in this job?
bool_1 = 1 * (interval_df$end > min(snp_positions$pos))
bool_2 = 1 * (interval_df$start <= max(snp_positions$pos))
interval_in_range = interval_df[which(bool_1 + bool_2 == 2),]
rownames(interval_in_range) = c(1:nrow(interval_in_range))

# Import MatrixEQTL
library(MatrixEQTL)
# source(paste0(proj_dir, '/Scripts/eqtls/MatrixEQTL_wrapper.R'))

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1e-5
useModel = modelLINEAR;
me_trans_list = list()
# me_trans_list_perm = list()
MAF_list = list()
total_ntests = 0
hist_counts = rep(0,100)
# hist_counts_perm = rep(0,100)

# i = 21905
print("Number of intervals:")
print(nrow(interval_in_range))
for (i in c(1:nrow(interval_in_range))) {
  print(i)

  # Gene filter #2 - Find genes that are cis to variant, and remove genes that are potentially cross-mapping
  expression_filtered = expression_matrix
  gene_pos_filtered = gene_positions
  cis_genes = unlist(strsplit(as.character(interval_in_range[i,]$gene_list), ','))
  if (length(cis_genes) > 0) {
    cis_genes = cis_genes[2:length(cis_genes)]
    cross_mapping_genes = unique(unlist(sapply(cis_genes, function(x) {strsplit(cross_mapping_df[x,]$dest, ',')[[1]]})))
    if (sum(!is.na(cross_mapping_genes)) > 0) {
      # Filter out cross-mapping genes
      filtered_gene_set = setdiff(rownames(expression_matrix), cross_mapping_genes)
      print(length(filtered_gene_set))
      expression_filtered = expression_matrix[filtered_gene_set,]
      gene_pos_filtered = gene_positions[filtered_gene_set,]
    }
  }

  # Select the set of SNPs to test
  inds = which((1 * (snp_positions$pos >= interval_in_range[i,]$start)) + (1 * (snp_positions$pos < interval_in_range[i,]$end)) == 2)
  if (length(inds) == 0) {next}

  genotype_matrix_subset = genotype_matrix[inds,]
  snp_positions_subset = snp_positions[inds,]

  snps = SlicedData$new()
  snps$CreateFromMatrix(as.matrix(genotype_matrix_subset));
  # snps_perm = SlicedData$new()
  # snps_perm$CreateFromMatrix(as.matrix(genotype_matrix_perm[i,inds_perm]));
  cvrt = SlicedData$new()
  cvrt$CreateFromMatrix(as.matrix(covars));
  # cvrt_perm = SlicedData$new()
  # cvrt_perm$CreateFromMatrix(as.matrix(covars[,inds_perm]));
  gene = SlicedData$new()
  gene$CreateFromMatrix(as.matrix(expression_filtered));
  # expression value not permuted, just to match the permuted SNP matrix
  # gene_perm = SlicedData$new()
  # gene_perm$CreateFromMatrix(as.matrix(expression_filtered[,inds_perm]));

  me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = NULL,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      verbose = FALSE,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = pvOutputThreshold,
      snpspos = snp_positions_subset,
      genepos = gene_pos_filtered,
      cisDist = cis_dist,
      pvalue.hist = 100,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)

  total_ntests = total_ntests + me$trans$ntests
  hist_counts = hist_counts + me$trans$hist.counts
  # hist_counts_perm = hist_counts_perm + me_perm$trans$hist.counts
  if (me$trans$neqtls > 0) {
    me_trans_list[[i]] = me$trans$eqtls
  }
  # if (me_perm$trans$neqtls > 0) {
  #   me_trans_list_perm[[i]] = me_perm$trans$eqtls
  # }
}

me_trans = do.call('rbind', me_trans_list)
# me_trans_perm = do.call('rbind', me_trans_list_perm)

snps_tested = nrow(snp_positions)
gene_tested = nrow(gene_positions)

if (nrow(me_trans) > 0) {
  me_trans[,c('gene_chr', 'gene_start', 'gene_end')] = gene_positions[as.character(me_trans$gene),c('chr', 'start', 'end')]
  me_trans$MAF = MAF_list[as.character(me_trans$snps)]
}

# if (nrow(me_trans_perm) > 0) {
#   me_trans_perm[,c('gene_chr', 'gene_start', 'gene_end')] = gene_positions[as.character(me_trans_perm$gene),c('chr', 'start', 'end')]
#   me_trans_perm$MAF = MAF_df[as.character(me_trans_perm$snps),'MAF']
# }

# Check for existence of summary file to see if job completed
# save(me_trans, me_trans_perm, total_ntests, hist_counts, hist_counts_perm, snps_tested, gene_tested, file = paste0(out_dir, '/summary/chr', chr_number, '_part', part_number, '.RData'))
save(me_trans, total_ntests, hist_counts, snps_tested, gene_tested, file = paste0(out_dir, '/summary/chr', chr_number, '_part', part_number, '.RData'))
