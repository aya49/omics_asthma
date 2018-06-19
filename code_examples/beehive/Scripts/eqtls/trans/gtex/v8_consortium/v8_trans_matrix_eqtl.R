##################################################
#  v8_trans_matrix_eqtl.R
#
#  $proj/Scripts/eqtls/trans/gtex/v8_consortium/v8_trans_matrix_eqtl.R
# 
#  This version performs the all-by-all trans-mapping analysis for v8 consortium.
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
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Cells_EBV-transformed_lymphocytes.v8.normalized_expression.bed.gz'
# args[2] = '5'
# args[3] = '1'
# args[4] = '25000'
# args[5] = '2.5e8'
# args[6] = 'Cells_EBV-transformed_lymphocytes'
# args[7] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
# args[8] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/v8/all-by-all/Cells_EBV-transformed_lymphocytes/'

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

# Filter by mappability of 0.8
mappability_list = read.table('/tigress/BEE/RNAseq/Output/processing/mappability/annotation/hg38_gene_mappability.txt', col.names = c('gene', 'mappability'), stringsAsFactors = FALSE)
rownames(mappability_list) = mappability_list$gene
expression_matrix = expression_matrix[sapply(rownames(expression_matrix), function(x) {(x %in% mappability_list$gene) && (mappability_list[x,2] >= 0.8)}),]

# Get the gene positions
gene_positions = expression_matrix[,c(4,1,2,3)]
colnames(gene_positions) = c('gene_id', 'chr', 'start', 'end')
expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]
# Save the list of genes tested for this tissue
if (!file.exists(paste0(out_dir, 'gene_positions.txt'))) {write.table(gene_positions, file = paste0(out_dir, 'gene_positions.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)}

# Read in the genotype positions - also generate permutations here with the same seed for consistency across runs
load(paste0(proj_dir, '/Data/Genotype/gtex/v8/allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_dosage_MAF_05_not_in_repeat.RData'))
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

# Get the SNP positions
snp_positions = data.frame(ID = rownames(genotype_matrix), chr = sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][1]}), pos = as.numeric(sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][2]})))

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

# Import MatrixEQTL
library(MatrixEQTL)
# source(paste0(proj_dir, '/Scripts/eqtls/MatrixEQTL_wrapper.R'))

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1e-5
cis_gene_threshold = 100000
useModel = modelLINEAR;

me_trans_list = list()
# me_trans_list_perm = list()
MAF_list = list()
total_ntests = 0
hist_counts = rep(0,100)
# hist_counts_perm = rep(0,100)
# i = 21905
for (i in c(1:nrow(genotype_matrix))) {
  print(i)
  inds = !is.na(genotype_matrix[i,])
  # inds_perm = !is.na(genotype_matrix_perm[i,])
  MAF_list[[i]] = sum(genotype_matrix[i,inds]) / (length(inds)*2)
  MAF_list[[i]] = min(MAF_list[[i]], 1 - MAF_list[[i]])
  # Skip if tissue-specific MAF is < 0.05
  if (MAF_list[[i]] < 0.05) {next}

  # Find genes that are cis to variant, and remove genes that are potentially cross-mapping
  expression_filtered = expression_matrix
  gene_pos_filtered = gene_positions
  cis_chr = gene_positions[gene_positions$chr == paste0('chr', chr_number),]
  cis_genes = cis_chr$gene[intersect(which(snp_positions$pos[i] >= (cis_chr$start - cis_gene_threshold)), which(snp_positions$pos[i] <= (cis_chr$end + cis_gene_threshold)))]
  if (length(cis_genes) > 0) {
    cross_mapping_genes = unique(unlist(sapply(cis_genes, function(x) {strsplit(cross_mapping_df[x,]$dest, ',')[[1]]})))
    if (sum(!is.na(cross_mapping_genes)) > 0) {
      # Filter out cross-mapping genes
      filtered_gene_set = setdiff(rownames(expression_matrix), cross_mapping_genes)
      print(length(filtered_gene_set))
      expression_filtered = expression_matrix[filtered_gene_set,]
      gene_pos_filtered = gene_positions[filtered_gene_set,]
    }
  }

  snps = SlicedData$new()
  snps$CreateFromMatrix(as.matrix(genotype_matrix[i,inds]));
  # snps_perm = SlicedData$new()
  # snps_perm$CreateFromMatrix(as.matrix(genotype_matrix_perm[i,inds_perm]));
  cvrt = SlicedData$new()
  cvrt$CreateFromMatrix(as.matrix(covars[,inds]));
  # cvrt_perm = SlicedData$new()
  # cvrt_perm$CreateFromMatrix(as.matrix(covars[,inds_perm]));
  gene = SlicedData$new()
  gene$CreateFromMatrix(as.matrix(expression_filtered[,inds]));
  # expression value not permuted, just to match the permuted SNP matrix
  # gene_perm = SlicedData$new()
  # gene_perm$CreateFromMatrix(as.matrix(expression_filtered[,inds_perm]));

  # This code necessary only for v8 - since in some cases the available genotypes will non-unique cov
  me = tryCatch({
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
      snpspos = snp_positions[i,],
      genepos = gene_pos_filtered,
      cisDist = cis_dist,
      pvalue.hist = 100,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)

    return(me)
  }, error = function(e) {
    temp_cvrt = covars[,inds]
    temp_cvrt = temp_cvrt[sapply(rownames(temp_cvrt), function(x) {length(unique(as.numeric(temp_cvrt[x,]))) > 1}),]
    cvrt$CreateFromMatrix(as.matrix(temp_cvrt));

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
      snpspos = snp_positions[i,],
      genepos = gene_pos_filtered,
      cisDist = cis_dist,
      pvalue.hist = 100,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)

    return(me)
  })

  # me_perm = Matrix_eQTL_main(
  #   snps = snps_perm,
  #   gene = gene_perm,
  #   cvrt = cvrt_perm,
  #   output_file_name = NULL,
  #   pvOutputThreshold = pvOutputThreshold,
  #   useModel = useModel,
  #   verbose = FALSE,
  #   output_file_name.cis = NULL,
  #   pvOutputThreshold.cis = pvOutputThreshold,
  #   snpspos = snp_positions[i,],
  #   genepos = gene_pos_filtered,
  #   cisDist = cis_dist,
  #   pvalue.hist = 100,
  #   min.pv.by.genesnp = FALSE,
  #   noFDRsaveMemory = FALSE)

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
MAF_df = data.frame(snp = row.names(genotype_matrix), MAF = unlist(MAF_list))
rownames(MAF_df) = MAF_df$snp
snps_tested = sum(MAF_df$MAF > 0.05)
gene_tested = nrow(gene_positions)

if (nrow(me_trans) > 0) {
  me_trans[,c('gene_chr', 'gene_start', 'gene_end')] = gene_positions[as.character(me_trans$gene),c('chr', 'start', 'end')]
  me_trans$MAF = MAF_df[as.character(me_trans$snps),'MAF']
}

# if (nrow(me_trans_perm) > 0) {
#   me_trans_perm[,c('gene_chr', 'gene_start', 'gene_end')] = gene_positions[as.character(me_trans_perm$gene),c('chr', 'start', 'end')]
#   me_trans_perm$MAF = MAF_df[as.character(me_trans_perm$snps),'MAF']
# }

# Check for existence of summary file to see if job completed
# save(me_trans, me_trans_perm, total_ntests, hist_counts, hist_counts_perm, snps_tested, gene_tested, file = paste0(out_dir, '/summary/chr', chr_number, '_part', part_number, '.RData'))
save(me_trans, total_ntests, hist_counts, snps_tested, gene_tested, file = paste0(out_dir, '/summary/chr', chr_number, '_part', part_number, '.RData'))
