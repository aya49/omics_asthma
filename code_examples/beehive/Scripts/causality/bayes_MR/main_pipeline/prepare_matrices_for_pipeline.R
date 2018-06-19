##################################################
#  prepare_matrices_for_pipeline.R
#
#  $proj/Scripts/causality/bayes_MR/prepare_matrices_for_pipeline.R
# 
#  This script prepares the expression and genotype matrices for the pipeline.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

# List of tissues to prepare
# tissue_list = c('Muscle_Skeletal', 'Adipose_Subcutaneous', 'Whole_Blood', 'Lung', 'Skin_Sun_Exposed_Lower_leg', 'Thyroid')
tissue_list = list.files('/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/', pattern = '.v8.normalized_expression.bed.gz.tbi')
tissue_list = as.character(sapply(tissue_list, function(x) {strsplit(x, '.v8')[[1]][1]}))
expression_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
cov_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
save_dir = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/prep_files/'

# Read in Gencode annotations
gencode = read.table('/tigress/BEE/RNAseq/Data/Resources/annotations/silver/gencode.v26.annotation.gtf', stringsAsFactors = F, sep = '\t', header = F)
gencode = gencode[gencode$V3 == 'gene',]
gene_ids = sapply(c(1:nrow(gencode)), function(x) {strsplit(strsplit(gencode[x,]$V9, ';')[[1]][1], ' ')[[1]][2]})
rownames(gencode) = gene_ids

# Read in the list of AA and EA individuals 
pop_clust_dir = '/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/support_files/pop_clust/'
AA_list_1 = read.table(paste0(pop_clust_dir, 'AA_1.txt'), stringsAsFactors = FALSE)
AA_list_2 = read.table(paste0(pop_clust_dir, 'AA_2.txt'), stringsAsFactors = FALSE)
EA_list_1 = read.table(paste0(pop_clust_dir, 'EA_1.txt'), stringsAsFactors = FALSE)
EA_list_2 = read.table(paste0(pop_clust_dir, 'EA_2.txt'), stringsAsFactors = FALSE)
EA_list_3 = read.table(paste0(pop_clust_dir, 'EA_3.txt'), stringsAsFactors = FALSE)
AA_list = c(AA_list_1$V1, AA_list_2$V1)
EA_list = c(EA_list_1$V1, EA_list_2$V1, EA_list_3$V1)

for (tissue in tissue_list) {
	expression_file_location = paste0(expression_dir, tissue, '.v8.normalized_expression.bed.gz')
	header = readLines(gzfile(expression_file_location), n = 1)
	header = strsplit(header, '\t')[[1]]
	expression_matrix = read.csv(gzfile(expression_file_location, 'r'), sep = '\t', stringsAsFactors = F)

	colnames(expression_matrix) = header
	rownames(expression_matrix) = expression_matrix$gene_id

	# Filter out genes with low mappability
	mappability_list = read.table('/tigress/BEE/RNAseq/Output/processing/mappability/annotation/hg38_gene_mappability.txt', col.names = c('gene', 'mappability'), stringsAsFactors = FALSE)
	rownames(mappability_list) = mappability_list$gene
	expression_matrix = expression_matrix[sapply(rownames(expression_matrix), function(x) {(x %in% mappability_list$gene) && (mappability_list[x,2] >= 0.8)}),]

	# Only take the genes that are in the filtered list:
	filtered_gene_list = read.table('/tigress/BEE/RNAseq/Data/Resources/annotations/silver/gene_filter/Bayes_MR_list.txt', stringsAsFactors = FALSE)
	expression_matrix = expression_matrix[sapply(rownames(expression_matrix), function(x) {x %in% filtered_gene_list$V1}),]

	gene_positions = expression_matrix[,c(4,1,2,3)]
	colnames(gene_positions) = c('gene_id', 'chr', 'start', 'end')
	expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]

	# Gene positions in the v8 GTEx gct file are incorrect - need to put in correct gene locations
	gene_positions$start = gencode[rownames(gene_positions), 'V4']
	gene_positions$end = gencode[rownames(gene_positions), 'V5']

	# Read in the covariates
	suffix = '.v8.covariates.txt'
	covars = read.csv(paste0(cov_dir, tissue, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	colnames(covars) = as.character(sapply(colnames(covars), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))
	rownames(covars) = covars$ID
	covars = covars[,colnames(expression_matrix)]

	# invert and mean-center covars
	inv_cov = t(covars)
	cov = data.frame(inv_cov)

	# Read in the list of AA and EA individuals 
	AA_list_in_tissue = AA_list[sapply(AA_list, function(x) {x %in% colnames(expression_matrix)})]
	EA_list_in_tissue = EA_list[sapply(EA_list, function(x) {x %in% colnames(expression_matrix)})]

	# Prepare the expression, gene position, and cov matrices for pipeline
	expression_matrix_AA = expression_matrix[,AA_list_in_tissue]
	expression_matrix_EA = expression_matrix[,EA_list_in_tissue]
	cov_AA = cov[AA_list_in_tissue,]
	cov_EA = cov[EA_list_in_tissue,]

	# remove covs columns that are not unique
	uniq = which(sapply(c(1:ncol(cov_AA)), function(x) {length(unique(cov_AA[,x]))}) == 1)
	if (length(uniq > 0)) {cov_AA = subset(cov_AA, select = -uniq)}
	uniq = which(sapply(c(1:ncol(cov_EA)), function(x) {length(unique(cov_EA[,x]))}) == 1)
	if (length(uniq > 0)) {cov_EA = subset(cov_EA, select = -uniq)}

	# save the matrices for each tissue
	save(expression_matrix, expression_matrix_EA, expression_matrix_AA, gene_positions, cov, cov_EA, cov_AA, file = paste0(save_dir, tissue, '.RData'))
}
