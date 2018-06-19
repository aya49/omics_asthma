##################################################
#  v8_ld_pruned_subset_allelic.R
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/v8_ld_pruned_subset_allelic.R
# 
#  This script saves a separate file of genotypes that are not LD-pruned
#
#  Author: Brian Jo
#
##################################################

proj_dir = Sys.getenv('proj')

for (chr_number in c(c(1:22), 'X')) {
	print(chr_number)
	genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_dosage_MAF_05.txt')
	genotype_matrix_master = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	# Add in the LD filter:
	pruned_list = as.character(read.table(paste0('/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/ld_prune/chr', chr_number, '_MAF_05.in'))$V1)
	rownames(genotype_matrix_master) = genotype_matrix_master$X
	genotype_matrix_master = genotype_matrix_master[pruned_list,]
	genotype_matrix_master = genotype_matrix_master[,c(2:ncol(genotype_matrix_master))]
	colnames(genotype_matrix_master) = as.character(sapply(colnames(genotype_matrix_master), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))
	save(genotype_matrix_master, file = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_05_ld_pruned.RData'))

	genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_dosage_MAF_01.txt')
	genotype_matrix_master = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	# Add in the LD filter:
	pruned_list = as.character(read.table(paste0('/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/ld_prune/chr', chr_number, '_MAF_01.in'))$V1)
	rownames(genotype_matrix_master) = genotype_matrix_master$X
	genotype_matrix_master = genotype_matrix_master[pruned_list,]
	genotype_matrix_master = genotype_matrix_master[,c(2:ncol(genotype_matrix_master))]
	colnames(genotype_matrix_master) = as.character(sapply(colnames(genotype_matrix_master), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))
	save(genotype_matrix_master, file = paste0(proj_dir, '/Data/Genotype/gtex/v8/ld_prune/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_MAF_01_ld_pruned.RData'))
}

