##################################################
#  collect_genotypes_not_in_repeat.R
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/collect_genotypes_not_in_repeat.R
# 
#  This script saves a separate file of genotypes that are not in repeat elements
#
#  Author: Brian Jo
#
##################################################

# Record the genotypes for both original and permuted data
in_path = '/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous'
out_path = '/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous'
in_path_p = '/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous_permute'
out_path_p = '/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous_permute'
# Take only SNPs that are not in a repeated element, and also pass 
for (i in c(1:22)) {
  genotype_metadata = read.table(paste('/tigress/BEE/eQTLs/Data/Genotype/GTEx/SNP_metadata/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr', i, '.metadata.expanded.txt', sep=''), header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  genotype_metadata = genotype_metadata[genotype_metadata$IN_REPEAT == FALSE,]
  # Used to have a MAF filter - discontinued to apply MAF at tissue level for future tests
  #genotype_metadata = genotype_metadata[genotype_metadata$maf05_FILTER == '1',c('dbSNPID', 'IN_REPEAT')]
  genotype_matrix = read.table(paste(in_path, '/GTEx_genotypes_maf05_continuous_Chr', i, '_Final.txt', sep=''), header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  genotype_matrix = genotype_matrix[(rownames(genotype_matrix) %in% genotype_metadata$dbSNPID),]
  write.table(genotype_matrix, paste(out_path, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', i, '_Final.txt', sep=''), quote = FALSE, sep = "\t")
  genotype_matrix = read.table(paste(in_path_p, '/GTEx_genotypes_maf05_continuous_Chr', i, '_Final.txt', sep=''), header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  genotype_matrix = genotype_matrix[(rownames(genotype_matrix) %in% genotype_metadata$dbSNPID),]
  write.table(genotype_matrix, paste(out_path_p, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', i, '_Final.txt', sep=''), quote = FALSE, sep = "\t")
  print(i)
}