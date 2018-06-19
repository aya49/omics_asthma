##################################################
#  v9_collect_genotypes_not_in_repeat.R
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/v9_collect_genotypes_not_in_repeat.R
# 
#  This script saves a separate file of genotypes that are not in repeat elements
#
#  Author: Brian Jo
#
##################################################

args <-commandArgs(TRUE)
i = args[1]

# RepeatMasker has been unpacked to Della:
repeat_masker_f = '/tigress/BEE/RNAseq/Data/Resources/RepeatMasker/repeat_masker_v38.RData'
if (!file.exists(repeat_masker_f)) {
  # Prepare the R version of repeatmasker file
  repeat_masker = read.table('/tigress/BEE/RNAseq/Data/Resources/RepeatMasker/hg38.fa.out', header=F, stringsAsFactors=F, sep='', skip=2, fill=T)
  colnames(repeat_masker)[5] = 'chr'
  colnames(repeat_masker)[6] = 'start'
  colnames(repeat_masker)[7] = 'end'
  save(repeat_masker, file = repeat_masker_f)
} else {
  load(repeat_masker_f)
}

# Record the genotypes for both original and permuted data
file_path = '/tigress/BEE/RNAseq/Data/Genotype/gtex/v9/allelic_dosage/'

# Take only SNPs that are not in a repeated element, and also pass 
# for (i in c(c(1:22), 'X')) {
  # First, for MAF 01
  genotype_mat = read.table(paste0(file_path, 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr',i,'_dosage_MAF_01.txt'), header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  locs = as.numeric(sapply(genotype_mat$X, function(x) {strsplit(x, '_')[[1]][2]}))
  repeat_subset = repeat_masker[repeat_masker$chr == paste0('chr', i),]
  not_in_repeat = sapply(locs, function(x) {sum(x >= repeat_subset$start) == sum(x > repeat_subset$end)})
  genotype_mat_not_in_repeat = genotype_mat[not_in_repeat,]
  rownames(genotype_mat_not_in_repeat) = genotype_mat_not_in_repeat$X
  genotype_mat_not_in_repeat = genotype_mat_not_in_repeat[,c(2:ncol(genotype_mat_not_in_repeat))]
  save(genotype_mat_not_in_repeat, file = paste0(file_path, 'GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_chr',i,'_dosage_MAF_01_not_in_repeat.RData'))

  # Next, for MAF 05
  genotype_mat = read.table(paste0(file_path, 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr',i,'_dosage_MAF_05.txt'), header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  locs = as.numeric(sapply(genotype_mat$X, function(x) {strsplit(x, '_')[[1]][2]}))
  repeat_subset = repeat_masker[repeat_masker$chr == paste0('chr', i),]
  not_in_repeat = sapply(locs, function(x) {sum(x >= repeat_subset$start) == sum(x > repeat_subset$end)})
  genotype_mat_not_in_repeat = genotype_mat[not_in_repeat,]
  rownames(genotype_mat_not_in_repeat) = genotype_mat_not_in_repeat$X
  genotype_mat_not_in_repeat = genotype_mat_not_in_repeat[,c(2:ncol(genotype_mat_not_in_repeat))]
  save(genotype_mat_not_in_repeat, file = paste0(file_path, 'GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_chr',i,'_dosage_MAF_05_not_in_repeat.RData'))
# }
