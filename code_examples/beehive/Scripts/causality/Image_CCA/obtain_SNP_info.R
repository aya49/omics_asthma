##################################################
#  obtain_SNP_info.R
#
#  $proj/Scripts/causality/Image_CCA/obtain_SNP_info.R
# 
#  Get the SNP information from the genotype file, and split them into N chunks for the MR pipeline.
#
#  Author: Brian Jo
#
##################################################

genotype_file_input = '/tigress/gdarnell/imagecca/MR_data/genotypesChrAll_snps_and_indv_of_interest.txt'
genotype_df = read.table(genotype_file_input, header=T, stringsAsFactors=F)

SNP_position_dir = '/tigress/BEE/RNAseq/Data/Genotype/gtex/SNP_positions_hg19/'