##################################################
#  extract_gtex_gct_rpkm_expression.py
#
#  $proj/Scripts/misc/gtex/extract_gtex_gct_rpkm_expression.py
# 
#  This script extracts the GTEx gct expression data in a usable format.
#
#  Author: Brian Jo
#
##################################################

# TODO: complete this portion
# Need to have processed the gtf files
# Read in the expression files, gene locations, list of subjects with genotypes and nonoverlapping gene list

import os
import glob
import numpy as np
import pandas as pd
import string

proj_dir = os.environ['proj']
subject_table = pd.read_csv(proj_dir + '/Data/Resources/gtex/tables/subject_table.txt', sep = '\t')
# Get the list of subjects with genotypes
genotype_subjects = subject_table.index[subject_table['genotype_avail']]

nonoverlap_gene_list = pd.read_csv(proj_dir + '/Data/Resources/gtex/annotations/gencode.v19.nonoverlap.genes_for_quantification_certain_and_uncertain.txt', header=None)
nonoverlap_gene_list = np.array(nonoverlap_gene_list[0])
nonoverlap_gene_list_certain = pd.read_csv(proj_dir + '/Data/Resources/gtex/annotations/gencode.v19.nonoverlap.genes_for_quantification_certain.txt', header=None)
nonoverlap_gene_list_certain = np.array(nonoverlap_gene_list_certain[0])

in_dir = proj_dir + '/Data/Resources/gtex/dbGaP/GTEx_phs000424/v6p_fastQTL_FOR_QC_ONLY/'
bedfiles = glob.glob(in_dir + '*normalized.expr.bed')

out_dir = proj_dir + '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'

for files in bedfiles:
	# String translation, the python 3 way:
	tissue_name = files.split('/')[-1].split('_Analysis')[0]
	tissue_name = tissue_name.translate(str.maketrans({key: None for key in string.punctuation})).lower()
	print(tissue_name)
	# read in expression data
	tissue_df = pd.read_csv(files, sep = '\t', compression = 'gzip')
	tissue_df.index = tissue_df['gene_id']
	# Divide up into X chromosome and autosomes
	tissue_df_xchr = tissue_df[tissue_df['#chr'] == 'X']
	tissue_df_auto = tissue_df[tissue_df['#chr'] != 'X']
	# Take only columns that have genotypes
	tissue_subjects = [x for x in genotype_subjects if x in tissue_df.columns]
	tissue_df_xchr = tissue_df_xchr[tissue_subjects]
	tissue_df_auto = tissue_df_auto[tissue_subjects]
	# Take only from the list of nonoverlapping genes
	tissue_df_xchr_all = tissue_df_xchr.loc[[x for x in nonoverlap_gene_list if x in tissue_df_xchr.index]]
	tissue_df_xchr_certain = tissue_df_xchr.loc[[x for x in nonoverlap_gene_list_certain if x in tissue_df_xchr.index]]
	tissue_df_auto_all = tissue_df_auto.loc[[x for x in nonoverlap_gene_list if x in tissue_df_auto.index]]
	tissue_df_auto_certain = tissue_df_auto.loc[[x for x in nonoverlap_gene_list_certain if x in tissue_df_auto.index]]
	# output files
	tissue_df_xchr_all.to_csv(out_dir + tissue_name + '_nonverlapping_all_xchr_normalized.txt', sep = '\t')
	tissue_df_xchr_certain.to_csv(out_dir + tissue_name + '_nonverlapping_certain_xchr_normalized.txt', sep = '\t')
	tissue_df_auto_all.to_csv(out_dir + tissue_name + '_nonverlapping_all_autosomes_normalized.txt', sep = '\t')
	tissue_df_auto_certain.to_csv(out_dir + tissue_name + '_nonverlapping_certain_autosomes_normalized.txt', sep = '\t')

