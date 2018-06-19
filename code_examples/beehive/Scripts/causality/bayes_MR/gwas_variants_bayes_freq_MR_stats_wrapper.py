##################################################
#  gwas_variants_bayes_freq_MR_stats_wrapper.py
#
#  $proj/Scripts/causality/bayes_MR/gwas_variants_bayes_freq_MR_stats_wrapper.py
# 
#  This script calculates the Wakefield Approximate Bayes Factors for the GTEx v8 data.
#
#  Author: Brian Jo
#
##################################################

import glob
import os
import os.path
import subprocess as sp
import math

proj_dir = os.environ['proj']

# Pre-set tissue list
tissue_list = ['Adipose_Subcutaneous', 'Whole_Blood', 'Lung', 'Skin_Sun_Exposed_Lower_leg', 'Thyroid']
tissue_list_gwas_file = ['adipose', 'blood', 'lung', 'skin', 'thyroid']

expression_input_dir = proj_dir + '/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
out_dir = proj_dir + '/Output/causality/gtex/bayes_MR/gwas/'
cov_dir = expression_input_dir + 'GTEx_Analysis_v8_eQTL_covariates/'

for tis in tissue_list:
  if not os.path.exists(out_dir + tis):
    os.mkdir(out_dir + tis)

master_script = proj_dir + '/Scripts/causality/bayes_MR/gwas_variants_bayes_freq_MR_stats.sh'
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

for i in range(5):
  tis = tissue_list[i]
  input_expr = expression_input_dir + tis + '.v8.normalized_expression.bed.gz'
  input_gwas_file = proj_dir + '/Data/Resources/gwas/gwas-association-downloaded-' + tissue_list_gwas_file[i] + '.tsv'
  out_file_prefix = out_dir + tis + '/gwas_bayes_freq_MR_stats'
  sbatchfile = proj_dir + '/Scripts/causality/bayes_MR/batch/gwas_variants_bayes_freq_MR_stats_' + tis + '.slurm'
  sbatchhandle=open(sbatchfile, 'w')
  cmd=r"""#!/bin/bash
#SBATCH -J gwas_MR_%s           # job name
#SBATCH --mem=32000             # 32 GB requested
#SBATCH -t 24:00:00
#SBATCH -e /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/gwas_variants_bayes_freq_MR_stats_%s.err           # err output directory
#SBATCH -o /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/gwas_variants_bayes_freq_MR_stats_%s.out           # out output directory        

Rscript /tigress/BEE/RNAseq/Scripts/causality/bayes_MR/gwas_variants_bayes_freq_MR_stats.R %s %s %s %s %s
"""%(tis, tis, tis, input_gwas_file, input_expr, tis, cov_dir, out_file_prefix)
  sbatchhandle.write(cmd)
  sbatchhandle.close()
  master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))

# args[1] = '/tigress/BEE/RNAseq/Data/Resources/gwas/gwas-association-downloaded-lung.tsv'
# args[2] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Lung.v8.normalized_expression.bed.gz'
# args[3] = 'Lung'
# args[4] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
# args[5] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/gwas/Lung/gwas_bayes_freq_MR_stats'
