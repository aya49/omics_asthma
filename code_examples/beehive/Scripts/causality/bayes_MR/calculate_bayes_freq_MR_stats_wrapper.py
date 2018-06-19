##################################################
#  calculate_bayes_freq_MR_stats_wrapper.py
#
#  $proj/Scripts/causality/bayes_MR/calculate_bayes_freq_MR_stats_wrapper.py
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

expression_input_dir = proj_dir + '/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
tissue_list = [x.split('/')[-1].split('.v8')[0] for x in glob.glob(expression_input_dir + '*' + '.v8.normalized_expression.bed.gz')]
wabf_raw_dir = proj_dir + '/Output/trans-mapping/gtex/WABF/raw/'
cov_dir = expression_input_dir + 'GTEx_Analysis_v8_eQTL_covariates/'
out_dir = proj_dir + '/Output/causality/gtex/bayes_MR/raw/'

for tis in tissue_list:
  if not os.path.exists(out_dir + tis):
    os.mkdir(out_dir + tis)

master_script = proj_dir + '/Scripts/causality/bayes_MR/calculate_bayes_freq_MR_stats.sh'
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

for tis in tissue_list:
  input_expr = expression_input_dir + tis + '.v8.normalized_expression.bed.gz'
  wabf_files = glob.glob(wabf_raw_dir + tis + '/' + '*')
  for f in wabf_files:
    file_prefix = f.split('wabf_raw_output_')[-1].split('.RData')[0]
    chr_number = file_prefix.split('_')[0][3:]
    out_file_prefix = out_dir + tis + '/bayes_freq_MR_stats_' + file_prefix + '_'
    sbatchfile = proj_dir + '/Scripts/causality/bayes_MR/batch/calculate_bayes_freq_MR_stats_' + tis + '_' + file_prefix + '.slurm'
    sbatchhandle=open(sbatchfile, 'w')
    cmd=r"""#!/bin/bash
#SBATCH -J %s      # job name
#SBATCH --mem=16000             # 16 GB requested
#SBATCH -t 24:00:00
#SBATCH -e /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/calculate_MR_stats_%s.err           # err output directory
#SBATCH -o /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/calculate_MR_stats_%s.out           # out output directory        

Rscript /tigress/BEE/RNAseq/Scripts/causality/bayes_MR/calculate_bayes_freq_MR_stats.R %s %s %s %s %s %s
"""%(file_prefix, file_prefix, file_prefix, input_expr, f, chr_number, tis, cov_dir, out_file_prefix)
    sbatchhandle.write(cmd)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))

# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz'
# args[2] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/WABF/raw/Whole_Blood/wabf_raw_output_chr10_part1_risk_1.5_pi1_1e-05.RData'
# args[3] = '10'
# args[4] = 'Whole_Blood'
# args[5] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
# args[6] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/raw/Whole_Blood/bayes_freq_MR_stats_chr10_part1_risk_1.5_pi1_1e-05_'