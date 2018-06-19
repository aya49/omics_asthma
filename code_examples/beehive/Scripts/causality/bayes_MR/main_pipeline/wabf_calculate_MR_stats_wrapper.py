##################################################
#  wabf_calculate_MR_stats_wrapper.py
#
#  $proj/Scripts/causality/bayes_MR/main_pipeline/wabf_calculate_MR_stats_wrapper.py
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
out_dir = proj_dir + '/Output/causality/gtex/bayes_MR/MR_stats/'
prep_dir = proj_dir + '/Output/causality/gtex/bayes_MR/prep_files/'
n_per_run = 20000

master_script = proj_dir + '/Scripts/causality/bayes_MR/batch/wabf_calculate_MR_stats.sh'
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

tissue_list = ['Muscle_Skeletal', 'Adipose_Subcutaneous', 'Whole_Blood', 'Thyroid']

for tis in tissue_list:
  if not os.path.exists(out_dir + tis):
    os.mkdir(out_dir + tis)

for tis in tissue_list:
  out_file_prefix = out_dir + tis + '/wabf_mrstats_'
  eqtl_dir = proj_dir + '/Output/causality/gtex/bayes_MR/eqtls/' + tis + '/'
  for i in [str(x+1) for x in range(22)] + ['X']:
    sbatchfile = proj_dir + '/Scripts/causality/bayes_MR/batch/wabf_calculate_MR_stats_' + tis + '_chr' + i + '.slurm'
    sbatchhandle=open(sbatchfile, 'w')
    cmd=r"""#!/bin/bash
#SBATCH -J wabf_%s_%s      # job name
#SBATCH --mem=24000             # 24 GB requested
#SBATCH -t 4:00:00
#SBATCH -e /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/wabf_calculate_MR_stats_%s_%s.err           # err output directory
#SBATCH -o /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/wabf_calculate_MR_stats_%s_%s.out           # out output directory        

Rscript /tigress/BEE/RNAseq/Scripts/causality/bayes_MR/main_pipeline/wabf_calculate_MR_stats.R %s %s %s %s %s
"""%(tis, i, tis, i, tis, i, tis, i, prep_dir, eqtl_dir, out_file_prefix)
    sbatchhandle.write(cmd)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))

# args[1] = 'Muscle_Skeletal'
# args[2] = '10'
# args[3] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/prep_files/'
# args[4] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/eqtls/Muscle_Skeletal/'
# args[5] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/MR_stats/Muscle_Skeletal/wabf_mrstats_'
