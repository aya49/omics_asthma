##################################################
#  wabf_calculate_bayes_factor_wrapper.py
#
#  $proj/Scripts/causality/bayes_MR/main_pipeline/wabf_calculate_bayes_factor_wrapper.py
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
out_dir = proj_dir + '/Output/causality/gtex/bayes_MR/eqtls/'
prep_dir = proj_dir + '/Output/causality/gtex/bayes_MR/prep_files/'
n_per_run = 20000

master_script = proj_dir + '/Scripts/causality/bayes_MR/batch/wabf_calculate_bayes_factor.sh'
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

tissue_list = ['Muscle_Skeletal', 'Adipose_Subcutaneous', 'Whole_Blood', 'Lung', 'Skin_Sun_Exposed_Lower_leg', 'Thyroid']

for tis in tissue_list:
  if not os.path.exists(out_dir + tis):
    os.mkdir(out_dir + tis)

for tis in tissue_list:
  out_file_prefix = out_dir + tis + '/wabf_eqtls_'
  for i in [str(x+1) for x in range(22)] + ['X']:
    z = sp.run(['wc', '-l', proj_dir + '/Data/Genotype/gtex/v8/ld_prune/chr' + i + '_MAF_05.in'], stdout = sp.PIPE)
    n_lines = int(z.stdout.decode("utf-8").split(' ')[0])
    z = sp.run(['wc', '-l', proj_dir + '/Data/Genotype/gtex/v8/ld_prune/chr' + i + '_MAF_01.in'], stdout = sp.PIPE)
    n_lines = n_lines + int(z.stdout.decode("utf-8").split(' ')[0])
    n_partitions = math.ceil(n_lines / n_per_run)
    for j in range(n_partitions):
      sbatchfile = proj_dir + '/Scripts/causality/bayes_MR/batch/wabf_calculate_bayes_factor_' + tis + '_chr' + i + '_part' + str(j+1) + '.slurm'
      sbatchhandle=open(sbatchfile, 'w')
      cmd=r"""#!/bin/bash
#SBATCH -J wabf_%s_%s_%s      # job name
#SBATCH --mem=24000             # 24 GB requested
#SBATCH -t 24:00:00
#SBATCH -e /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/wabf_calculate_bayes_factor_%s_%s_%s.err           # err output directory
#SBATCH -o /tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/joblogs/wabf_calculate_bayes_factor_%s_%s_%s.out           # out output directory        

Rscript /tigress/BEE/RNAseq/Scripts/causality/bayes_MR/main_pipeline/wabf_calculate_bayes_factor.R %s %s %s %s %s %s
"""%(tis, i, str(j+1), tis, i, str(j+1), tis, i, str(j+1), tis, i, str(j+1), str(n_per_run), prep_dir, out_file_prefix)
      sbatchhandle.write(cmd)
      sbatchhandle.close()
      master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))

# args[1] = 'Muscle_Skeletal'
# args[2] = '10'
# args[3] = '1'
# args[4] = '20000'
# args[5] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/prep_files/'
# args[6] = '/tigress/BEE/RNAseq/Output/causality/gtex/bayes_MR/eqtls/Muscle_Skeletal/wabf_eqtls_'
