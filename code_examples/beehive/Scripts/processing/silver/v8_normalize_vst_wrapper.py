##################################################
#  v8_normalize_vst_wrapper.py
#
#  $proj/Scripts/processing/silver/v8_normalize_vst_wrapper.py
# 
#  Matrix normalization script using vst in DESeq
#
#  Author: Brian Jo
#
##################################################
import glob
import os

proj_dir = os.environ['proj']

expression_input_dir = proj_dir + '/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
tissue_list = [x.split('/')[-1].split('.v8')[0] for x in glob.glob(expression_input_dir + '*' + '.v8.normalized_expression.bed.gz')]

master_script = proj_dir + '/Scripts/processing/silver/batch/v8_normalize_vst.sh'
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

for tis in tissue_list:
	sbatchfile = proj_dir + '/Scripts/processing/silver/batch/v8_normalize_vst_' + tis + '.slurm'
	sbatchhandle=open(sbatchfile, 'w')
	cmd=r"""#!/bin/bash
#SBATCH -J vst_%s      # job name
#SBATCH --mem=30000             # 16 GB requested
#SBATCH -t 24:00:00
#SBATCH -e /tigress/BEE/RNAseq/Output/joblogs/processing/v8_normalize_vst_%s.err           # err output directory
#SBATCH -o /tigress/BEE/RNAseq/Output/joblogs/processing/v8_normalize_vst_%s.out           # out output directory        

Rscript /tigress/BEE/RNAseq/Scripts/processing/silver/v8_normalize_vst.R %s
"""%(tis, tis, tis, tis)
	sbatchhandle.write(cmd)
	sbatchhandle.close()
	master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()
print('sh %s'%(master_script))

