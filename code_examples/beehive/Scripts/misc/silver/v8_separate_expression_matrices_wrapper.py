##################################################
#  v8_separate_expression_matrices_wrapper.py
#
#  $proj/Scripts/misc/silver/v8_separate_expression_matrices_wrapper.py
# 
#  This script separates the GTEx expression matrices by tissues
#
#  Author: Brian Jo
#
##################################################

import os

proj_dir = os.environ['proj']

tissue_table = proj_dir + '/Data/Resources/gtex/tables/v8/tissue_table.txt'
f = open(tissue_table)
f.readline()

master_script = "/tigress/BEE/RNAseq/Scripts/misc/silver/batch/separate_expression_matrix.sh"
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

for line in f.readlines():
	tissue = str.split(line.strip(), '\t')[0]
	sbatchfile = "/tigress/BEE/RNAseq/Scripts/misc/silver/batch/separate_expression_matrix_" + tissue + ".sh"
	sbatchhandle=open(sbatchfile, 'w')
	script=r"""#!/bin/bash
#SBATCH -J %s      			# job name
#SBATCH -o /tigress/BEE/gtex/results/group_general/joblogs/quant/%s_v8
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

python /tigress/BEE/RNAseq/Scripts/misc/silver/v8_separate_expression_matrices.py %s

"""%(tissue, tissue, tissue)
	sbatchhandle.write(script)
	sbatchhandle.close()
	master_handle.write("sbatch " + sbatchfile + '\n')

print("sh " + master_script)
master_handle.close()