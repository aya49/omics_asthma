##################################################
#  v8_normalize_matrices_wrapper.py
#
#  $proj/Scripts/processing/silver/v8_normalize_matrices_wrapper.py
#
#  Matrix normalization script
#
#  Authors: Brian Jo
#
##################################################
import pandas as pd

tables_dir = '/tigress/BEE/RNAseq/Data/Resources/gtex/tables/v8/'
tissue_table = pd.read_csv(tables_dir + 'tissue_table.txt', sep='\t')

master_script = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/v8_normalize_matrices.sh"
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

tissue_table.index = tissue_table['tissue_name']

for tissue in tissue_table.index:
	if tissue_table.loc[tissue]['num_samples'] > 10:
		sbatchfile = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/v8_normalize_matrices_" + tissue + ".sh"
		sbatchhandle=open(sbatchfile, 'w')
		header=r"""#!/bin/bash
#SBATCH -J %s_norm      # job name
#SBATCH -o /tigress/BEE/gtex/results/group_general/joblogs/quant/%s_normalize
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

python /tigress/BEE/RNAseq/Scripts/processing/silver/v8_normalize_matrices.py %s

"""%(tissue, tissue, tissue)
		sbatchhandle.write(header)
		sbatchhandle.close()
		master_handle.write("sbatch " + sbatchfile + '\n')

print("sh " + master_script)
master_handle.close()
