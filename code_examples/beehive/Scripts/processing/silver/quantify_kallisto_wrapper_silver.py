##################################################
#  quantify_kallisto_silver_wrapper.py
#
#  $proj/Scripts/processing/silver/quantify_kallisto_wrapper_silver.py
#
#  Create expression matrix for kallisto
#
#  Authors: Brian Jo
#
##################################################
import pandas as pd

tables_dir = '/tigress/BEE/RNAseq/Data/Resources/gtex/tables/'
tissue_table = pd.read_csv(tables_dir + 'tissue_table.txt', sep='\t')

master_script = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/quantify_kallisto_silver.sh"
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

for tissue in tissue_table.index:
	if tissue_table.loc[[tissue]]['num_samples_with_geno'][0] >= 10:
		sbatchfile = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/quantify_kallisto_silver_" + tissue + ".sh"
		sbatchhandle=open(sbatchfile, 'w')
		header=r"""#!/bin/bash
#SBATCH -J %s_kal_quant      # job name
#SBATCH -o /tigress/BEE/gtex/results/group_general/joblogs/quant/%s_kallisto
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

python /tigress/BEE/RNAseq/Scripts/processing/silver/quantify_kallisto_silver.py %s

"""%(tissue, tissue, tissue)
		sbatchhandle.write(header)
		sbatchhandle.close()
		master_handle.write("sbatch " + sbatchfile + '\n')

print("sh " + master_script)
master_handle.close()