##################################################
#  v8_sfamix_incremental_wrapper.py
#
#  $proj/Scripts/processing/covariates_exploratory/v8_sfamix_incremental_wrapper.py
#
#  Trying out various normalization methods/factor analysis in v8 expression data, without the known covariates
#
#  Authors: Brian Jo
#
##################################################
from sys import argv
import glob
import os

proj_dir = os.environ['proj']
n_factors = 500
n_itr = 500
out_itr = 50

master_script = proj_dir + "/Scripts/processing/covariates_exploratory/batch/v8_sfamix_incremental.sh"
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

tissue_list = glob.glob(proj_dir + '/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/' + '*' + '.v8.normalized_expression.bed.gz')
tissue_list = [x.split('/')[-1] for x in tissue_list]
tissue_list = [x.split('.')[0] for x in tissue_list]

# Script for SFAmix
out_dir = proj_dir + '/Output/processing/exploratory/v8/sfamix/GTEx_Analysis_v8/'
exp_dir = proj_dir + '/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'

for tissue in tissue_list:
	directory = out_dir + tissue + '/'
	if not os.path.exists(directory):
		os.makedirs(directory)
	sbatchfile = proj_dir + '/Scripts/processing/covariates_exploratory/batch/v8_sfamix_incremental_' + tissue + '_' + str(n_factors) + '.sh'
	sbatchhandle=open(sbatchfile, 'w')
	script=r"""#!/bin/bash
#SBATCH -J %s_sfamix      # job name
#SBATCH -o %s/Scripts/processing/covariates_exploratory/joblogs/%s_%s_sfamix
#SBATCH --get-user-env
#SBATCH --mem=50000           # 50 GB requested
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

Rscript %s/Scripts/processing/covariates_exploratory/v8_sfamix_incremental.R %s %s %s %s %s %s

"""%(tissue, proj_dir, tissue, str(n_factors), proj_dir, exp_dir, directory, tissue, str(n_factors), str(n_itr), str(out_itr))
	sbatchhandle.write(script)
	sbatchhandle.close()
	master_handle.write("sbatch " + sbatchfile + '\n')

print("sh " + master_script)
master_handle.close()
