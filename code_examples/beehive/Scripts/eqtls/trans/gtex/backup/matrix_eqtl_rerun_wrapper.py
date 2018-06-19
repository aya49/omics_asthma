##################################################
#  matrix_eqtl_rerun_wrapper.py
#
#  $proj/Scripts/eqtls/trans/gtex/matrix_eqtl_rerun_wrapper.py
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

import glob
from sys import argv
import os.path
import os
import math

proj_dir = os.environ['proj']

# Examples
argv = [str(i) for i in range(9)]
argv[1] = '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
argv[2] = '_nonverlapping_certain_autosomes_normalized.txt'
argv[3] = '/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/'
argv[4] = '/Output/joblogs/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/'
argv[5] = 'trans_matrix_eqtl_PEER_increments.R'
argv[6] = '2.5e8'
argv[7] = '20'
argv[8] = '30'

in_path = proj_dir + argv[1]
in_suffix = argv[2]
out_dir = proj_dir + argv[3]
joblog_dir = proj_dir + argv[4]
Rscript = argv[5]
cis_dist = argv[6]
# Split up each chromosome into how many parts?
num_split = argv[7]
# How many scripts to run in a job?
scripts_per_run = argv[8]

# Make the job log directories
if not os.path.exists(joblog_dir):
    os.makedirs(joblog_dir)

master_script = proj_dir + '/Scripts/eqtls/trans/gtex/batch/matrix_eqtl_wrapper_' + Rscript + '.sh'
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

matrices = glob.glob(in_path + '*' + in_suffix)
print(len(matrices))

num_jobs = math.ceil(22 * int(num_split) / int(scripts_per_run))

for i, matrix in enumerate(matrices):
    filename = matrix.split('/')[-1]
    tissue_name = str.split(filename, in_suffix)[0]
    out_tissue_dir = out_dir + tissue_name + '/'
    if not os.path.exists(out_tissue_dir):
        os.makedirs(out_tissue_dir)
    out_file = out_tissue_dir + filename.replace('.txt','_MatrixEQTL')
    # Check if the tissue ran successfully with all files:
    if len(glob.glob(out_file + '*' + '_summary.txt')) == 22 * int(num_split):
        continue
    # Iterate through 22*int(num_split) parts - we will have num_jobs jobs that process int(scripts_per_run) parts each. User can further customize this construction.
    for k in range(num_jobs):
        sbatchfile = proj_dir + '/Scripts/eqtls/trans/gtex/batch/matrix_eqtl' + Rscript + '_' + tissue_name + '_part' + str(k+1) + '.slurm'
        job_outfile = 'trans_matrix_eqtl' + Rscript + '_' + tissue_name + '_part' + str(k+1)
        sbatchhandle=open(sbatchfile, 'w')
        cmd=r"""#!/bin/bash
#SBATCH -J meqtl_%s_%s     # job name
#SBATCH --mem=24000           # 24 GB requested
#SBATCH -t 24:00:00           # 24-hour short jobs
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

umask 002
"""%(tissue_name, str(k+1), joblog_dir+job_outfile+'.err', joblog_dir+job_outfile+'.out')
        sbatchhandle.write(cmd)
        for n in range(int(scripts_per_run)):
            j = k*int(scripts_per_run) + n + 1
            # This file temporarily created for the PEER increment runs : the only difference is to check for the existence of summary file in this case
            if not os.path.isfile(out_file + '_part' + ("%03d" % j) + '_summary.txt'):
                cmd_r = r"""/usr/bin/Rscript %s/Scripts/eqtls/trans/gtex/%s \
%s %s %s %s %s %s %s
"""%(proj_dir, Rscript, matrix, str(j), cis_dist, tissue_name, in_path + 'covariates/', out_file, num_split)
                sbatchhandle.write(cmd_r)
            if j == 22 * int(num_split):
                break
        sbatchhandle.close()
        master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))
