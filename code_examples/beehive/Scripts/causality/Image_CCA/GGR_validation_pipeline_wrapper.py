##################################################
#  GGR_validation_pipeline_wrapper.py
#
#  $proj/Scripts/causality/GGR/gtex/GGR_validation_pipeline_wrapper.py
# 
#  GGR validation pipeline with trans-eQTLs and MR implementations.
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
# argv = [str(i) for i in range(11)]
# argv[1] = '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
# argv[2] = '_v6p_consortium_autosomes_normalized.txt'
# argv[3] = '/Data/Networks/GGR/retrofitted/'
# argv[4] = '/Output/causality/GGR/output/'
# argv[5] = '/Output/joblogs/causality/GGR/'
# argv[6] = '0mean-1var_g-null_l-fdr'
# argv[7] = '2.5e8'
# argv[8] = '1000'
# argv[9] = 'lasso-2'
# argv[10] = 'lasso-1'

in_path = proj_dir + argv[1]
in_suffix = argv[2]
network_dir = proj_dir + argv[3]
output_dir = proj_dir + argv[4]
joblog_dir = proj_dir + argv[5]
mode = argv[6]
cis_dist = argv[7]
pairs_per_script = int(argv[8])
regression_list = argv[9:]

# list of tissues to process
tissues = ['lung', 'adiposesubcutaneous', 'cellstransformedfibroblasts', 'arterytibial', 'thyroid']

# Make the job log directories
if not os.path.exists(joblog_dir):
    os.makedirs(joblog_dir)

if not os.path.exists(output_dir + mode + '/'):
    os.makedirs(output_dir + mode + '/')

master_script = proj_dir + '/Scripts/causality/GGR/gtex/batch/GGR_validation_wrapper_' + mode + '_' + '_'.join(regression_list) + '.sh'
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

files = [glob.glob(network_dir + mode + '/' + '*' + x + '*') for x in regression_list]
files = [item for sublist in files for item in sublist]
print(len(files))

for f in files:
    num_lines = sum(1 for line in open(f)) - 1
    num_parts = math.ceil(num_lines/pairs_per_script)
    filename = f.split('/')[-1]
    for k in range(num_parts):
        sbatchfile = proj_dir + '/Scripts/causality/GGR/gtex/batch/GGR_validation_' + filename + '_part' + str(k+1) + '.slurm'
        job_outfile = 'GGR_validation_' + filename + '_part' + str(k+1)
        sbatchhandle=open(sbatchfile, 'w')
        cmd=r"""#!/bin/bash
#SBATCH -J GGR_%s     # job name
#SBATCH --mem=24000           # 24 GB requested
#SBATCH -t 24:00:00           # 24-hour short jobs
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

umask 002
"""%(str(k+1), joblog_dir+job_outfile+'.err', joblog_dir+job_outfile+'.out')
        sbatchhandle.write(cmd)
        for tis in tissues:
            matrix = in_path + tis + in_suffix
            cmd_r = r"""/usr/bin/Rscript %s/Scripts/causality/GGR/gtex/GGR_validation_pipeline_redux.R \
%s %s %s %s %s %s %s %s
"""%(proj_dir, matrix, str(k+1), cis_dist, tis, in_path + 'covariates/', output_dir + mode + '/', f, str(num_parts))
            sbatchhandle.write(cmd_r)
        sbatchhandle.close()
        master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()
print('sh %s'%(master_script))
