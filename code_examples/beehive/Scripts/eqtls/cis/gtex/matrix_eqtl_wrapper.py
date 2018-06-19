##################################################
#  matrix_eqtl_wrapper.py
#
#  $proj/Scripts/eqtls/cis/gtex/matrix_eqtl_wrapper.py
# 
#  This version is the most up-to-date version for cis- pipeline.
#
#  Author: Brian Jo
#
##################################################

import glob
from sys import argv
import os.path
import os

proj_dir = os.environ['proj']

# Examples
# argv = [str(i) for i in range(8)]
# argv[1] = '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
# argv[2] = '_nonverlapping_certain_autosomes_normalized.txt'
# argv[3] = '/Output/cis-mapping/gtex/MatrixEQTL/gct-normalized-1M/'
# argv[4] = '/Output/joblogs/cis-mapping/gtex/MatrixEQTL/gct-normalized-1M/'
# argv[5] = 'cis_matrix_eqtl.R'
# argv[6] = '1000000'

in_path = proj_dir + argv[1]
in_suffix = argv[2]
out_dir = proj_dir + argv[3]
joblog_dir = proj_dir + argv[4]
Rscript = argv[5]
cis_dist = argv[6]

# Make the job log directories
if not os.path.exists(joblog_dir):
    os.makedirs(joblog_dir)

master_script = proj_dir + '/Scripts/eqtls/cis/gtex/batch/matrix_eqtl_wrapper_' + Rscript + '.sh'
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

matrices = glob.glob(in_path + '*' + in_suffix)
print(len(matrices))

for i, matrix in enumerate(matrices):
    filename = matrix.split('/')[-1]
    tissue_name = str.split(filename, in_suffix)[0]
    out_tissue_dir = out_dir + tissue_name + '/'
    if not os.path.exists(out_tissue_dir):
        os.makedirs(out_tissue_dir)
    out_file = out_tissue_dir + filename.replace('.txt','_MatrixEQTL')
    # Create one batch job for each tissue
    sbatchfile = proj_dir + '/Scripts/eqtls/cis/gtex/batch/matrix_eqtl' + Rscript + '_' + tissue_name + '.slurm'
    job_outfile = 'cis_matrix_eqtl' + Rscript + '_' + tissue_name
    sbatchhandle=open(sbatchfile, 'w')
    cmd=r"""#!/bin/bash
#SBATCH -J meqtl_%s     # job name
#SBATCH --mem=32000           # 32 GB requested
#SBATCH -t 24:00:00           # 24-hour short jobs
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

umask 002
"""%(tissue_name, joblog_dir+job_outfile+'.err', joblog_dir+job_outfile+'.out')
    sbatchhandle.write(cmd)
    for n in range(22):
        cmd_r = r"""/usr/bin/Rscript %s/Scripts/eqtls/cis/gtex/%s \
%s %s %s %s %s %s
"""%(proj_dir, Rscript, matrix, str(n+1), cis_dist, tissue_name, in_path + 'covariates/', out_file)
        sbatchhandle.write(cmd_r)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))
