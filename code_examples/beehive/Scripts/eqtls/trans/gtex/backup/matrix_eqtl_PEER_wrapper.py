##################################################
#  matrix_eqtl_PEER_wrapper.py
#
#  $proj/Scripts/eqtls/trans/gtex/matrix_eqtl_PEER_wrapper.py
# 
#  This script runs genome-wide association test between genotypes and covariates (PEER factors and genotype PCs)
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
# argv = [str(i) for i in range(7)]
# argv[1] = '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
# argv[2] = '_Analysis.covariates.txt'
# argv[3] = '/Output/trans-mapping/gtex/MatrixEQTL/cov_association/'
# argv[4] = '/Output/joblogs/trans-mapping/gtex/MatrixEQTL/cov_association/'
# argv[5] = 'trans_matrix_eqtl_PEER.R'
# argv[6] = '2.5e8'

in_path = proj_dir + argv[1]
in_suffix = argv[2]
out_dir = proj_dir + argv[3]
joblog_dir = proj_dir + argv[4]
Rscript = argv[5]
cis_dist = argv[6]

# Make the job log directories
if not os.path.exists(joblog_dir):
    os.makedirs(joblog_dir)

master_script = proj_dir + '/Scripts/eqtls/trans/gtex/batch/matrix_eqtl_wrapper_PEER_' + Rscript + '.sh'
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
    sbatchfile = proj_dir + '/Scripts/eqtls/trans/gtex/batch/matrix_eqtl' + Rscript + '_' + tissue_name + '.slurm'
    job_outfile = 'trans_matrix_eqtl' + Rscript + '_' + tissue_name
    sbatchhandle=open(sbatchfile, 'w')
    cmd=r"""#!/bin/bash
#SBATCH -J meqtl_%s     # job name
#SBATCH --mem=16000           # 16 GB requested
#SBATCH -t 24:00:00           # 24-hour short jobs
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

umask 002
"""%(tissue_name, joblog_dir+job_outfile+'.err', joblog_dir+job_outfile+'.out')
    sbatchhandle.write(cmd)
    for k in range(22):
        cmd_r = r"""/usr/bin/Rscript %s/Scripts/eqtls/trans/gtex/%s \
%s %s %s %s %s
"""%(proj_dir, Rscript, matrix, str(k+1), cis_dist, tissue_name, out_file)
        sbatchhandle.write(cmd_r)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))