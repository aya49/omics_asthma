##################################################
#  matrix_eqtl_wrapper_preset_list.py
#
#  $proj/Scripts/eqtls/trans/gtex/matrix_eqtl_wrapper_preset_list.py
# 
#  This script sets up MatrixEQTL runs for a preset list of variant-gene pairs
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
# argv = [str(i) for i in range(10)]
# argv[1] = '/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/Final_trans_eQTL_list_0.5.txt'
# argv[2] = '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
# argv[3] = '_nonverlapping_certain_autosomes_normalized.txt'
# argv[4] = '/Output/trans-mapping/gtex/MatrixEQTL/seg_runs/'
# argv[5] = '/Output/joblogs/trans-mapping/gtex/MatrixEQTL/seg_runs/'
# argv[6] = 'trans_matrix_eqtl_seg_runs.R'
# argv[7] = '2.5e8'
# argv[8] = '1'
# argv[9] = 'all'

input_list = proj_dir + argv[1]
in_path = proj_dir + argv[2]
in_suffix = argv[3]
out_dir = proj_dir + argv[4]
joblog_dir = proj_dir + argv[5]
Rscript = argv[6]
cis_dist = argv[7]
num_split = argv[8]
tissue_list = argv[9:]

# Make the job log directories
if not os.path.exists(joblog_dir):
    os.makedirs(joblog_dir)

master_script = proj_dir + '/Scripts/eqtls/trans/gtex/batch/matrix_eqtl_wrapper_' + Rscript + '.sh'
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

if tissue_list[0] == 'all':
    matrices = glob.glob(in_path + '*' + in_suffix)
else:
    matrices = [in_path+x+in_suffix for x in tissue_list]

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
    for n in range(int(num_split)):
        cmd_r = r"""/usr/bin/Rscript %s/Scripts/eqtls/trans/gtex/%s \
%s %s %s %s %s %s %s
"""%(proj_dir, Rscript, input_list, matrix, str(n+1), cis_dist, tissue_name, in_path + 'covariates/', out_file)
        sbatchhandle.write(cmd_r)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))
