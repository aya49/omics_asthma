##################################################
#  v8_trans_matrix_eqtl_wrapper.py
#
#  $proj/Scripts/eqtls/trans/gtex/v8_consortium/v8_trans_matrix_eqtl_wrapper.py
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
# argv = [str(i) for i in range(9)]
# argv[1] = '/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
# argv[2] = '.v8.normalized_expression.bed.gz'
# argv[3] = '/Output/trans-mapping/gtex/MatrixEQTL/v8/all-by-all/'
# argv[4] = '/Output/joblogs/trans-mapping/gtex/MatrixEQTL/v8/'
# argv[5] = 'v8_trans_matrix_eqtl.R'
# argv[6] = '2.5e8'
# argv[7] = '20000'

in_path = proj_dir + argv[1]
in_suffix = argv[2]
out_dir = proj_dir + argv[3]
joblog_dir = proj_dir + argv[4]
Rscript = argv[5]
cis_dist = argv[6]
# How many SNPs to handle in a job?
partition_size = argv[7]

# Make the job log directories
if not os.path.exists(joblog_dir):
    os.makedirs(joblog_dir)

master_script = proj_dir + '/Scripts/eqtls/trans/gtex/v8_consortium/batch/matrix_eqtl_wrapper_' + Rscript + '.sh'
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

matrices = glob.glob(in_path + '*' + in_suffix)
print(len(matrices))

# custom tissue list for now
tissue_list = ['Adipose_Subcutaneous', 'Skin_Sun_Exposed_Lower_leg', 'Testis', 'Thyroid', 'Cells_EBV-transformed_lymphocytes']
variants_per_chr = [int(x.strip()) for x in open(proj_dir + '/Data/Genotype/gtex/v8/allelic_dosage/num_snps_per_chr.txt').readlines()]

num_parts_per_chr = [math.ceil(x / int(partition_size)) for x in variants_per_chr]
print(num_parts_per_chr)

# num_jobs = math.ceil(22 * int(num_split) / int(scripts_per_run))

for tissue in tissue_list:
    filename = in_path + tissue + in_suffix
    out_tissue_dir = out_dir + tissue + '/'
    summary_dir = out_tissue_dir + 'summary/'
    if not os.path.exists(out_tissue_dir):
        os.makedirs(out_tissue_dir)
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)
    for chr_num in range(22):
        num_parts = math.ceil(variants_per_chr[chr_num] / int(partition_size))
        # Iterate through 22*int(num_split) parts - we will have num_jobs jobs that process int(scripts_per_run) parts each. User can further customize this construction.
        for k in range(num_parts):
            # check that output exists
            if os.path.exists(summary_dir + 'chr' + str(chr_num + 1) + '_part' + str(k + 1) + '.RData'):
                continue
            sbatchfile = proj_dir + '/Scripts/eqtls/trans/gtex/v8_consortium/batch/matrix_eqtl' + Rscript + '_' + tissue + '_chr' + str(chr_num + 1) + '_part' + str(k + 1) + '.slurm'
            job_outfile = 'trans_matrix_eqtl' + Rscript + '_' + tissue + '_chr' + str(chr_num + 1) + '_part' + str(k + 1)
            sbatchhandle=open(sbatchfile, 'w')
            cmd=r"""#!/bin/bash
#SBATCH -J %s_%s_%s     # job name
#SBATCH --mem=24000           # 24 GB requested
#SBATCH -t 24:00:00           # 24-hour short jobs
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

umask 002

Rscript %s/Scripts/eqtls/trans/gtex/v8_consortium/%s %s %s %s %s %s %s %s %s

"""%(tissue, str(chr_num + 1), str(k + 1), joblog_dir+job_outfile+'.err', joblog_dir+job_outfile+'.out', proj_dir, Rscript, filename, str(k + 1), str(chr_num + 1), partition_size, cis_dist, tissue, in_path + 'GTEx_Analysis_v8_eQTL_covariates/', out_tissue_dir)
            sbatchhandle.write(cmd)
            sbatchhandle.close()
            master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))

# args = c(1:7)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz'
# args[2] = '1'
# args[3] = '2'
# args[4] = '100'
# args[5] = '2.5e8'
# args[6] = 'Whole_Blood'
# args[7] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
# args[8] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/v8/all-by-all/Whole_Blood/'
