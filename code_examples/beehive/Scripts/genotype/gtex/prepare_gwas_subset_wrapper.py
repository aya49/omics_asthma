##################################################
#  prepare_gwas_subset_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/prepare_gwas_subset_wrapper.py
# 
#  This script saves a separate genotype file for the gwas-subsetting run, by taking the best cis SNP, in a tissue-specific fashion
#
#  Author: Brian Jo
#
##################################################

import glob
import os
import string
import os.path
from sys import argv

exclude = set(string.punctuation)

joblog_dir = '/tigress/BEE/RNAseq/Output/joblogs/misc/'

snpgene_input_repo = '/tigress/BEE/RNAseq/Data/References/GTEx/v6p_fastQTL_FOR_QC_ONLY/'
os.chdir(snpgene_input_repo)
snpgene_files = glob.glob('*' + '_Analysis.v6p.FOR_QC_ONLY.snpgenes.txt.gz')

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/prepare_gwas_subset.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

# 44 snpgene files, all in our mapping set
# exclude = set(string.punctuation)

# for item in [str.split(x, '_Analysis')[0] for x in snpgene_files]:
#       tissue = ''.join(ch for ch in item if ch not in exclude).lower()
#       print(tissue in tissue_list)

for files in snpgene_files:
        tissue = ''.join(ch for ch in str.split(files, '_Analysis')[0] if ch not in exclude).lower()
        sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/prepare_gwas_subset_" + tissue + '.slurm'
        sbatchhandle=open(sbatchfile, 'w')
        cmd=r"""#!/bin/bash
#SBATCH -J %s      # job name
#SBATCH --mem=24000     # 24 GB requested
#SBATCH -t 8:00:00     
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

/usr/bin/Rscript /tigress/BEE/RNAseq/Scripts/genotype/gtex/prepare_gwas_subset.R %s
"""%(tissue, joblog_dir+'err/prepare_gwas_subset_'+tissue+'.err', joblog_dir+'out/prepare_gwas_subset_'+tissue+'.out', tissue)
        sbatchhandle.write(cmd)
        sbatchhandle.close()
        master_handle.write("sbatch " + sbatchfile  + " \n")

print('sh ' + master_script)
master_handle.close()