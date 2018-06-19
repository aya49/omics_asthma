
# quantification_script_wrapper.py
# This script prepares the kallisto quantifiation pipeline, with 100 bootstrap samples, distributed over jobs
# Author: Brian Jo
# Date: 10/28/16

import glob
from sys import argv
import os.path
from pandas import DataFrame
import numpy as np
import string
import random

# Fixed params

# Location of trimmed read source
trimmed_reads_source = '/tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/'
# The suffix of read files in each sample read directory
read_file_suffix = '_2.trimmed.P.fastq.bz2'
# Location of samples table:
samp_table = "/tigress/BEE/RNAseq/Data/Resources/gtex/tables/sample_table.txt"
# Location of kallisto quantification output directory
out_dir = '/tigress/BEE/RNAseq/Output/kallisto/'
# Location of batch script output
batch_script_dir = '/tigress/BEE/RNAseq/Scripts/kallisto/batch/'

# Tweakable params

# Number of samples to include per run can be customized
num_samp_per_run = int(argv[1])
# Number of people running the job (or number of sessions to split the runs over)
num_session_per_run = int(argv[2])
# Memory in gigabytes
mem = argv[3]
# List of tissues to process
tissue_list = argv[4:]

# Example param set

# num_samp_per_run = 3
# num_session_per_run = 1
# mem = 50
# tissue_list = ['pancreas']

# Read in the list of samples that correspond to the tissue list:
sample_df = DataFrame.from_csv(samp_table, sep="\t")

sample_list = []

for tissue in tissue_list:
	sample_list = np.append(sample_list, DataFrame.as_matrix(sample_df['Run_s'][sample_df['tissue_name'] == tissue]))

# Exclude samples that have already been processed:
completed_list = []
already_processed = glob.glob(out_dir + 'cdna/' + '*')
for directory in already_processed:
	# Tentative: just check for existence of abundance.h5
	if os.path.isfile(directory + '/abundance.h5'):
	# if os.path.isdir(directory + '/bootstrap/'):
	# 	# All bootstrap samples completed?
	# 	if len(glob.glob(directory + '/bootstrap/bs_*')) == 100:
		completed_list.append(str.split(directory, '/')[-1])

sample_list = [x for x in sample_list if x not in completed_list]

# Now sample_list has the list of all samples to process - split them up into appropriate jobs
print('Samples to process: ' + str(len(sample_list)))
sample_partition = []
# Number of partitions
for i in range(num_session_per_run):
	partial_list = [sample_list[x] for x in range(len(sample_list)) if x % num_session_per_run == i]
	# Split up the list into jobs that contain the number of samples indicated:
	session_jobs = []
	for j in range(int(len(partial_list)/num_samp_per_run) + 1):
		session_jobs.append(partial_list[num_samp_per_run*j : num_samp_per_run*(j+1)])
	# Append to the master list
	sample_partition.append(session_jobs)

# Now put the partitions into processing scripts:
for i in range(len(sample_partition)):
	rand_key = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
	# Sessions split up
	master_script = batch_script_dir + 'ka_quant_' + rand_key + '.sh'
	master_handle = open(master_script, 'w')
	# Each job assigned the number of samples specified in num_samp_per_run
	for j in range(len(sample_partition[i])):
		sbatch_script = batch_script_dir + 'ka_quant_' + rand_key + 'part_' + str(j+1) + '.sh'
		cmd="""#!/bin/bash
#SBATCH -J ka_%s     # job name
#SBATCH --mem=%s000           # NN GB requested
#SBATCH -t 24:00:00           # 24-hour short jobs
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory
    
umask 002

"""%(str(j+1), mem, batch_script_dir + 'joblog/ka_quant_' + rand_key + '_' + str(j+1) + '.err', batch_script_dir + 'joblog/ka_quant_' + rand_key + '_' + str(j+1) + '.out')
		sbatch_script = batch_script_dir + 'ka_quant_' + rand_key + '_part_' + str(j+1) + '.slurm'
		sbatch_handle = open(sbatch_script, 'w')
		sbatch_handle.write(cmd)
		# Note: still optimizing: start with number of threads 4
		for sample in sample_partition[i][j]:
			cmd = """echo Starting %s

/tigress/BEE/tools_group/bin/kallisto quant -t 4 -i /tigress/BEE/RNAseq/Data/kallisto/GRCh37.p13.transcriptome_ncrna.idx \
-o /tigress/BEE/RNAseq/Output/kallisto/cdna_ncrna/%s -b 100 \
<(bzip2 -dkc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/%s_1.trimmed.P.fastq.bz2) \
<(bzip2 -dkc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/%s_2.trimmed.P.fastq.bz2)

echo Finished %s

"""%(sample, sample, sample, sample, sample)
			sbatch_handle.write(cmd)
			# Also process the version without ncrna (comment out if only doing the cdna_ncrna version)
			cmd = """echo Starting %s

/tigress/BEE/tools_group/bin/kallisto quant -t 4 -i /tigress/BEE/RNAseq/Data/kallisto/GRCh37.p13.transcriptome.idx \
-o /tigress/BEE/RNAseq/Output/kallisto/cdna/%s -b 100 \
<(bzip2 -dkc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/%s_1.trimmed.P.fastq.bz2) \
<(bzip2 -dkc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/%s_2.trimmed.P.fastq.bz2)

echo Finished %s

"""%(sample, sample, sample, sample, sample)
			sbatch_handle.write(cmd)
		sbatch_handle.close()
		master_handle.write('sbatch ' + sbatch_script + '\n')
	master_handle.close()
	print('sh ' + master_script)

print("Don't forget to change the permissions for the master scripts before distributing to other users!")
