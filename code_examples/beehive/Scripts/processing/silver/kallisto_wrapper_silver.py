##################################################
#  kallisto_wrapper_silver.py
#
#  $proj/Scripts/processing/silver/kallisto_wrapper_silver.py
#
#  For pseudo-alignment kallisto runs with bootstraps.
#
#  Author: Brian Jo
#
##################################################
from sys import argv
import math
import glob
import os
import subprocess

split = int(argv[1])
contributors = int(argv[2])
# A little less than an hour per sample, 24 hour job
# split = 15

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/v6p_samples.txt', 'r') as f:
    samples_to_do = [line.strip() for line in f]
    f.close()

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/kallisto_hg38_bootstrap_samples.txt', 'r') as f:
    bootstrap_samples_done = [line.strip() for line in f]
    f.close()

samples_to_do = [sample for sample in samples_to_do if sample not in set(bootstrap_samples_done)]
samples_to_do_split = [samples_to_do[x:x+split] for x in range(0, len(samples_to_do), split)]
jobs_per_contrib = math.ceil(len(samples_to_do_split) / contributors)
print(jobs_per_contrib)

counter = 0
# job names are in numbers - can adjust to the list of contributors
for contributor in range(contributors):
    n = str(contributor + 1)
    master_script = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/kallisto_wrapper_" + n + ".sh"
    master_handle=open(master_script, 'w')
    master_handle.write("#!/bin/bash\n\n")
    for i in range(jobs_per_contrib):
        sbatchfile = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/kallisto_" + n + "_part_" + (str(i+1)) + ".sh"
        sbatchhandle=open(sbatchfile, 'w')
        header=r"""#!/bin/bash
#SBATCH -J %s_%s_kal      # job name
#SBATCH -o /tigress/BEE/gtex/results/group_general/joblogs/kallisto_hg38/%s_%s_kal.txt
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

"""%(n, str(i+1), n, str(i+1))
        sbatchhandle.write(header)
        if counter >= len(samples_to_do_split):
            break
        for sample in samples_to_do_split[counter]:
            sbatch="""

SAMPLE=%s

/tigress/BEE/bin/kallisto quant -t 4 -i /tigress/BEE/gtex/data/silver/kallisto/GRCh38_combined_transcripts.idx \
-o /tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/kallisto_hg38/${SAMPLE} -b 100 \
<(bzip2 -dkc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/${SAMPLE}_1.trimmed.P.fastq.bz2) \
<(bzip2 -dkc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/${SAMPLE}_2.trimmed.P.fastq.bz2)

            """%(sample)
            sbatchhandle.write(sbatch)
        sbatchhandle.close()
        counter = counter + 1
        master_handle.write("sbatch " + sbatchfile  + " \n")
    master_handle.close()
    print('sh ' + master_script)
