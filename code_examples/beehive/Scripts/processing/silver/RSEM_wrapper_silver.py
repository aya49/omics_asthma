##################################################
#  RSEM_wrapper_silver.py
#
#  $proj/Scripts/processing/silver/RSEM_wrapper_silver.py
#
#  After the 2nd-pass of STAR, quantify isoform expression levels
#
#  Authors: Ian McDowell, Brian Jo
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

# with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/RSEM_hg38_samples.txt', 'r') as f:
#     samples_done = [line.strip() for line in f]

# Check for completed samples:
RSEM_dir = '/tigress/BEE/gtex/data/phenotype/expression/quantified_rna_seq_reads/silver/RSEM_hg38/'
directories = glob.glob(RSEM_dir + '*')
samples_done = [x for x in directories if str.split(x, '/')[-1]+'.temp' not in os.listdir(x + '/')]
samples_done = [str.split(x, '/')[-1] for x in samples_done]

samples_not_done = [x for x in directories if str.split(x, '/')[-1]+'.temp' in os.listdir(x + '/')]
for unfinished in samples_not_done:
    subprocess.call(['rm', '-r', unfinished])

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/RSEM_hg38_samples.txt', 'w') as f:
    [f.write(x + '\n') for x in samples_done]
    f.close()

samples_to_do = [sample for sample in samples_to_do if sample not in set(samples_done)]

samples_to_do_split = [samples_to_do[x:x+split] for x in range(0, len(samples_to_do), split)]
jobs_per_contrib = math.ceil(len(samples_to_do_split) / contributors)
print(jobs_per_contrib)

counter = 0
# job names are in numbers - can adjust to the list of contributors
for contributor in range(contributors):
    n = str(contributor + 1)
    master_script = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/RSEM_wrapper_" + n + ".sh"
    master_handle=open(master_script, 'w')
    master_handle.write("#!/bin/bash\n\n")
    for i in range(jobs_per_contrib):
        sbatchfile = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/RSEM_" + n + "_part_" + (str(i+1)) + ".sh"
        sbatchhandle=open(sbatchfile, 'w')
        header=r"""#!/bin/bash
#SBATCH -J %s_%s_RSEM      # job name
#SBATCH -o /tigress/BEE/gtex/results/group_general/joblogs/RSEM_hg38/%s_%s_RSEM.txt
#SBATCH --mem=24000
#SBATCH --ntasks-per-node=4 
#SBATCH --get-user-env
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

QUANT_DIR=/tigress/BEE/gtex/data/phenotype/expression/quantified_rna_seq_reads/silver/RSEM_hg38
MAPPED_DIR=/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/STAR_2pass_hg38
RSEM_REF=/tigress/BEE/gtex/data/silver/RSEM_hg38/rsem_hg38
"""%(n, str(i+1), n, str(i+1))
        sbatchhandle.write(header)
        if counter >= len(samples_to_do_split):
            break
        for sample in samples_to_do_split[counter]:
            sbatch="""

SAMPLE=%s

mkdir -m 775 -p ${QUANT_DIR}/${SAMPLE}
cd ${QUANT_DIR}/${SAMPLE}

/tigress/BEE/bin/RSEM-1.3.0/rsem-calculate-expression --quiet -p 4 \
--no-bam-output --bam --paired-end --seed 1234 --append-names \
${MAPPED_DIR}/${SAMPLE}/Aligned.toTranscriptome.out.bam \
$RSEM_REF $SAMPLE &> ${SAMPLE}.err.txt

"""%(sample)
            sbatchhandle.write(sbatch)
        sbatchhandle.close()
        counter = counter + 1
        master_handle.write("sbatch " + sbatchfile  + " \n")
    master_handle.close()
    print('sh ' + master_script)
