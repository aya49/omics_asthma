##################################################
#  STAR_1pass_wrapper_silver.py
#
#  $proj/Scripts/processing/silver/STAR_1pass_wrapper_silver.py
#
#  For the 1st-pass of STAR, map reads using default parameters.
#
#  Authors: Ian McDowell, Brian Jo
#
##################################################
from sys import argv
import math
import glob
import os

split = int(argv[1])
contributors = int(argv[2])
# A little less than an hour per sample, 24 hour job
# split = 15

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/v6p_samples.txt', 'r') as f:
    samples_to_do = [line.strip() for line in f]

# with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/STAR_1pass_hg38_samples.txt', 'r') as f:
#     samples_done = [line.strip() for line in f]

# Check for completed samples:
star_1pass_dir = '/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/STAR_1pass_hg38/'
directories = glob.glob(star_1pass_dir + '*')
samples_done = [x for x in directories if 'Log.final.out' in os.listdir(x + '/')]
samples_done = [str.split(x, '/')[-1] for x in samples_done]

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/STAR_1pass_hg38_samples.txt', 'w') as f:
    [f.write(x + '\n') for x in samples_done]
    f.close()

samples_to_do = [sample for sample in samples_to_do if sample not in set(samples_done)]
samples_to_do_split = [samples_to_do[x:x+split] for x in range(0, len(samples_to_do), split)]
jobs_per_contrib = math.ceil(len(samples_to_do_split) / len(contributors))

counter = 0
for contributor in range(contributors):
    n = str(contributor + 1)
    master_script = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/STAR_1pass_wrapper_" + n + ".sh"
    master_handle=open(master_script, 'w')
    master_handle.write("#!/bin/bash\n\n")
    for i in range(jobs_per_contrib):
        sbatchfile = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/STAR_1pass_" + n + "_part_" + (str(i+1)) + ".sh"
        sbatchhandle=open(sbatchfile, 'w')
        header=r"""#!/bin/bash
#SBATCH -J %s_%s_1pass      # job name
#SBATCH -o /tigress/BEE/gtex/results/group_general/joblogs/STAR_1pass_hg38/%s_%s_1pass.txt
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

# PASS 1 - alignReads
"""%(n, str(i+1), n, str(i+1))
        sbatchhandle.write(header)
        if counter >= len(samples_to_do_split):
            break
        for sample in samples_to_do_split[counter]:
            sbatch="""

SAMPLE=%s

STAR_genomeDir=/tigress/BEE/gtex/data/silver/STAR_hg38_1pass/

runDir=/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/STAR_1pass_hg38/${SAMPLE}
rm -rf $runDir
mkdir -p -m 775 $runDir
cd $runDir

bzip2 -dc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/${SAMPLE}_1.trimmed.P.fastq.bz2 > /scratch/gpfs/${USER}/${SAMPLE}_1.fastq
bzip2 -dc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/${SAMPLE}_2.trimmed.P.fastq.bz2 > /scratch/gpfs/${USER}/${SAMPLE}_2.fastq

/tigress/BEE/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--genomeLoad LoadAndKeep \
--genomeDir $STAR_genomeDir \
--readFilesIn /scratch/gpfs/${USER}/${SAMPLE}_1.fastq /scratch/gpfs/${USER}/${SAMPLE}_2.fastq \
--runThreadN 1 \

# remove temporary unzipped read files
rm -f /scratch/gpfs/${USER}/${SAMPLE}_1.fastq
rm -f /scratch/gpfs/${USER}/${SAMPLE}_2.fastq

# remove first alignment file (alignments will be saved after 2-pass mapping completed)
# only splice junction database is of interest for the first round of mappings
rm -f Aligned.out.sam

            """%(sample)
            sbatchhandle.write(sbatch)
        sbatchhandle.close()
        counter = counter + 1
        master_handle.write("sbatch " + sbatchfile  + " \n")
    master_handle.close()
    print('sh ' + master_script)
