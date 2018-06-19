##################################################
#  STAR_2pass_wrapper_silver.py
#
#  $proj/Scripts/processing/silver/STAR_2pass_wrapper_silver.py
#
#  For the 2nd-pass of STAR, map reads using default parameters.
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

# with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/STAR_2pass_hg38_samples.txt', 'r') as f:
#     samples_done = [line.strip() for line in f]

# Check for completed samples:
star_2pass_dir = '/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/STAR_2pass_hg38/'
directories = glob.glob(star_2pass_dir + '*')
samples_done = [x for x in directories if 'Log.final.out' in os.listdir(x + '/')]
samples_done = [str.split(x, '/')[-1] for x in samples_done]

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/STAR_2pass_hg38_samples.txt', 'w') as f:
    [f.write(x + '\n') for x in samples_done]
    f.close()

samples_to_do = [sample for sample in samples_to_do if sample not in set(samples_done)]

f = open('/tigress/BEE/RNAseq/Scripts/processing/silver/STAR_cleanup_silver.sh', 'w')
for sample in samples_to_do:
    f.write("yes|rm -r " + star_2pass_dir + sample + "\n")

f.close()

samples_to_do_split = [samples_to_do[x:x+split] for x in range(0, len(samples_to_do), split)]
jobs_per_contrib = math.ceil(len(samples_to_do_split) / contributors)
print(jobs_per_contrib)

counter = 0
# job names are in numbers - can adjust to the list of contributors
for contributor in range(contributors):
    n = str(contributor + 1)
    master_script = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/STAR_2pass_wrapper_" + n + ".sh"
    master_handle=open(master_script, 'w')
    master_handle.write("#!/bin/bash\n\n")
    for i in range(jobs_per_contrib):
        sbatchfile = "/tigress/BEE/RNAseq/Scripts/processing/silver/batch/STAR_2pass_" + n + "_part_" + (str(i+1)) + ".sh"
        sbatchhandle=open(sbatchfile, 'w')
        header=r"""#!/bin/bash
#SBATCH -J %s_%s_2pass      # job name
#SBATCH -o /tigress/BEE/gtex/results/group_general/joblogs/STAR_2pass_hg38/%s_%s_2pass.txt
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

# PASS 2 - alignReads
"""%(n, str(i+1), n, str(i+1))
        sbatchhandle.write(header)
        if counter >= len(samples_to_do_split):
            break
        for sample in samples_to_do_split[counter]:
            sbatch="""

SAMPLE=%s

STAR_genomeDir=/tigress/BEE/gtex/data/silver/STAR_hg38_2pass/

runDir=/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/STAR_2pass_hg38/${SAMPLE}
rm -rf $runDir
mkdir -p -m 775 $runDir
cd $runDir

bzip2 -dc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/${SAMPLE}_1.trimmed.P.fastq.bz2 > /scratch/gpfs/${USER}/${SAMPLE}_1.fastq
bzip2 -dc /tigress/BEE/gtex/data/phenotype/expression/trimmed_rna_seq_reads/${SAMPLE}_2.trimmed.P.fastq.bz2 > /scratch/gpfs/${USER}/${SAMPLE}_2.fastq

/tigress/BEE/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--genomeLoad LoadAndKeep \
--genomeDir $STAR_genomeDir \
--readFilesIn /scratch/gpfs/${USER}/${SAMPLE}_1.fastq /scratch/gpfs/${USER}/${SAMPLE}_2.fastq \
--outSAMattributes NH HI AS NM MD \
--outFilterMismatchNoverReadLmax 0.04 \
--sjdbScore 1 \
--outFilterType BySJout \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--quantMode TranscriptomeSAM \
--outReadsUnmapped Fastx

rm -f Aligned.out.sam
rm -f /scratch/gpfs/${USER}/${SAMPLE}_1.fastq
rm -f /scratch/gpfs/${USER}/${SAMPLE}_2.fastq

            """%(sample)
            sbatchhandle.write(sbatch)
        sbatchhandle.close()
        counter = counter + 1
        master_handle.write("sbatch " + sbatchfile  + " \n")
    master_handle.close()
    print('sh ' + master_script)

# description of parameters
# --genomeLoad LoadAndKeep \
# --genomeDir $STAR_genomeDir \
# --readFilesIn /scratch/gpfs/${USER}/${SAMPLE}_1.fastq /scratch/gpfs/${USER}/${SAMPLE}_2.fastq \
# --outSAMattributes NH HI AS NM MD \
# # --outFilterIntronMotifs RemoveNoncanonical \
# --outFilterMismatchNoverReadLmax 0.04 \
# --sjdbScore 1 \
# # ENCODE standard options
# --outFilterType BySJout \
# --alignIntronMin 20 \
# --alignIntronMax 1000000 \
# --alignMatesGapMax 1000000 \
# --alignSJoverhangMin 8 \
# --alignSJDBoverhangMin 1 \
# --outFilterMultimapNmax 20 \
# --outFilterMismatchNmax 999 \
# # needed for RSEM \
# --quantMode TranscriptomeSAM \
# # needed for rop \
# --outReadsUnmapped Fastx
