##################################################
#  kallisto_check_finished_silver.py
#
#  $proj/Scripts/processing/silver/kallisto_check_finished_silver.py
#
#  For checking whether the bootstraps are finished
#
#  Author: Brian Jo
#
##################################################
import glob
import os
import subprocess

# Check for completed samples:
kallisto_dir = '/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/kallisto_hg38/'
directories = glob.glob(kallisto_dir + '*')
# modify later
samples_done = [x for x in directories if 'abundance.tsv' in os.listdir(x + '/')]
samples_done = [str.split(x, '/')[-1] for x in samples_done]

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/kallisto_hg38_samples.txt', 'w') as f:
    [f.write(x + '\n') for x in samples_done]
    f.close()

bootstrap_samples_done = [x for x in directories if subprocess.call(['h5debug', x + '/abundance.h5'], stdout=open(os.devnull, 'wb')) == 0]
bootstrap_samples_done = [str.split(x, '/')[-1] for x in bootstrap_samples_done]

with open('/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/kallisto_hg38_bootstrap_samples.txt', 'w') as f:
    [f.write(x + '\n') for x in bootstrap_samples_done]
    f.close()