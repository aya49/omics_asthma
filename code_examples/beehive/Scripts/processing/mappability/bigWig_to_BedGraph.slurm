#!/bin/bash
#SBATCH -J gem_mappability      # job name
#SBATCH --mem=40000            # 40 GB requested
#SBATCH -t 24:00:00            # to be placed in the short queue
#SBATCH -D /tigress/BEE/RNAseq/Output/processing/mappability   # set working directory

export PATH=/tigress/BEE/bin/:$PATH

# wigToBigWig hg38_75mer.wig_mod hg38_75mer.sizes_mod hg38_75mer.bigWig
# bigWigToBedGraph hg38_75mer.bigWig hg38_75mer.bed
wigToBigWig hg38_36mer.wig_mod hg38_36mer.sizes_mod hg38_36mer.bigWig
bigWigToBedGraph hg38_36mer.bigWig hg38_36mer.bed
