#!/bin/bash
#SBATCH -J gem_mappability      # job name
#SBATCH --mem=40000            # 40 GB requested
#SBATCH -t 24:00:00            # to be placed in the short queue
#SBATCH -D /tigress/BEE/RNAseq/Output/processing/mappability   # set working directory

export PATH=/tigress/BEE/bin/:$PATH

wig2bed < hg38_75mer.wig > hg38_75mer.bed
