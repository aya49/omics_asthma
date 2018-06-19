#!/bin/bash
#SBATCH -J gem_mappability      # job name
#SBATCH --mem=40000            # 40 GB requested
#SBATCH -t 24:00:00            # to be placed in the short queue
#SBATCH -D /tigress/BEE/RNAseq/Output/processing/mappability   # set working directory

export PATH=/tigress/BEE/bin/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin:$PATH

# gem-2-wig -I hg38_gem_index.gem -i hg38_75mer.mappability -o hg38_75mer
gem-2-wig -I hg38_gem_index.gem -i hg38_36mer.mappability -o hg38_36mer
