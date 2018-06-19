#!/bin/bash
#SBATCH -J gem_mappability      # job name
#SBATCH --mem=40000            # 40 GB requested
#SBATCH -t 06:00:00            # to be placed in the short queue
#SBATCH -n 16             # to be placed in the short queue
#SBATCH -D /tigress/BEE/RNAseq/Output/processing/mappability   # set working directory

export PATH=/tigress/BEE/bin/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin:$PATH

# gem-indexer -T 16 -c dna -i /tigress/BEE/gtex/dbGaP-7716/57610/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta -o hg38_gem_index
# gem-mappability -T 16 -I hg38_gem_index.gem -l 75 -o hg38_75mer
gem-mappability -T 16 -I hg38_gem_index.gem -l 36 -o hg38_36mer
