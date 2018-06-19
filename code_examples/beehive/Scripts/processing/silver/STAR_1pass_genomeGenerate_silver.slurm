#!/bin/bash
#SBATCH -J create_STAR_genome
#SBATCH --get-user-env
#SBATCH -N 4                  # 1 is default
#SBATCH --ntasks-per-node=8   # 1 is default
#SBATCH --mem=35000           # 35 GB requested
#SBATCH -t 3:00:00            # 3 hours walltime

STAR_genomeDir=/tigress/BEE/gtex/data/silver/STAR_hg38_1pass
cd $STAR_genomeDir

/tigress/BEE/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runMode genomeGenerate \
--sjdbGTFfile /tigress/BEE/gtex/external_sources/hg38_silver/gencode.v26.annotation.mod.gtf \
--sjdbOverhang 75 \
--genomeDir $STAR_genomeDir \
--genomeFastaFiles /tigress/BEE/gtex/external_sources/hg38_silver/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--runThreadN 32
