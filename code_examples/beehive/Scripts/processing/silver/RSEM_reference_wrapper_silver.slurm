#!/bin/bash
#SBATCH -J RSEM_ref
#SBATCH --get-user-env
#SBATCH --mem=35000           # 35 GB requested
#SBATCH -t 24:00:00            # 24 hours walltime

/tigress/BEE/bin/RSEM-1.3.0/rsem-prepare-reference \
--gtf /tigress/BEE/gtex/external_sources/hg38_silver/gencode.v26.annotation.mod.gtf \
/tigress/BEE/gtex/external_sources/hg38_silver/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
/tigress/BEE/gtex/data/silver/RSEM_hg38