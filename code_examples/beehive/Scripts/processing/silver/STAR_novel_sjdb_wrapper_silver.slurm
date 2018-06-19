#!/bin/bash
#SBATCH -J SJDB
#SBATCH --get-user-env
#SBATCH --mem=35000           # 35 GB requested
#SBATCH -t 24:00:00            # 24 hours walltime

python /tigress/BEE/RNAseq/Scripts/processing/silver/STAR_novel_sjdb_silver.py \
/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/STAR_1pass_hg38/ \
/tigress/BEE/gtex/data/silver/STAR_hg38_1pass/STAR_1pass_sjdb.p \
/tigress/BEE/gtex/data/silver/STAR_hg38_1pass/Number_of_splice_sites_by_number_of_samples_hg19.png \
/tigress/BEE/gtex/data/silver/STAR_hg38_1pass/SJ.out.tab.Pass1.one.percent.sjdb