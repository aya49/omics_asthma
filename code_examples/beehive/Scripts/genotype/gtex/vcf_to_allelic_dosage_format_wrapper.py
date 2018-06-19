##################################################
#  vcf_to_allelic_dosage_format_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_to_allelic_dosage_format_wrapper.py
# 
#  This script will convert the genotypes from vcf format to allelic dosage format for eQTL analyses
#
#  Author: Ian McDowell
#
##################################################

#!/usr/bin/env python
from sys import argv

maf = argv[1]

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_to_allelic_dosage_format_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

chroms = range(1,23) + ['X','Y']
# chroms = ['X','Y']

for i, chrom in enumerate(chroms):
    sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_to_allelic_dosage_format_" + str(chrom) + ".slurm"
    sbatchhandle=open(sbatchfile, 'w')
    sbatch=r"""#!/bin/bash
#SBATCH -J %s_vcf_AD
#SBATCH -t 00:63:00
#SBATCH --mem=10000
#SBATCH --get-user-env

cd /tigress/BEE/eQTLs/Data/Genotype/GTEx

CHROM=%s
MAF=%s

plink \
--vcf imputed_genotypes/vcf/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}.vcf.gz \
--recode A \
--out imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}

/tigress/BEE/RNAseq/Scripts/genotype/gtex/rowsToCols \
imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}.raw \
imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}.raw.T

python /tigress/BEE/RNAseq/Scripts/genotype/gtex/plink_recode_A_to_allelic_dosage_format.py \
imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}.raw.T \
imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}.allelic_dosage.txt 

rm imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}.raw.T

"""% (chrom,chrom,maf)
    sbatchhandle.write(sbatch)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + "\n")

master_handle.close()

print 'Number of jobs', i+1
print 'sh %s'%(master_script)