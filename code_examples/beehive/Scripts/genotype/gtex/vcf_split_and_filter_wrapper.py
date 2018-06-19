##################################################
#  vcf_split_and_filter_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_split_and_filter_wrapper.py
# 
#  This script remove indels, filters SNPs that did not PASS filters, and splits the vcf by chromosome,
#  which will be convenient for later processing.
#
#  Author: Ian McDowell
#
##################################################

#!/usr/bin/env python

from sys import argv

# Input is either maf01 or maf05
maf = argv[1]

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_split_and_filter_" + maf + "_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

chroms = range(1,23)

for i, chrom in enumerate(chroms):
    sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_split_and_filter_" + maf + "_" + str(chrom) + ".slurm"
    sbatchhandle=open(sbatchfile, 'w')
    sbatch=r"""#!/bin/bash
#SBATCH -J chr_%s_vcf_split
#SBATCH -t 03:00:00
#SBATCH --mem=10000
#SBATCH --get-user-env

cd /tigress/BEE/eQTLs/Data/Genotype/GTEx

CHROM=%s

vcftools \
--recode \
--remove-filtered HWEPVALLT1E6 \
--remove-filtered %s \
--remove-filtered info4 \
--recode-INFO-all  \
--remove-indels \
--gzvcf imputed_genotypes/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_%s_HWEp1E6_ConstrVarIDs.vcf.gz \
--chr ${CHROM} \
--stdout | gzip -c > imputed_genotypes/vcf/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_%s_HWEp1E6_ConstrVarIDs.chr${CHROM}.vcf.gz

"""% (chrom,chrom,maf,maf,maf)
    sbatchhandle.write(sbatch)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + "\n")

master_handle.close()

print 'Number of jobs', i+1
print 'sh %s'%(master_script)