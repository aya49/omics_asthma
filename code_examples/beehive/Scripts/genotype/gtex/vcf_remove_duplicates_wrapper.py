##################################################
#  vcf_remove_duplicates_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_remove_duplicates_wrapper.py
# 
#  It was discovered that at least one SNP is duplicated in one imputed genotype file, this script ensures that
#  there are no duplicates in all genotype files (duplicates are removed as it is unreasonable to speculate which is "correct").
#
#  Author: Ian McDowell
#
##################################################

#!/usr/bin/env python
from sys import argv

maf = argv[1]

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_remove_duplicates_wrapper.%s.sh"%(maf)
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

chroms = range(1,23)

for i, chrom in enumerate(chroms):
    sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_remove_duplicates.%s.%s.slurm"%(chrom,maf)
    sbatchhandle=open(sbatchfile, 'w')
    sbatch=r"""#!/bin/bash
#SBATCH -J chr_%s_vcf_rm_dup
#SBATCH -t 02:00:00
#SBATCH --mem=20000
#SBATCH --get-user-env

CHROM=%s
MAF=%s

cd /tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/vcf

python /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_remove_duplicate_snps.py \
GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_ConstrVarIDs.chr${CHROM}.vcf.gz \
GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_ConstrVarIDs_rm_dups.chr${CHROM}.vcf.gz \
GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_ConstrVarIDs_dup_report.chr${CHROM}.txt


"""% (chrom,chrom,maf)
    sbatchhandle.write(sbatch)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + "\n")

master_handle.close()

print 'Number of jobs', i+1
print 'sh %s'%(master_script)