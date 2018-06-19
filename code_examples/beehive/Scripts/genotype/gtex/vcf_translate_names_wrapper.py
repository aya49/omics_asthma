##################################################
#  vcf_translate_names_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_translate_names_wrapper.py
# 
#  This script will translate GTEx custom names to conventional dbSNP (v.142) names
#
#  Author: Ian McDowell
#
##################################################

#!/usr/bin/env python
from sys import argv

maf = argv[1]

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_translate_names_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

chroms = range(1,23)

for i, chrom in enumerate(chroms):
    sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_translate_names_" + str(chrom) + ".slurm"
    sbatchhandle=open(sbatchfile, 'w')
    sbatch=r"""#!/bin/bash
#SBATCH -J chr_%s_vcf_translate
#SBATCH -t 02:00:00
#SBATCH --mem=10000
#SBATCH --get-user-env

cd /tigress/BEE/eQTLs/Data/Genotype/GTEx

CHROM=%s
MAF=%s

python /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_translate_names.py \
imputed_genotypes/vcf/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_ConstrVarIDs.chr${CHROM}.vcf.gz \
auxiliary/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.chr${CHROM}.dict \
imputed_genotypes/vcf/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_${MAF}_HWEp1E6_dbSNPIDs.chr${CHROM}.vcf.gz 

"""% (chrom,chrom,maf)
    sbatchhandle.write(sbatch)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + "\n")

master_handle.close()

print 'Number of jobs', i+1
print 'sh %s'%(master_script)