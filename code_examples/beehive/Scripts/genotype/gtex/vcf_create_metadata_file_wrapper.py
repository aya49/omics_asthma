##################################################
#  vcf_create_metadata_file_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_create_metadata_file_wrapper.py
# 
#  This script creates the metadata files that store additional information about each genotype which may come in handy
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/env python
from sys import argv

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_create_metadata_file.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

chroms = range(1,23)

for i, chrom in enumerate(chroms):
    sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/vcf_create_metadata_" + str(chrom) + ".slurm"
    sbatchhandle=open(sbatchfile, 'w')
    sbatch=r"""#!/bin/bash
#SBATCH -J chr_%s_metadata
#SBATCH -t 03:00:00
#SBATCH --mem=12000
#SBATCH --get-user-env

cd /tigress/BEE/eQTLs/Data/Genotype/GTEx/

CHROM=%s

/usr/bin/python /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_create_metadata_file.py \
imputed_genotypes/vcf/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.chr${CHROM}.vcf.gz \
auxiliary/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.chr${CHROM}.dict \
SNP_metadata/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr${CHROM}.metadata.txt

"""% (chrom,chrom)
    sbatchhandle.write(sbatch)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + "\n")

master_handle.close()

print 'Number of jobs', i+1
print 'sh %s'%(master_script)