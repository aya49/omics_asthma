##################################################
#  metadata_add_addl_info_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/metadata_add_addl_info_wrapper.py
# 
#  This script creates the metadata files that also have the information about the cis genes and repeat regions
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/env python
from sys import argv

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/metadata_add_addl_info.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

chroms = range(1,23)

for i, chrom in enumerate(chroms):
    sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/metadata_add_addl_info_" + str(chrom) + ".slurm"
    sbatchhandle=open(sbatchfile, 'w')
    sbatch=r"""#!/bin/bash
#SBATCH -J chr_%s_metadata
#SBATCH -t 10:00:00
#SBATCH --mem=12000
#SBATCH --get-user-env

cd /tigress/BEE/eQTLs/Data/Genotype/GTEx/

CHROM=%s

Rscript /tigress/BEE/RNAseq/Scripts/genotype/gtex/metadata_add_addl_info.R \
${CHROM} \
GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr${CHROM}.metadata.txt \
GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr${CHROM}.metadata.expanded.txt \

"""% (chrom,chrom)
    sbatchhandle.write(sbatch)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + "\n")

master_handle.close()

print 'Number of jobs', i+1
print 'sh %s'%(master_script)