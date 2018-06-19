##################################################
#  obtain_imputed_genotype_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/obtain_imputed_genotype_wrapper.py
# 
#  This script takes the vcf files and returns the allelic dosage format file in continuous genotypes
#
#  Author: Brian Jo
#
##################################################

vcf_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/vcf/'
output_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous/'

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/obtain_imputed_genotype_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

k = 0
for i in range(22):
    k+=1
    chr_num = str(i+1)
    sbatchfile = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/batch/obtain_imputed_genotype_Chr" + chr_num + ".slurm"
    sbatchhandle=open(sbatchfile, 'w')
    cmd=r"""#!/bin/bash
#SBATCH -J geno_%s_obtain      # job name
#SBATCH --mem=8000             # 8 GB requested
#SBATCH -t 02:00:00            # to be placed in the short queue
#SBATCH -D %s                  # set working directory

/usr/bin/python /tigress/BEE/RNAseq/Scripts/genotype/gtex/obtain_imputed_genotype.py \
%s %s %s

cat genotypesChr%s_temp.txt | grep -ve '^\.' > GTEx_genotypes_maf05_continuous_Chr%s_Final.txt;
rm genotypesChr%s_temp.txt

"""%(i, output_dir, vcf_dir, output_dir, chr_num, chr_num, chr_num, chr_num)
    sbatchhandle.write(cmd)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print 'sh %s'%(master_script)
print 'Number of jobs:', k+1
