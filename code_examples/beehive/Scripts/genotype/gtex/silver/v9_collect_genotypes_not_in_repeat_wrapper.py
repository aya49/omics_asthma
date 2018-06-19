##################################################
#  v9_collect_genotypes_not_in_repeat_wrapper.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/v9_collect_genotypes_not_in_repeat_wrapper.py
# 
#  This script saves a separate file of genotypes that are not in repeat elements
#
#  Author: Brian Jo
#
##################################################

master_script = "/tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/v9_collect_genotypes_not_in_repeat.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

for i in [str(x+1) for x in range(22)] + ['X']:
  sbatchfile = '/tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/batch/v9_collect_genotypes_not_in_repeat_' + i + '.slurm'
  sbatchhandle=open(sbatchfile, 'w')
  cmd=r"""#!/bin/bash
#SBATCH -J geno_filter_%s      # job name
#SBATCH --mem=32000             # 32 GB requested
#SBATCH -t 24:00:00
#SBATCH -e /tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/logs/v9_collect_genotypes_not_in_repeat_%s.err           # err output directory
#SBATCH -o /tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/logs/v9_collect_genotypes_not_in_repeat_%s.out           # out output directory        

Rscript /tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/v9_collect_genotypes_not_in_repeat.R %s
"""%(i, i, i, i)
  sbatchhandle.write(cmd)
  sbatchhandle.close()
  master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))