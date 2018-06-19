
master_script = '/tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/write_genotypes_amish_comparison.sh'
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

for i in range(22):
  sbatchfile = '/tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/write_genotypes_amish_comparison_' + str(i+1) + '.slurm'
  sbatchhandle=open(sbatchfile, 'w')
  cmd=r"""#!/bin/bash
#SBATCH -J %s           # job name
#SBATCH --mem=32000             # 32 GB requested
#SBATCH -t 24:00:00
#SBATCH -e /tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/write_genotypes_amish_comparison_%s.err
#SBATCH -o /tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/write_genotypes_amish_comparison_%s.out        

python /tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/write_genotypes_amish_comparison.py %s
"""%(str(i+1), str(i+1), str(i+1), str(i+1))
  sbatchhandle.write(cmd)
  sbatchhandle.close()
  master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()
print('sh %s'%(master_script))