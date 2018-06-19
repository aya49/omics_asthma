##################################################
#  cis_eqtls_amish_comparison_wrapper.py
#
#  $proj/Scripts/eqtls/amish_pipeline/cis_eqtls_amish_comparison_wrapper.py
#
#  Script for generating comparison stats for Amish and GTEx cis-eQTLs
#
#  Author: Brian Jo
#
##################################################

# args = c(1:8)
# args[1] = '/tigress/BEE/amish/analyses/ciseqtl/genomewide/cis_eqtls_1mb_chr22.txt'
# args[2] = '22'
# args[3] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
# args[4] = 'cellsebvtransformedlymphocytes_v6p_consortium_autosomes_normalized.txt'
# args[5] = 'v6p'
# args[6] = '/tigress/BEE/RNAseq/Data/Genotype/gtex/amish_comparison/chr22.txt'
# args[7] = 'cellsebvtransformedlymphocytes'
# args[8] = '/tigress/BEE/RNAseq/Output/cis-mapping/amish/'

from sys import argv

tissue = 'cellsebvtransformedlymphocytes'
tissue = argv[1]

master_script = '/tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/cis_eqtls_amish_comparison.sh'
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

amish_eqtl_prefix = '/tigress/BEE/amish/analyses/ciseqtl/genomewide/cis_eqtls_1mb_chr'
expr_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
expr_name = tissue + '_v6p_consortium_autosomes_normalized.txt'
version = 'v6p'
genotype_file_prefix = '/tigress/BEE/RNAseq/Data/Genotype/gtex/amish_comparison/chr'
out_dir = '/tigress/BEE/RNAseq/Output/cis-mapping/amish/'

for i in range(22):
  sbatchfile = '/tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/cis_eqtls_amish_comparison_' + tissue + '_' + str(i+1) + '.slurm'
  sbatchhandle=open(sbatchfile, 'w')
  cmd=r"""#!/bin/bash
#SBATCH -J %s           # job name
#SBATCH --mem=32000             # 32 GB requested
#SBATCH -t 24:00:00
#SBATCH -e /tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/cis_eqtls_amish_comparison_%s_%s.err
#SBATCH -o /tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/batch/cis_eqtls_amish_comparison_%s_%s.out        

Rscript /tigress/BEE/RNAseq/Scripts/eqtls/amish_pipeline/cis_eqtls_amish_comparison.R %s %s %s %s %s %s %s %s
"""%(str(i+1), tissue, str(i+1), tissue, str(i+1), amish_eqtl_prefix + str(i+1) + '.txt', str(i+1), expr_dir, expr_name, version, genotype_file_prefix + str(i+1) + '.txt', tissue, out_dir)
  sbatchhandle.write(cmd)
  sbatchhandle.close()
  master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()
print('sh %s'%(master_script))