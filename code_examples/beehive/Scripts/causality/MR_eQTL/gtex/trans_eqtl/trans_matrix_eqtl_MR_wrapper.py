##################################################
#  trans_matrix_eqtl_MR_prep.R
#
#  $proj/Scripts/causality/MR_eQTL/gtex/trans_eqtl/trans_matrix_eqtl_MR_wrapper.py
# 
#  Prepare the list of associations for cis and trans genes for each trans-eQTL pair
#
#  Author: Brian Jo
#
##################################################

from sys import argv
import glob
import os
import os.path

# Examples
argv = ['' for i in range(6)]
argv[1] = '/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/Final_trans_eQTL_list_0.5.txt'
argv[2] = '/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
argv[3] = '/Output/trans-mapping/gtex/MR_eQTL/FDR_0.5/'
argv[4] = '500'
argv[5] = 'all'

proj_dir = os.environ['proj']

input_list = proj_dir + argv[1]
input_dir = proj_dir + argv[2]
output_dir = proj_dir + argv[3]
part_size = int(argv[4])
tissues = argv[5:len(argv)]

joblog_dir = proj_dir + '/Output/joblogs/causality/MR_eQTL/gtex/'

# args = c(1:7)
# args[1] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all/Final_trans_eQTL_list_0.5.txt'
# args[2] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/thyroid_nonverlapping_certain_autosomes_normalized.txt'
# args[3] = '1'
# args[4] = '2.5e8'
# args[5] = 'thyroid'
# args[6] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
# args[7] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MR_eQTL/FDR_0.5/thyroid_MR_results.txt'

master_script = proj_dir + "/Scripts/causality/MR_eQTL/gtex/trans_eqtl/batch/trans_matrix_eqtl_MR_wrapper_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

os.chdir(input_dir)

if not os.path.exists(out_path):
    os.makedirs(out_path)

if tissues == ['all']:
    tissues = set([x.split('_')[0] for x in glob.glob('*' + 'autosomes_normalized.txt')])

for tissue in tissues:
    sbatchfile = proj_dir + "/Scripts/causality/MR_eQTL/gtex/trans_eqtl/batch/trans_matrix_eqtl_MR_wrapper_wrapper.slurm"
    sbatchhandle=open(sbatchfile, 'w')
    cmd=r"""#!/bin/bash
#SBATCH -J FDR_%s      # job name
#SBATCH --mem=12000     # 12 GB requested
#SBATCH -t 24:00:00     
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

/usr/bin/Rscript %s/Scripts/eqtls/trans/gtex/%s \
%s %s %s
    """%(tissue, joblog_dir+'FDR_'+tissue+'.err', joblog_dir+'FDR_'+tissue+'.out', proj_dir, script, in_path, tissue, out_path)
    sbatchhandle.write(cmd)
    sbatchhandle.close()
    master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))