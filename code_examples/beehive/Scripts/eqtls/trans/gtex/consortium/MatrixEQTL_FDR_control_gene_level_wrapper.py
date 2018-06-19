##################################################
#  MatrixEQTL_FDR_control_wrapper.py
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/MatrixEQTL_FDR_control_wrapper.py
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

import glob
import os
import os.path

proj_dir = os.environ['proj']

in_path = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/gene_level_FDR/'
joblog_dir = proj_dir + '/Output/joblogs/trans-mapping/gtex/FDR_control/'
out_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/gene_level_FDR/tissues/'

master_script = proj_dir + "/Scripts/eqtls/trans/gtex/batch/MatrixEQTL_gene_FDR_control_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

os.chdir(in_path)
tissues = glob.glob('*')
for tissue in tissues:
	if tissue != 'genes_represented':
		sbatchfile = proj_dir + "/Scripts/eqtls/trans/gtex/batch/MatrixEQTL_gene_FDR_control_" + tissue + "_wrapper.slurm"
		sbatchhandle=open(sbatchfile, 'w')
		cmd=r"""#!/bin/bash
#SBATCH -J FDR_%s      # job name
#SBATCH --mem=12000     # 12 GB requested
#SBATCH -t 24:00:00     
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

/usr/bin/Rscript %s/Scripts/eqtls/trans/gtex/consortium/MatrixEQTL_FDR_control_gene_level.R \
%s %s %s
"""%(tissue, joblog_dir+'gene_FDR_'+tissue+'.err', joblog_dir+'gene_FDR_'+tissue+'.out', proj_dir, in_path, tissue, out_dir)
		sbatchhandle.write(cmd)
		sbatchhandle.close()
		master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))