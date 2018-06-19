##################################################
#  Final_trans_list_wrapper.py
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/Final_trans_list_wrapper.py
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

from sys import argv
import glob
import os
import os.path

proj_dir = os.environ['proj']

in_path = proj_dir + argv[1]
dist_thresh = argv[2]
tissue_list_file = proj_dir + argv[3]

joblog_dir = proj_dir + '/Output/joblogs/trans-mapping/gtex/FDR_control/'

# Examples
# argv[1] = '/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/'
# argv[2] = '1e5'
# argv[3] = '/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/tissue_list.txt'

master_script = proj_dir + "/Scripts/eqtls/trans/gtex/batch/Final_trans_list_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

tissue_list = open(tissue_list_file)

for line in tissue_list.readlines():
	tissue = line.strip()
	print(tissue)
	sbatchfile = proj_dir + "/Scripts/eqtls/trans/gtex/batch/Final_trans_list_" + tissue + ".slurm"
	sbatchhandle=open(sbatchfile, 'w')
	cmd=r"""#!/bin/bash
#SBATCH -J %s      # job name
#SBATCH --mem=16000     # 16 GB requested
#SBATCH -t 24:00:00     
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory

/usr/bin/Rscript %s/Scripts/eqtls/trans/gtex/consortium/Final_trans_list_PEER.R \
%s %s %s
	"""%(tissue, joblog_dir+'Final_list_'+tissue+'.err', joblog_dir+'Final_list_'+tissue+'.out', proj_dir, in_path, dist_thresh, tissue)
	sbatchhandle.write(cmd)
	sbatchhandle.close()
	master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))