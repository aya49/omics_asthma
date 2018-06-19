##################################################
#  MatrixEQTL_FDR_control_wrapper.py
#
#  $proj/Scripts/eqtls/trans/gtex/MatrixEQTL_FDR_control_wrapper.py
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
out_path = proj_dir + argv[2]
method = argv[3]
script = argv[4]
tissue_list_file = proj_dir + argv[5]

joblog_dir = proj_dir + '/Output/joblogs/trans-mapping/gtex/FDR_control/'

# Examples
# argv[1] = '/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/'
# argv[2] = '/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/'
# argv[3] = 'all-by-all' - this parameter just for the separation of different runs
# argv[4] = 'MatrixEQTL_FDR_control_PEER_increments.R'
# argv[5] = '/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/tissue_list.txt'

master_script = proj_dir + "/Scripts/eqtls/trans/gtex/batch/MatrixEQTL_FDR_control_" + method + "_wrapper.sh"
master_handle=open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

# os.chdir(in_path)

if not os.path.exists(out_path):
	os.makedirs(out_path)

tissue_list = open(tissue_list_file)

for line in tissue_list.readlines():
	tissue = line.strip()
	print(tissue)
	tissue_dir = in_path + tissue + '/'
	PEER_increments = [y.split('_part')[0] for y in [x.split('PEER')[-1] for x in glob.glob(tissue_dir + '*' + '_part440.RData')]]
	sbatchfile = proj_dir + "/Scripts/eqtls/trans/gtex/batch/MatrixEQTL_FDR_control_" + method + "_" + tissue + "_wrapper.slurm"
	sbatchhandle=open(sbatchfile, 'w')
	cmd=r"""#!/bin/bash
#SBATCH -J %s      # job name
#SBATCH --mem=16000     # 16 GB requested
#SBATCH -t 24:00:00     
#SBATCH -e %s           # err output directory
#SBATCH -o %s           # out output directory
	"""%(tissue, joblog_dir+'FDR_'+tissue+'.err', joblog_dir+'FDR_'+tissue+'.out')
	sbatchhandle.write(cmd)
	for PEER in PEER_increments:	
		cmd=r"""

/usr/bin/Rscript %s/Scripts/eqtls/trans/gtex/%s \
%s %s %s %s
		"""%(proj_dir, script, in_path, tissue, out_path, PEER)
		sbatchhandle.write(cmd)
	sbatchhandle.close()
	master_handle.write("sbatch " + sbatchfile  + " \n")

master_handle.close()

print('sh %s'%(master_script))
