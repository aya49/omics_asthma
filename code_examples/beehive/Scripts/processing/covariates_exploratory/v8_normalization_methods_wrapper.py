##################################################
#  v8_normalization_methods_wrapper.py
#
#  $proj/Scripts/processing/covariates_exploratory/v8_normalization_methods_wrapper.py
#
#  Trying out various normalization methods/factor analysis in v8 expression data, without the known covariates
#
#  Authors: Brian Jo
#
##################################################
import pandas as pd
import subprocess as sp
from sys import argv
import os

proj_dir = os.environ['proj']

tables_dir = proj_dir + '/Data/Resources/gtex/tables/v8/'
tissue_table = pd.read_csv(tables_dir + 'tissue_table.txt', sep='\t')

master_script = proj_dir + "/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods.sh"
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

tissue_table.index = tissue_table['tissue_name']
# Only process for tissues with samples over 10
tissue_table = tissue_table.loc[tissue_table['num_samples'] > 10]

n_factors_array = [500]

methods = argv[1:]


# Script for SFAmix
# We will try out gene qn and count matrices
if 'sfamix' in methods:
	out_dir = proj_dir + '/Output/processing/exploratory/v8/sfamix/no_cov/'
	exp_dir = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/'
	# Populate checkpoint tables
	log_checkpoint_table = out_dir + 'log_transform/checkpoint_table.txt'
	qn_checkpoint_table = out_dir + 'quantile_norm/sample_norm/checkpoint_table.txt'
	# write headers
	log_f = open(log_checkpoint_table, 'w')
	log_f.write('\t' + '\t'.join([str(x) for x in n_factors_array]) + '\n')
	qn_f = open(qn_checkpoint_table, 'w')
	qn_f.write('\t' + '\t'.join([str(x) for x in n_factors_array]) + '\n')
	for tissue in tissue_table.index:
		log_f.write(tissue + '\t')
		qn_f.write(tissue + '\t')
		file_checks = [os.path.exists(out_dir + 'log_transform/' + str(x) + '/' + tissue + '/final') for x in n_factors_array]
		log_f.write('\t'.join(str(x) for x in file_checks) + '\n')
		file_checks = [os.path.exists(out_dir + 'quantile_norm/sample_norm/' + str(x) + '/' + tissue + '/final') for x in n_factors_array]
		qn_f.write('\t'.join(str(x) for x in file_checks) + '\n')
	log_f.close()
	qn_f.close()
	log_table = pd.read_csv(log_checkpoint_table, sep='\t', index_col=0)
	qn_table = pd.read_csv(qn_checkpoint_table, sep='\t', index_col=0)		
	# now make a list of jobs to run
	for tissue in tissue_table.index:
		# Is the inverted matrix in scratch?
		exp_mat = 'v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt'
		if not os.path.exists('/scratch/gpfs/bj5/' + exp_mat):
			# Script for inverting the matrix
			sbatchhandle.write('Rscript ' + proj_dir + '/Scripts/processing/covariates_exploratory/v8_sfamix_invert_matrix.R ' + exp_dir + '/quantile_norm/ ' + exp_mat + '\n')
		# Is the inverted matrix in scratch?
		exp_mat = 'v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt'
		if not os.path.exists('/scratch/gpfs/bj5/' + exp_mat):
			# Script for inverting the matrix
			sbatchhandle.write('Rscript ' + proj_dir + '/Scripts/processing/covariates_exploratory/v8_sfamix_invert_matrix.R ' + exp_dir + '/log_transform/ ' + exp_mat + '\n')
		for n_factors in n_factors_array:
			# process qn files - Create directory
			directory = out_dir + 'quantile_norm/sample_norm/' + str(n_factors)
			if not os.path.exists(directory):
				os.makedirs(directory)
			directory = directory + '/' + tissue
			if not os.path.exists(directory):
				os.makedirs(directory)
			# If the output is not there:
			if not qn_table.loc[tissue][str(n_factors)]:
				sbatchfile = proj_dir + '/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods_sfamix_qn_' + tissue + '_' + str(n_factors) + '.sh'
				sbatchhandle=open(sbatchfile, 'w')
				header=r"""#!/bin/bash
#SBATCH -J %s_sfamix      # job name
#SBATCH -o %s/Scripts/processing/covariates_exploratory/joblogs/qn_%s_%s_sfamix
#SBATCH --get-user-env
#SBATCH --mem=50000           # 50 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 6-00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

"""%(tissue, proj_dir, tissue, str(n_factors))
				sbatchhandle.write(header)
				exp_mat = '/scratch/gpfs/bj5/v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt'
				param_array = ['/tigress/BEE/bin/SFAmix_code_documentation/SFAmix', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
				sbatchhandle.close()
				master_handle.write("sbatch " + sbatchfile + '\n')
			# process log_transform files - Create directory
			directory = out_dir + 'log_transform/' + str(n_factors)
			if not os.path.exists(directory):
				os.makedirs(directory)
			directory = directory + '/' + tissue
			if not os.path.exists(directory):
				os.makedirs(directory)
			# If the output is not there:
			if not log_table.loc[tissue][str(n_factors)]:
				sbatchfile = proj_dir + '/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods_sfamix_log_' + tissue + '_' + str(n_factors) + '.sh'
				sbatchhandle=open(sbatchfile, 'w')
				header=r"""#!/bin/bash
#SBATCH -J %s_sfamix      # job name
#SBATCH -o %s/Scripts/processing/covariates_exploratory/joblogs/log_%s_%s_sfamix
#SBATCH --get-user-env
#SBATCH --mem=50000           # 50 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 6-00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

"""%(tissue, proj_dir, tissue, str(n_factors))
				sbatchhandle.write(header)
				exp_mat = '/scratch/gpfs/bj5/v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt'
				param_array = ['/tigress/BEE/bin/SFAmix_code_documentation/SFAmix', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
				sbatchhandle.close()
				master_handle.write("sbatch " + sbatchfile + '\n')


# Script for PEER factors
if 'peer' in methods:
	out_dir = proj_dir + '/Output/processing/exploratory/v8/peer/'
	exp_dir = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/'
	# Populate checkpoint tables
	log_checkpoint_table = out_dir + 'log_transform/checkpoint_table.txt'
	qn_checkpoint_table = out_dir + 'quantile_norm/sample_norm/checkpoint_table.txt'
	# write headers
	log_f = open(log_checkpoint_table, 'w')
	log_f.write('\t' + '\t'.join([str(x) for x in n_factors_array]) + '\n')
	qn_f = open(qn_checkpoint_table, 'w')
	qn_f.write('\t' + '\t'.join([str(x) for x in n_factors_array]) + '\n')
	for tissue in tissue_table.index:
		log_f.write(tissue + '\t')
		qn_f.write(tissue + '\t')
		file_checks = [os.path.exists(out_dir + 'log_transform/' + str(x) + '/' + tissue + '/final') for x in n_factors_array]
		log_f.write('\t'.join(str(x) for x in file_checks) + '\n')
		file_checks = [os.path.exists(out_dir + 'quantile_norm/sample_norm/' + str(x) + '/' + tissue + '/final') for x in n_factors_array]
		qn_f.write('\t'.join(str(x) for x in file_checks) + '\n')
	log_f.close()
	qn_f.close()
	log_table = pd.read_csv(log_checkpoint_table, sep='\t', index_col=0)
	qn_table = pd.read_csv(qn_checkpoint_table, sep='\t', index_col=0)		
	# now make a list of jobs to run
	for tissue in tissue_table.index:
		# Is the inverted matrix in scratch?
		qn_exp_mat = 'v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt'
		log_exp_mat = 'v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt'
		for n_factors in n_factors_array:
			# process qn files - Create directory
			directory = out_dir + 'quantile_norm/sample_norm/' + str(n_factors)
			if not os.path.exists(directory):
				os.makedirs(directory)
			directory = directory + '/' + tissue
			if not os.path.exists(directory):
				os.makedirs(directory)
			# If the output is not there:
			if not qn_table.loc[tissue][str(n_factors)]:
				sbatchfile = proj_dir + '/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods_peer_qn_' + tissue + '_' + str(n_factors) + '.sh'
				sbatchhandle=open(sbatchfile, 'w')
				header=r"""#!/bin/bash
#SBATCH -J %s_peer      # job name
#SBATCH -o %s/Scripts/processing/covariates_exploratory/joblogs/qn_%s_%s_peer
#SBATCH --get-user-env
#SBATCH --mem=50000           # 50 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

"""%(tissue, proj_dir, tissue, str(n_factors))
				sbatchhandle.write(header)
				exp_mat = exp_dir + 'quantile_norm/' + qn_exp_mat
				# param_array = ['', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
				sbatchhandle.close()
				master_handle.write("sbatch " + sbatchfile + '\n')
			# process log_transform files - Create directory
			directory = out_dir + 'log_transform/' + str(n_factors)
			if not os.path.exists(directory):
				os.makedirs(directory)
			directory = directory + '/' + tissue
			if not os.path.exists(directory):
				os.makedirs(directory)
			# If the output is not there:
			if not log_table.loc[tissue][str(n_factors)]:
				sbatchfile = proj_dir + '/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods_peer_log_' + tissue + '_' + str(n_factors) + '.sh'
				sbatchhandle=open(sbatchfile, 'w')
				header=r"""#!/bin/bash
#SBATCH -J %s_peer      # job name
#SBATCH -o %s/Scripts/processing/covariates_exploratory/joblogs/log_%s_%s_peer
#SBATCH --get-user-env
#SBATCH --mem=50000           # 50 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

"""%(tissue, proj_dir, tissue, str(n_factors))
				sbatchhandle.write(header)
				exp_mat = exp_dir + 'log_transform/' + log_exp_mat
				# param_array = ['/tigress/BEE/bin/SFAmix_code_documentation/SFAmix', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
				sbatchhandle.close()
				master_handle.write("sbatch " + sbatchfile + '\n')


print("sh " + master_script)
master_handle.close()
