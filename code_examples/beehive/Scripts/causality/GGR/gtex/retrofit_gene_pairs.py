##################################################
#  retrofit_gene_pairs.py
#
#  $proj/Scripts/causality/GGR/gtex/retrofit_gene_pairs.py
# 
#  Changes gene names from hg38 (gencode v.22) to hg19 (gencode v.19), and also adds information about gene positions.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

from sys import argv
import glob
import os
import os.path
import pandas as pd
import json
import numpy as np

proj_dir = os.environ['proj']

GGR_home_dir = proj_dir + argv[1]
condition = argv[2]
output_dir = proj_dir + argv[3]
gene_metadata_v19 = proj_dir + argv[4]

# Examples
# GGR_home_dir = '/tigress/BEE/RNAseq/Data/Networks/GGR/'
# condition = '0mean-1var_g-null_l-fdr'
# output_dir = '/tigress/BEE/RNAseq/Data/Networks/GGR/retrofitted/'
# gene_metadata_v19 = '/tigress/BEE/RNAseq/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt'

files_to_modify = glob.glob(GGR_home_dir + condition + '/networks/' + '*')

# Given the gene id, return the chromosome and start position of the gene
# def return_chr_pos(gene_id, gene_name, gene_metadata):
# 	geneID = 'NA'
# 	geneChr = 'NA'
# 	startPos = 0
# 	if gene_id != None and gene_id in np.array(gene_metadata['ensg']):
# 		geneID = str(np.array(gene_metadata[gene_metadata['ensg'] == gene_id]['ensg'])[0])
# 		geneChr = str(np.array(gene_metadata[gene_metadata['ensg'] == gene_id]['chr'])[0])
# 		startPos = int(gene_metadata[gene_metadata['ensg'] == gene_id]['start'])
# 	elif gene_name in np.array(gene_metadata['gene_name']):
# 		geneID = str(np.array(gene_metadata[gene_metadata['gene_name'] == gene_name]['ensg'])[0])
# 		geneChr = str(np.array(gene_metadata[gene_metadata['gene_name'] == gene_name]['chr'])[0])
# 		startPos = int(gene_metadata[gene_metadata['gene_name'] == gene_name]['start'])
# 	return((geneID, geneChr, startPos))

def return_chr_pos(gene_id, v22_map, gene_metadata):
	gene_id_trunc = gene_id.split('.')[0]
	geneID = 'NA'
	geneChr = 'NA'
	startPos = 0	
	# If there the ensg id is consistent:
	if gene_id_trunc in np.array(gene_metadata['ensg']):
		geneID = str(np.array(gene_metadata[gene_metadata['ensg'] == gene_id_trunc]['ensg'])[0])
		geneChr = str(np.array(gene_metadata[gene_metadata['ensg'] == gene_id_trunc]['chr'])[0])
		startPos = int(np.array(gene_metadata[gene_metadata['ensg'] == gene_id_trunc]['start'])[0])
	# If not, search by gene name:
	elif v22_map[gene_id] in np.array(gene_metadata['gene_name']):
		gene_name = v22_map[gene_id]
		geneID = str(np.array(gene_metadata[gene_metadata['gene_name'] == gene_name]['ensg'])[0])
		geneChr = str(np.array(gene_metadata[gene_metadata['gene_name'] == gene_name]['chr'])[0])
		startPos = int(np.array(gene_metadata[gene_metadata['gene_name'] == gene_name]['start'])[0])
	return((geneID, geneChr, startPos))

f = open('/tigress/BEE/RNAseq/Data/Resources/annotations/gencode.v22.gene_id_to_gene_name.json.txt')
v22_map = json.load(f)

gene_metadata = pd.read_csv(gene_metadata_v19, sep = '\t')
gene_metadata['ensg'] = [x.split('.')[0] for x in np.array(gene_metadata['gene_id'])]

for network_f in files_to_modify:
	filename = str.split(network_f, '/')[-1]
	print(filename)

	f = open(network_f)
	if not os.path.exists(output_dir + condition):
	    os.makedirs(output_dir + condition)

	f_out = open(output_dir + condition + '/' + filename, 'w')
	header = f.readline()
	mod_header = header.strip() + '\tCause_v19\tgeneChr_cis\tstartPos_cis\tEffect_v19\tgeneChr_trans\tstartPos_trans\n'
	f_out.write(mod_header)

	for line in f.readlines():
		entry = line.strip().split('\t')
		cause = entry[1]
		effect = entry[2]
		(geneID_cis, geneChr_cis, startPos_cis) = return_chr_pos(cause, v22_map, gene_metadata)
		(geneID_trans, geneChr_trans, startPos_trans) = return_chr_pos(effect, v22_map, gene_metadata)
		# if (geneChr_cis == 'NA' or geneChr_trans == 'NA'):
		# 	print('Unable to map : ')
		# 	print(cause + '\t' + geneID_cis)
		# 	print(effect + '\t' + geneID_trans)
		# 	continue
		f_out.write(line.strip() + '\t' + '\t'.join([geneID_cis, geneChr_cis, str(startPos_cis), geneID_trans, geneChr_trans, str(startPos_trans)]) + '\n')

	f_out.close()
	f.close()
