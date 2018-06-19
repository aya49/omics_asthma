##################################################
#  write_amish_ciseqtls.py
#
#  $proj/Scripts/eqtls/amish_pipeline/write_amish_ciseqtls.py
#
#  Script for generating comparison stats for Amish and GTEx cis-eQTLs
#
#  Author: Brian Jo
#
##################################################

import glob
import pandas as pd
import gzip
import os.path

gemma_results = '/tigress/BEE/amish/analyses/ciseqtl/gemma_output_per_gene/'
output_f = '/tigress/BEE/amish/analyses/ciseqtl/genomewide/cis_eqtls_1mb_chr'

# Gene annotation for hg19
gene_metadata_f = '/tigress/BEE/RNAseq/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt'
gene_metadata = pd.read_csv(gene_metadata_f, sep='\t')
dist_thresh = 1000000

# for each gene file:
# get the gene location, read the file, and save all SNPs that are within 1MB of the gene TSS
files = glob.glob(gemma_results + '*')
for f in files:
	gene_f = gzip.open(f, 'rb')
	gene_id = f.split('/')[-1].split('.')[0] + '.' + f.split('/')[-1].split('.')[1]
	print(gene_id)
	ind = [x for x in range(gene_metadata.shape[0]) if gene_metadata['gene_id'][x] == gene_id]
	if len(ind) == 0:
		gene_f.close()
		continue
	chrom = gene_metadata['chr'][ind[0]]
	start = gene_metadata['start'][ind[0]]
	end = gene_metadata['end'][ind[0]]
	# Header
	header = gene_f.readline()
	# Separate output by chromosomes
	chr_file = output_f + chrom + '.txt'
	# if this is a new file -
	if not os.path.isfile(chr_file):
		out_f = open(chr_file, 'w')
		out_f.write('gene_id' + '\t' + header.decode('utf-8'))
	else:
		out_f = open(chr_file, 'a')
	for line in gene_f.readlines():
		entry = line.decode('utf-8').split('\t')
		# check for cis distance (changed to 1MB from TSS)
		# if ((int(entry[2]) >= (start - dist_thresh)) and (int(entry[2]) <= (end + dist_thresh))):
		if ((int(entry[2]) >= (start - dist_thresh)) and (int(entry[2]) <= (start + dist_thresh))):
			out_f.write(gene_id + '\t' + line.decode('utf-8'))
	out_f.close()