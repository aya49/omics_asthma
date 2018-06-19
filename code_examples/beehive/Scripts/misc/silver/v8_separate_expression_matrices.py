##################################################
#  v8_separate_expression_matrices.py
#
#  $proj/Scripts/misc/silver/v8_separate_expression_matrices.py
# 
#  This script separates the GTEx expression matrices by tissues
#
#  Author: Brian Jo
#
##################################################

import os
from sys import argv
import gzip
proj_dir = os.environ['proj']
tissue = argv[1]
# tissue = 'arterytibial'

# sample table
sample_table = proj_dir + '/Data/Resources/gtex/tables/v8/sample_table.txt'
f = open(sample_table)
f.readline()

tissue_samps = []
tissue_subjs = []
ind = 0
for line in f.readlines():
	samp_tissue = str.split(line.strip(), '\t')[-1]
	if samp_tissue == tissue:
		tissue_samps.append(str.split(line, '\t')[0])
		tissue_subjs.append(str.split(line, '\t')[2])
	ind = ind + 1

# gene expression files - tpm
RSEM_gene_tpm = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/dbGaP/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz'

f = gzip.open(RSEM_gene_tpm, 'rb')
f.readline()
f.readline()

header = f.readline()
samples = header.strip().decode("utf-8").split('\t')
samples = samples[2:]

tissue_inds = [x for x in range(len(samples)) if samples[x] in tissue_samps]
tissue_subjects = [samples[x].split('-')[0] + '-' + samples[x].split('-')[1] for x in tissue_inds]

out_file = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/raw/v8_RSEMv1.3.0_gene_tpm_' + tissue + '.txt'
f_out = open(out_file, 'w')
# write header
f_out.write('gene\t' + '\t'.join(tissue_subjects) + '\n')

for line in f.readlines():
	entry = line.strip().decode("utf-8").split('\t')
	gene = entry[0]
	transcripts = entry[1]
	exp_vals = entry[2:]
	tissue_vals = [exp_vals[x] for x in tissue_inds]
	f_out.write(gene + '\t' + '\t'.join(tissue_vals) + '\n')

f_out.close()


# gene expression files - counts
RSEM_gene_count = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/dbGaP/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count.gct.gz'

f = gzip.open(RSEM_gene_count, 'rb')
f.readline()
f.readline()

header = f.readline()
samples = header.strip().decode("utf-8").split('\t')
samples = samples[2:]

tissue_inds = [x for x in range(len(samples)) if samples[x] in tissue_samps]
tissue_subjects = [samples[x].split('-')[0] + '-' + samples[x].split('-')[1] for x in tissue_inds]

out_file = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/raw/v8_RSEMv1.3.0_gene_count_' + tissue + '.txt'
f_out = open(out_file, 'w')
# write header
f_out.write('gene\t' + '\t'.join(tissue_subjects) + '\n')

for line in f.readlines():
	entry = line.strip().decode("utf-8").split('\t')
	gene = entry[0]
	transcripts = entry[1]
	exp_vals = entry[2:]
	tissue_vals = [exp_vals[x] for x in tissue_inds]
	f_out.write(gene + '\t' + '\t'.join(tissue_vals) + '\n')

f_out.close()


# transcript expression files - tpm
RSEM_transcript_tpm = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/dbGaP/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz'

f = gzip.open(RSEM_transcript_tpm, 'rb')
f.readline()
f.readline()

header = f.readline()
samples = header.strip().decode("utf-8").split('\t')
samples = samples[2:]

tissue_inds = [x for x in range(len(samples)) if samples[x] in tissue_samps]
tissue_subjects = [samples[x].split('-')[0] + '-' + samples[x].split('-')[1] for x in tissue_inds]

out_file = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/raw/v8_RSEMv1.3.0_transcript_tpm_' + tissue + '.txt'
f_out = open(out_file, 'w')
# write header
f_out.write('transcript\tgene\t' + '\t'.join(tissue_subjects) + '\n')

for line in f.readlines():
	entry = line.strip().decode("utf-8").split('\t')
	transcript = entry[0]
	gene = entry[1]
	exp_vals = entry[2:]
	tissue_vals = [exp_vals[x] for x in tissue_inds]
	f_out.write(transcript + '\t' + gene + '\t' + '\t'.join(tissue_vals) + '\n')

f_out.close()


# transcript expression files - counts
RSEM_transcript_count = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/dbGaP/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz'

f = gzip.open(RSEM_transcript_count, 'rb')
f.readline()
f.readline()

header = f.readline()
samples = header.strip().decode("utf-8").split('\t')
samples = samples[2:]

tissue_inds = [x for x in range(len(samples)) if samples[x] in tissue_samps]
tissue_subjects = [samples[x].split('-')[0] + '-' + samples[x].split('-')[1] for x in tissue_inds]

out_file = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/raw/v8_RSEMv1.3.0_transcript_count_' + tissue + '.txt'
f_out = open(out_file, 'w')
# write header
f_out.write('transcript\tgene\t' + '\t'.join(tissue_subjects) + '\n')

for line in f.readlines():
	entry = line.strip().decode("utf-8").split('\t')
	transcript = entry[0]
	gene = entry[1]
	exp_vals = entry[2:]
	tissue_vals = [exp_vals[x] for x in tissue_inds]
	f_out.write(transcript + '\t' + gene + '\t' + '\t'.join(tissue_vals) + '\n')

f_out.close()

