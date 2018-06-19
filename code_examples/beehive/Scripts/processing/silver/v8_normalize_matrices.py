##################################################
#  v8_normalize_matrices.py
#
#  $proj/Scripts/processing/silver/v8_normalize_matrices.py
#
#  Matrix normalization script
#
#  Authors: Brian Jo
#
##################################################
import pickle
import pandas as pd
import numpy as np
from sys import argv
from scipy.stats import norm
from scipy.stats import rankdata
from collections import OrderedDict

# exp_dir = '/Users/brian_jo/Desktop/Project/RNAseq_pipeline/Data/Expression/gtex/hg38/v8/raw/'
exp_dir = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/'
# file_prefix = 'v8_RSEMv1.3.0_transcript_tpm_'
# tissue = 'muscleskeletal'
tissue = argv[1]

# Read in the matrices
gene_count = pd.read_csv(exp_dir + 'raw/v8_RSEMv1.3.0_gene_count_' + tissue + '.txt', sep='\t')
transcript_count = pd.read_csv(exp_dir + 'raw/v8_RSEMv1.3.0_transcript_count_' + tissue + '.txt', sep='\t')
gene_tpm = pd.read_csv(exp_dir + 'raw/v8_RSEMv1.3.0_gene_tpm_' + tissue + '.txt', sep='\t')
transcript_tpm = pd.read_csv(exp_dir + 'raw/v8_RSEMv1.3.0_transcript_tpm_' + tissue + '.txt', sep='\t')

# load transcript annotation
transcript_anno_file = '/tigress/BEE/RNAseq/Data/Resources/annotations/silver/gencode.v26.transcripts'
transcript_dict = pickle.load(open(transcript_anno_file, 'rb'))

print('Read in input files')

# filter gene and transcript tables based on minimal gene expression
# our threshold is at least 0.1 TPM for at least 10 individuals
gene_tpm_filter = pd.DataFrame()
transcript_tpm_filter = pd.DataFrame()
gene_count_filter = pd.DataFrame()
transcript_count_filter = pd.DataFrame()

temp_transcript_tpm = pd.DataFrame()
temp_transcript_count = pd.DataFrame()
all_chrs = [str(x) for x in range(1,23)] + ['X','Y','M']

filtered_gene_list = []
filtered_transcript_list = []

# With ordered dict, transcripts come in gencode order
cur_gene = ''
for cur_chr in all_chrs:
	print(cur_chr)
	for k,v in transcript_dict[cur_chr].items():
		# is the transcript in the isoform table?
		if k not in transcript_tpm.index:
			continue
		# new gene
		if v['gene_id'] != cur_gene:
			if cur_gene != '':
				# Condition: At least 10 individuals with TPM at least 0.1
				if (sum(gene_tpm.loc[cur_gene] >= 0.1) >= 10):
					filtered_gene_list.append(cur_gene)
					filtered_transcript_list = filtered_transcript_list + temp_transcript_tpm.index.tolist()
					gene_tpm_filter = gene_tpm_filter.append(gene_tpm.loc[[cur_gene]])
					gene_count_filter = gene_count_filter.append(gene_count.loc[[cur_gene]])
					transcript_tpm_filter = transcript_tpm_filter.append(temp_transcript_tpm)
					transcript_count_filter = transcript_count_filter.append(temp_transcript_count)
			cur_gene = v['gene_id']
			temp_transcript_tpm = pd.DataFrame()
			temp_transcript_count = pd.DataFrame()
		temp_transcript_tpm = temp_transcript_tpm.append(transcript_tpm.loc[[k]])
		temp_transcript_count = temp_transcript_count.append(transcript_count.loc[[k]])

# Add the last row
if (sum(gene_tpm.loc[cur_gene] >= 0.1) >= 10):
	filtered_gene_list.append(cur_gene)
	filtered_transcript_list = filtered_transcript_list + temp_transcript_tpm.index.tolist()
	gene_tpm_filter = gene_tpm_filter.append(gene_tpm.loc[[cur_gene]])
	gene_count_filter = gene_count_filter.append(gene_count.loc[[cur_gene]])
	transcript_tpm_filter = transcript_tpm_filter.append(temp_transcript_tpm)
	transcript_count_filter = transcript_count_filter.append(temp_transcript_count)

print('Finished filtering genes/transcripts')

f = open(exp_dir + 'gene_list/' + tissue + '_filtered_gene_list.txt', 'w')
f.write('\n'.join(filtered_gene_list))
f.close()
f = open(exp_dir + 'gene_list/' + tissue + '_filtered_transcript_list.txt', 'w')
f.write('\n'.join(filtered_transcript_list))
f.close()


# For QN across the entire expression profile:
ranks = rankdata(gene_tpm_filter)
table = pd.DataFrame(ranks.reshape(len(gene_tpm_filter.index), len(gene_tpm_filter.columns)), gene_tpm_filter.index, gene_tpm_filter.columns)

# Now run QN for each gene, breaking ties randomly:
for x in table.index:
	# Add a random vector, smaller than the minimum precision, to randomize tie breaking:
	rand_vector = pd.Series(np.random.choice(range(1,table.shape[1]+1),table.shape[1],replace=False) * 0.0005)
	ranks = rankdata(np.array(table.loc[x]) + rand_vector)
	new_qn = norm.ppf(ranks / (table.shape[1] + 1))
	table.loc[x] = new_qn

# Write the normalized gene tpm table:
table.to_csv(exp_dir + 'quantile_norm/v8_RSEMv1.3.0_gene_tpm_' + tissue + '_profile_norm.txt', sep='\t', float_format='%.4f')

# Repeat for transcript table:
ranks = rankdata(transcript_tpm_filter)
table = pd.DataFrame(ranks.reshape(len(transcript_tpm_filter.index), len(transcript_tpm_filter.columns)), transcript_tpm_filter.index, transcript_tpm_filter.columns)

# Now run QN for each gene, breaking ties randomly:
for x in table.index:
	# Add a random vector, smaller than the minimum precision, to randomize tie breaking:
	rand_vector = pd.Series(np.random.choice(range(1,table.shape[1]+1),table.shape[1],replace=False) * 0.0005)
	ranks = rankdata(np.array(table.loc[x]) + rand_vector)
	new_qn = norm.ppf(ranks / (table.shape[1] + 1))
	table.loc[x] = new_qn

# Write the normalized gene tpm table:
table.to_csv(exp_dir + 'quantile_norm/v8_RSEMv1.3.0_transcript_tpm_' + tissue + '_profile_norm.txt', sep='\t', float_format='%.4f')


# Also try QN for each sample, then for gene:
table = gene_tpm_filter.copy()

for sample in table.columns:
	table[sample] = rankdata(table[sample])

# Now run QN for each gene, breaking ties randomly:
for x in table.index:
	# Add a random vector, smaller than the minimum precision, to randomize tie breaking:
	rand_vector = pd.Series(np.random.choice(range(1,table.shape[1]+1),table.shape[1],replace=False) * 0.0005)
	ranks = rankdata(np.array(table.loc[x]) + rand_vector)
	new_qn = norm.ppf(ranks / (table.shape[1] + 1))
	table.loc[x] = new_qn

# Write the normalized gene tpm table:
table.to_csv(exp_dir + 'quantile_norm/v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt', sep='\t', float_format='%.4f')

# Repeat for transcript table:
table = transcript_tpm_filter.copy()

for sample in table.columns:
	table[sample] = rankdata(table[sample])

# Now run QN for each gene, breaking ties randomly:
for x in table.index:
	# Add a random vector, smaller than the minimum precision, to randomize tie breaking:
	rand_vector = pd.Series(np.random.choice(range(1,table.shape[1]+1),table.shape[1],replace=False) * 0.0005)
	ranks = rankdata(np.array(table.loc[x]) + rand_vector)
	new_qn = norm.ppf(ranks / (table.shape[1] + 1))
	table.loc[x] = new_qn

# Write the normalized gene tpm table:
table.to_csv(exp_dir + 'quantile_norm/v8_RSEMv1.3.0_transcript_tpm_' + tissue + '_sample_norm.txt', sep='\t', float_format='%.4f')

print('Finished normalizing TPM matrices')

# Now log-transform the count matrices:
np.log2(gene_count_filter+1).to_csv(exp_dir + 'log_transform/v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt', sep='\t', float_format='%.4f')
np.log2(transcript_count_filter+1).to_csv(exp_dir + 'log_transform/v8_RSEMv1.3.0_transcript_count_' + tissue + '_log_transform.txt', sep='\t', float_format='%.4f')

print('Finished log-transforming count matrices')