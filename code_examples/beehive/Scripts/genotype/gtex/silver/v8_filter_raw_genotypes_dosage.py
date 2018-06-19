##################################################
#  v8_filter_raw_genotypes_dosage.py
#
#  $proj/Scripts/genotype/gtex/silver/v8_filter_raw_genotypes_dosage.py
# 
#  This script process the raw WGS vcf file and splits them into chromosomes, while filtering for MAF and indels
#
#  Author: Brian Jo
#
##################################################

# function for processing each line
def process_vcf_line(line):
	# We will apply the two following filters:
	# Filter based on MAF
	# Only take SNPs - filter out indels
	entry = line.strip().decode("utf-8").split('\t')
	processed_line = {}
	processed_line['chr'] = entry[0][3:]
	processed_line['pos'] = entry[1]
	# don't process for MAF < 0.01
	processed_line['MAF'] = float(entry[7].split(';')[1][3:])
	processed_line['ID'] = entry[2]
	if processed_line['MAF'] < 0.01:
		return processed_line
	processed_line['SNP'] = True
	# is it an indel?
	if len(entry[3]) > 1 or len(entry[4]) > 1:
		processed_line['SNP'] = False
	dosage_arr = [str(int(entry[x][0]) + int(entry[x][2])) if entry[x][0] != '.' else '-' for x in range(9, len(entry))]
	processed_line['dosage'] = dosage_arr
	return processed_line

import gzip
import os

vcf_raw = '/tigress/BEE/gtex/dbGaP-7716/57610/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz'
f = gzip.open(vcf_raw, 'rb')

while True:
	line = f.readline()
	# traverse until the suject IDs are found
	if line.strip().decode("utf-8").split('\t')[0] == '#CHROM':
		break

# This line contains the subject IDs
vcf_header = line.strip().decode("utf-8")
subject_IDs = vcf_header.split('\t')[9:]
# len(subject_IDs) = 838
chrom_list = [str(x+1) for x in range(22)] + ['X', 'Y', 'M']

# We will make two files - for MAF 5% and 1%
out_dir = '/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/'
current_chr = ''
prev_pos = '0'
prev_id = ''
multiallelic_list = []
# dummy files
dosage_out_f_01 = open(out_dir + 'allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_01.txt', 'w')
dosage_out_f_05 = open(out_dir + 'allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_05.txt', 'w')
vcf_out_f_01 = open(out_dir + 'vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_MAF_01.txt', 'w')
vcf_out_f_05 = open(out_dir + 'vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_MAF_05.txt', 'w')
multiallele_f = open(out_dir + 'multiallelic/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_multiallelic.txt', 'w')
indel_f = open(out_dir + 'indels/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_indels.txt', 'w')

for line in f:
	processed_line = process_vcf_line(line)
	if processed_line['chr'] not in chrom_list:
		continue
	# Multiallelic?
	if processed_line['pos'] == prev_pos:
		multiallelic_list = multiallelic_list + [prev_id, processed_line['ID']]
	prev_pos = processed_line['pos']
	prev_id = processed_line['ID']
	if current_chr != processed_line['chr']:
		# Make new set of files
		current_chr = processed_line['chr']
		print(current_chr)
		# dosage files
		dosage_out_f_01.close()
		dosage_out_f_05.close()
		dosage_out_f_01 = open(out_dir + 'allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_01.txt', 'w')
		dosage_out_f_05 = open(out_dir + 'allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_05.txt', 'w')
		# write header
		dosage_out_f_01.write('\t' + '\t'.join(subject_IDs) + '\n')
		dosage_out_f_05.write('\t' + '\t'.join(subject_IDs) + '\n')
		# vcf files
		vcf_out_f_01.close()
		vcf_out_f_05.close()
		vcf_out_f_01 = open(out_dir + 'vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_vcf_MAF_01.txt', 'w')
		vcf_out_f_05 = open(out_dir + 'vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_vcf_MAF_05.txt', 'w')
		# write header
		vcf_out_f_01.write(vcf_header + '\n')
		vcf_out_f_05.write(vcf_header + '\n')
		# Other files:
		multiallelic_set = set(multiallelic_list)
		multiallele_f.write('\n'.join(sorted(multiallelic_set)))
		multiallele_f.close()
		multiallele_f = open(out_dir + 'multiallelic/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_multiallelic.txt', 'w')
		multiallelic_list = []
		indel_f.close()
		indel_f = open(out_dir + 'indels/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_indels.txt', 'w')
	if processed_line['MAF'] >= 0.05:
		dosage_out_f_05.write(processed_line['ID'] + '\t' + '\t'.join(processed_line['dosage']) + '\n')
		vcf_out_f_05.write(line.strip().decode("utf-8") + '\n')
	elif processed_line['MAF'] >= 0.01:
		dosage_out_f_01.write(processed_line['ID'] + '\t' + '\t'.join(processed_line['dosage']) + '\n')
		vcf_out_f_01.write(line.strip().decode("utf-8") + '\n')
	else:
		continue
	if not processed_line['SNP']:
		indel_f.write(processed_line['ID'] + '\n')

# Close the remaining files
dosage_out_f_01.close()
dosage_out_f_05.close()
vcf_out_f_01.close()
vcf_out_f_05.close()
multiallelic_set = set(multiallelic_list)
multiallele_f.write('\n'.join(sorted(multiallelic_set)))
multiallele_f.close()
indel_f.close()

# remove dummy files
os.remove(out_dir + 'allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_dosage_MAF_01.txt')
os.remove(out_dir + 'allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_dosage_MAF_01.txt')
os.remove(out_dir + 'vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_vcf_MAF_01.txt')
os.remove(out_dir + 'vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_vcf_MAF_05.txt')
os.remove(out_dir + 'multiallelic/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_multiallelic.txt')
os.remove(out_dir + 'indels/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_indels.txt')
