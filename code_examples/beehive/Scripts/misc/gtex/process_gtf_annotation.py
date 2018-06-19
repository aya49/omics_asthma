##################################################
#  process_gtf_annotation.py
#
#  $proj/Scripts/misc/gtex/process_gtf_annotation.py
# 
#  This script reorganizes gtf file into R format for easy loading
#
#  Author: Brian Jo
#
##################################################

import os
proj_dir = os.environ['proj']

def process_gtf_entry(entry):
	info = {}
	annotations = entry.strip().split("; ")
	for item in annotations:
		item = item.strip(';')
		item = item.strip(' ')
		kv_pair = item.split(" ")
		(key, value) = (kv_pair[0], kv_pair[1])
		value = value.strip('"')
		info[key] = value
	return info

gtf_file = open(proj_dir + '/Data/Resources/annotations/gencode.v19.annotation.gtf')

for i in range(5):
	gtf_file.readline()

orig_fields = ['chr', 'src', 'type', 'start', 'end', 'strand']

# Split up by chromosome
gene_entries = {}
exon_entries = {}
transcript_entries = {}

cur_chr = '1'
gene_entries[cur_chr] = []
exon_entries[cur_chr] = []
transcript_entries[cur_chr] = []

for line in gtf_file.readlines():
	entry = line.split('\t')
	if entry[0][3:] != cur_chr:
		cur_chr = entry[0][3:]
		print(cur_chr)
		gene_entries[cur_chr] = []
		exon_entries[cur_chr] = []
		transcript_entries[cur_chr] = []
	anno = process_gtf_entry(entry[8])
	anno['chr'] = cur_chr
	anno['src'] = entry[1]
	anno['type'] = entry[2]
	anno['start'] = entry[3]
	anno['end'] = entry[4]
	anno['strand'] = entry[6]
	if entry[2] == 'gene':
		gene_entries[cur_chr].append(anno)
	if entry[2] == 'exon':
		exon_entries[cur_chr].append(anno)
	if entry[2] == 'transcript':
		transcript_entries[cur_chr].append(anno)

gtf_file.close()

gene_fields = orig_fields + ['gene_id', 'gene_type', 'gene_status', 'gene_name']
transcript_fields = gene_fields + ['transcript_id', 'transcript_type', 'transcript_status', 'transcript_name']
exon_fields = transcript_fields + ['exon_id', 'exon_number']

# Output annotation files
out_dir = proj_dir + '/Data/Expression/gene_metadata_hg19/'
chr_list = [str(i+1) for i in range(22)] + ['X', 'Y', 'M']
for cur_chr in chr_list:
	print(cur_chr)
	# write gene file
	out_file = out_dir + 'gene_metadata_chr' + cur_chr + '.txt'
	f = open(out_file, 'w')
	# write header
	f.write('\t'.join(gene_fields) + '\n')
	annotation_list = gene_entries[cur_chr]
	for item in annotation_list:
		f.write('\t'.join([item[x] for x in gene_fields]) + '\n')
	# write exon file
	f.close()
	out_file = out_dir + 'transcript_metadata_chr' + cur_chr + '.txt'
	f = open(out_file, 'w')
	# write header
	f.write('\t'.join(transcript_fields) + '\n')
	annotation_list = transcript_entries[cur_chr]
	for item in annotation_list:
		f.write('\t'.join([item[x] for x in transcript_fields]) + '\n')
	# write gene file
	f.close()
	out_file = out_dir + 'exon_metadata_chr' + cur_chr + '.txt'
	f = open(out_file, 'w')
	# write header
	f.write('\t'.join(exon_fields) + '\n')
	annotation_list = exon_entries[cur_chr]
	for item in annotation_list:
		f.write('\t'.join([item[x] for x in exon_fields]) + '\n')
	f.close()
