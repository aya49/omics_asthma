##################################################
#  process_gtf_annotation.py
#
#  $proj/Scripts/misc/silver/process_gtf_annotation.py
# 
#  This script reorganizes gtf file for easy loading
#
#  Author: Brian Jo
#
##################################################

import os
import pickle
from collections import OrderedDict

proj_dir = os.environ['proj']

# subroutine for processing a gtf annotation
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

# gtf_file = open(proj_dir + '/Data/Resources/annotations/gencode.v19.annotation.gtf')
gtf_file = open('/Users/brian_jo/Desktop/Project/scratch/gencode/gencode.v26.annotation.gtf')

for i in range(5):
	gtf_file.readline()

orig_fields = ['chr', 'src', 'type', 'start', 'end', 'strand']

# Split up by chromosome
gene_entries = OrderedDict()
exon_entries = OrderedDict()
transcript_entries = OrderedDict()

cur_chr = '1'
gene_entries[cur_chr] = OrderedDict()
exon_entries[cur_chr] = OrderedDict()
transcript_entries[cur_chr] = OrderedDict()

for line in gtf_file.readlines():
	entry = line.split('\t')
	if entry[0][3:] != cur_chr:
		cur_chr = entry[0][3:]
		print(cur_chr)
		gene_entries[cur_chr] = OrderedDict()
		exon_entries[cur_chr] = OrderedDict()
		transcript_entries[cur_chr] = OrderedDict()
	anno = process_gtf_entry(entry[8])
	anno['chr'] = cur_chr
	anno['src'] = entry[1]
	anno['type'] = entry[2]
	anno['start'] = entry[3]
	anno['end'] = entry[4]
	anno['strand'] = entry[6]
	if entry[2] == 'gene':
		gene_entries[cur_chr][anno['gene_id']] = anno
	if entry[2] == 'exon':
		exon_entries[cur_chr][anno['exon_id']] = anno
	if entry[2] == 'transcript':
		transcript_entries[cur_chr][anno['transcript_id']] = anno

gtf_file.close()

pickle.dump(gene_entries, open('/Users/brian_jo/Desktop/Project/scratch/gencode/gencode.v26.genes', 'wb'))
pickle.dump(exon_entries, open('/Users/brian_jo/Desktop/Project/scratch/gencode/gencode.v26.exons', 'wb'))
pickle.dump(transcript_entries, open('/Users/brian_jo/Desktop/Project/scratch/gencode/gencode.v26.transcripts', 'wb'))

# import pickle
# temp = pickle.load(open('/Users/brian_jo/Desktop/Project/scratch/gencode/gencode.v26.genes', 'rb'))
# for key, value in temp['M'].items():

