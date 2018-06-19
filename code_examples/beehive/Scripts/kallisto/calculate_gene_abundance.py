
# calculate_gene_abundance.py
# This script prepares the three tables with necessary information (sample-level, subject-level and tissue-level) for downstream processing.
# Author: Brian Jo
# Date: 10/28/16

gencode_file = open('/tigress/BEE/eQTLs/Data/References/Annotations/gencode.v19.annotation.gtf')

for i in range(5):
	gencode_file.readline()

# create inverted index

transcript_dict = {}
for line in gencode_file.readlines():
	entry = str.split(line.strip(), '\t')
	if entry[2] == 'transcript':
		annotations = str.split(entry[8], '"; ')
		gene_id = annotations[0][9:]
		transcript_id = annotations[1][15:]
		transcript_dict[transcript_id] = 

gencode_file.close()
# Finished making inverted index

input_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125/abundance.tsv')
output_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125/gene_abundance.tsv', 'w')
output_file2 = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125/misc_abundance.tsv', 'w')

input_file.readline()

gene_abundance_dict = {}
gene_not_found_dict = {}
for line in input_file.readlines():
	entry = str.split(line.strip(), '\t')
	tpm = entry[4]
	transcript_id = entry[0]
	if float(tpm) > 0:
		if transcript_id in transcript_dict:
			gene_id = transcript_dict[transcript_id]
			if gene_id not in gene_abundance_dict:
				gene_abundance_dict[gene_id] = 0
			gene_abundance_dict[gene_id] = gene_abundance_dict[gene_id] + float(tpm)
		else:
			gene_not_found_dict[transcript_id] = float(tpm)

# write the gene quants
for item in gene_abundance_dict:
	output_file.write(item + '\t' + ("%.6f" % gene_abundance_dict[item]) + '\n');

for item in gene_not_found_dict:
	output_file2.write(item + '\t' + ("%.6f" % gene_not_found_dict[item]) + '\n');

output_file.close()
output_file2.close()

gene_bs_abundance_dict = {}

for i in range(100):
	print(i)
	input_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125/bootstrap/bs_abundance_'+str(i)+'.tsv')
	input_file.readline()
	for line in input_file.readlines():
		entry = str.split(line.strip(), '\t')
		tpm = entry[4]
		transcript_id = entry[0]
		if float(tpm) > 0:
			if transcript_id in transcript_dict:
				gene_id = transcript_dict[transcript_id]
				if gene_id not in gene_bs_abundance_dict:
					gene_bs_abundance_dict[gene_id] = {}
				if str(i) not in gene_bs_abundance_dict[gene_id]:
					gene_bs_abundance_dict[gene_id][str(i)] = 0
				gene_bs_abundance_dict[gene_id][str(i)] = gene_bs_abundance_dict[gene_id][str(i)] + float(tpm)
	input_file.close()

output_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125/bootstrap/gene_abundance_table.tsv', 'w')

for item in gene_bs_abundance_dict:
	output_file.write(item);
	for i in range(99):
		if str(i) in gene_bs_abundance_dict[item]:
			output_file.write('\t' + ("%.6f" % gene_bs_abundance_dict[item][str(i)]));
		else:
			output_file.write('\t' + "0.0");
	output_file.write('\n');

output_file.close()


# repeat for _ext


input_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125_ext/abundance.tsv')
output_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125_ext/gene_abundance.tsv', 'w')
output_file2 = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125_ext/misc_abundance.tsv', 'w')

input_file.readline()

gene_abundance_dict = {}
gene_not_found_dict = {}
for line in input_file.readlines():
	entry = str.split(line.strip(), '\t')
	tpm = entry[4]
	transcript_id = entry[0]
	if float(tpm) > 0:
		if transcript_id in transcript_dict:
			gene_id = transcript_dict[transcript_id]
			if gene_id not in gene_abundance_dict:
				gene_abundance_dict[gene_id] = 0
			gene_abundance_dict[gene_id] = gene_abundance_dict[gene_id] + float(tpm)
		else:
			gene_not_found_dict[transcript_id] = float(tpm)

# write the gene quants
for item in gene_abundance_dict:
	output_file.write(item + '\t' + ("%.6f" % gene_abundance_dict[item]) + '\n');

for item in gene_not_found_dict:
	output_file2.write(item + '\t' + ("%.6f" % gene_not_found_dict[item]) + '\n');

output_file.close()
output_file2.close()

gene_bs_abundance_dict = {}

for i in range(100):
	print(i)
	input_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125_ext/bootstrap/bs_abundance_'+str(i)+'.tsv')
	input_file.readline()
	for line in input_file.readlines():
		entry = str.split(line.strip(), '\t')
		tpm = entry[4]
		transcript_id = entry[0]
		if float(tpm) > 0:
			if transcript_id in transcript_dict:
				gene_id = transcript_dict[transcript_id]
				if gene_id not in gene_bs_abundance_dict:
					gene_bs_abundance_dict[gene_id] = {}
				if str(i) not in gene_bs_abundance_dict[gene_id]:
					gene_bs_abundance_dict[gene_id][str(i)] = 0
				gene_bs_abundance_dict[gene_id][str(i)] = gene_bs_abundance_dict[gene_id][str(i)] + float(tpm)
	input_file.close()

output_file = open('/tigress/BEE/RNAseq/Output/kallisto/SRR1476125_ext/bootstrap/gene_abundance_table.tsv', 'w')

for item in gene_bs_abundance_dict:
	output_file.write(item);
	for i in range(99):
		if str(i) in gene_bs_abundance_dict[item]:
			output_file.write('\t' + ("%.6f" % gene_bs_abundance_dict[item][str(i)]));
		else:
			output_file.write('\t' + "0.0");
	output_file.write('\n');

output_file.close()








