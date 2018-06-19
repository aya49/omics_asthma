
# This file is obsolete - going to directly compare with the SNP IDs

import glob

gemma_cis_eqtls = '/tigress/BEE/amish/analyses/ciseqtl/genomewide/'
# gemma_cis_eqtls = '/tigress/BEE/RNAseq/Output/cis-mapping/amish/temp/'
files = glob.glob(gemma_cis_eqtls + '*')
out_file_prefix = '/tigress/BEE/RNAseq/Output/cis-mapping/amish/cis_eqtls_150kb_chr'
for f in files:
	# Get chrom ID
	chrom = f.split('chr')[1][:-4]
	# Location of SNP data files
	snp_metadata_f = open('/tigress/BEE/RNAseq/Data/Genotype/gtex/SNP_metadata/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr' + chrom + '.metadata.txt', 'r')
	while True:
		# This breaks out when it reads the header line
		if snp_metadata_f.readline()[0:2] != '##':
			break
	# Populate the rsid dictionary for the chromosome
	rsid_dict = {}
	for line in snp_metadata_f.readlines():
		entry = line.strip().split('\t')
		rsid_dict[entry[2]] = entry[4]
	# Now go through the cis-eQTL list and write those that have a corresponding GTEx rsID:
	out_f = open(out_file_prefix + chrom + '.txt', 'w')
	ciseqtl_f = open(f, 'r')
	header = ciseqtl_f.readline()
	new_header = ['rsid', 'gene_id', 'snpid'] + header.split('\t')[7:]
	out_f.write('\t'.join(new_header))
	for line in ciseqtl_f.readlines():
		# Check if in rsid dictionary:
		entry = line.split('\t')
		if entry[2] in rsid_dict:
			# Write the line
			new_entry = [rsid_dict[entry[2]], entry[0], entry[2]] + entry[7:]
			out_f.write('\t'.join(new_entry))
	out_f.close()
