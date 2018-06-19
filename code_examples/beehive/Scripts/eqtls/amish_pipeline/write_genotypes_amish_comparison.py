
import gzip
import pandas as pd
from sys import argv

chrom = argv[1]

gemma_cis_eqtls = '/tigress/BEE/amish/analyses/ciseqtl/genomewide/cis_eqtls_150kb_chr'
# gemma_cis_eqtls = '/tigress/BEE/RNAseq/Output/cis-mapping/amish/temp/'
f = gemma_cis_eqtls + chrom + '.txt'
out_file_prefix = '/tigress/BEE/RNAseq/Data/Genotype/gtex/amish_comparison/chr'

# Get chrom ID
amish_ciseqtls = pd.read_csv(f, sep='\t')
snps = pd.Series(amish_ciseqtls['rs']).unique()
# Go through the original vcf file to prepare the genotype file:
chr_f = gzip.open('/tigress/BEE/RNAseq/Data/Genotype/gtex/imputed_genotypes/vcf/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.chr' + chrom + '.vcf.gz', 'rb')
while True:
	# This breaks out when it reads the header line
	line = chr_f.readline()
	if line.decode('utf-8')[0:2] != '##':
		header = line.decode('utf-8')
		break
# Parse header
samps = [x.split('-')[0] + '-' + x.split('-')[1] for x in header.strip().split('\t')[9:]]
out_f = open(out_file_prefix + chrom + '.txt', 'w')
# write header
out_f.write('\t'.join(samps) + '\n')
for line in chr_f.readlines():
	# Check if the SNP needs to be written
	entry = line.decode('utf-8').strip().split('\t')
	if entry[2] in snps:
		# Parse line
		genotypes = [x.split(':')[-1] for x in entry[9:]]
		out_f.write(entry[2] + '\t' + '\t'.join(genotypes) + '\n')
out_f.close()
