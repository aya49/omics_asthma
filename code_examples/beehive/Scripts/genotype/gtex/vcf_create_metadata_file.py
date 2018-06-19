##################################################
#  vcf_create_metadata_file.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_create_metadata_file.py
# 
#  This script creates the metadata files that store additional information about each genotype which may come in handy
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/env python
import pickle
from sys import argv
import gzip

in_file = argv[1]
SNP_dict = argv[2]
out_file = argv[3]

# import os
# os.chdir('/tigress/BEE/eQTLs/Data/Genotype/GTEx/')
# in_file = 'imputed_genotypes/vcf/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.chr1.vcf.gz'
# SNP_dict = 'auxiliary/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.chr1.dict'
# out_file = 'SNP_metadata/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr1.metadata.txt'

in_file = gzip.open(in_file, 'rb')
out_file = open(out_file, 'wb')

SNP_dict = pickle.load(open(SNP_dict, 'rb'))

while True:
    line = in_file.readline()
    if '#CHROM' in line:
        break
    elif '##INFO'in line:
        out_file.write(line)

# Define the header
metadata_elements = ['CHROM','POS','ID','IN_dbSNP','dbSNPID','REF','ALT','maf05_FILTER','A1_FREQ','IMPINFO','CERTAINTY','IMPUTED','MISS','HW_Pval']
header = '\t'.join(metadata_elements)
out_file.write(header)
out_file.write('\n')

for line in in_file.readlines():
    # Parse each line:
    split_line = line.split('\t')
    new_entry = split_line[0:3]
    if split_line[2] not in SNP_dict:
        new_entry += ['0','']
    else:
        new_entry += ['1',SNP_dict[split_line[2]]]

    new_entry += [split_line[3],split_line[4]]
    if split_line[6] == 'PASS':
        new_entry.append('1')
    elif split_line[6] == 'maf05':
        new_entry.append('0')
    else:
        print(split_line[6])
        raise ValueError('MAF Filter value should be either PASS or maf05!')

    new_entry += [x.split('=')[1] for x in split_line[7].split(';')]
    if new_entry[11] == '0':
        new_entry[11] = '1'
    elif new_entry[11] == '2':
        new_entry[11] = '0'
    else:
        print(new_entry)
        raise ValueError('IMPUTED value should be either 0 or 2!')

    out_file.write('\t'.join(new_entry))
    out_file.write('\n')

out_file.close()
in_file.close()