##################################################
#  vcf_translate_names.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_translate_names.py
# 
#  This script will translate GTEx custom names to conventional dbSNP (v.142) names
#
#  Author: Ian McDowell
#
##################################################

#!/usr/bin/env python
import pickle
from sys import argv
import gzip

in_file = argv[1]
SNP_dict = argv[2]
out_file = argv[3]

in_file = gzip.open(in_file, 'rb')
out_file = gzip.open(out_file, 'wb')

SNP_dict = pickle.load(open(SNP_dict, 'rb'))

for line in in_file:
    if "#" not in line:
        split_line = line.split('\t')
        newID = SNP_dict[split_line[2]]
        new_split_line = split_line[0:2] + [newID] + split_line[3:]
        line = '\t'.join(new_split_line)
    out_file.write(line)

out_file.close()
in_file.close()