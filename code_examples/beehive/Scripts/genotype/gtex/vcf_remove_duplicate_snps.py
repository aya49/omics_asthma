##################################################
#  vcf_remove_duplicate_snps.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/vcf_remove_duplicate_snps.py
# 
#  It was discovered that at least one SNP is duplicated in one imputed genotype file, this script ensures that
#  there are no duplicates in all genotype files (duplicates are removed as it is unreasonable to speculate which is "correct").
#
#  Author: Ian McDowell
#
##################################################

from sys import argv
import gzip
from collections import OrderedDict

in_file = argv[1]
out_file = argv[2]
report = argv[3]

in_file = gzip.open(in_file, 'rb')
out_file = gzip.open(out_file, 'wb')
report = open(report, 'w')

duplicate_counter = 0

out_lines_dict = OrderedDict()
snp_set = set()

for line in in_file:
    if line[0] =='#': # print VCF header lines to output directly
        out_file.write(line)
        continue
    snp = line.split()[2] # parse out snp name
    if snp in snp_set: # if SNP is already in snp_set that means it is a duplicate
        del out_lines_dict[snp] # remove duplicate snp from lines to be written
        duplicate_counter += 1 # count duplicates
        report.write(snp + '\n')
    else:
        snp_set.add(snp)
        out_lines_dict[snp] = line

for v in out_lines_dict.values():
    out_file.write(v)

print 'Number of duplicates =', duplicate_counter
in_file.close()
out_file.close()
report.close()