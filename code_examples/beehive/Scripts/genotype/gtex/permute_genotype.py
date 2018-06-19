##################################################
#  permute_genotype.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/permute_genotype.py
# 
#  This script permutes the genotypes for null distribution of data
#
#  Author: Brian Jo
#
##################################################

genotype_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous/'
genotype_perm_output_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous_permute/'

import random

# Only for autosomes
# Can also run this for maf01 if needed
for i in range(22):
    print(i)
    in_f = open(genotype_dir + 'GTEx_genotypes_maf05_continuous_Chr' + str(i+1) + '_Final.txt')
    out_f = open(genotype_perm_output_dir + 'GTEx_genotypes_maf05_continuous_Chr' + str(i+1) + '_Final.txt', 'w')
    header = in_f.readline()
    out_f.write(header)
    for line in in_f.readlines():
        entry = str.split(line.strip(), '\t')
        genotypes = entry[1:]
        random.shuffle(genotypes)
        out_f.write(entry[0] + '\t' + '\t'.join(genotypes) + '\n')
    out_f.close()
    in_f.close()