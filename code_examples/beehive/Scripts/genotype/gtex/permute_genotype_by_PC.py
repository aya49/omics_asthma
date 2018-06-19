##################################################
#  permute_genotype_by_PC.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/permute_genotype_by_PC.py
# 
#  This script permutes the genotypes for null distribution of data, by the first PC in order to decrease the FDR rate
#
#  Author: Brian Jo
#
##################################################

genotype_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous/'
genotype_perm_output_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous_permute_by_PC/'

import random
import pandas

genotype_pc_df = pandas.read_table('/tigress/BEE/eQTLs/Data/References/GTEx/GTEx_Genotype_PCs_top3.txt')
# Take only columns for which the second PC is over -0.3:
cols_to_take = genotype_pc_df.C2 > -0.3
cols_to_take_list = [genotype_pc_df.index[x] for x in range(450) if cols_to_take[x]]
# Population segregation by first genotype PC value of 0.04
pop_group_1 = genotype_pc_df.C1 > 0.04
pop_group_1_list = [genotype_pc_df.index[x] for x in range(450) if pop_group_1[x]]
pop_group_2_list = [genotype_pc_df.index[x] for x in range(450) if cols_to_take[x] and not pop_group_1[x]]

# Only for autosomes
for i in range(22):
    print(i)
    in_f = open(genotype_dir + 'GTEx_genotypes_maf05_continuous_Chr' + str(i+1) + '_Final.txt')
    out_f = open(genotype_perm_output_dir + 'GTEx_genotypes_maf05_continuous_Chr' + str(i+1) + '_Final.txt', 'w')
    header = in_f.readline()
    sample_list = str.split(header.strip(), '\t')
    pop_group_1_ind = [x for x in range(450) if sample_list[x] in pop_group_1_list]
    pop_group_2_ind = [x for x in range(450) if sample_list[x] in pop_group_2_list]
    out_f.write('\t'.join(sample_list[x] for x in pop_group_1_ind) + '\t' + '\t'.join(sample_list[x] for x in pop_group_2_ind))
    out_f.write('\n')
    for line in in_f.readlines():
        entry = str.split(line.strip(), '\t')
        genotypes = entry[1:]
        random.shuffle(pop_group_1_ind)
        random.shuffle(pop_group_2_ind)
        out_f.write(entry[0] + '\t' + '\t'.join(genotypes[x] for x in pop_group_1_ind) + '\t' + '\t'.join(genotypes[x] for x in pop_group_2_ind))
        out_f.write('\n')
    out_f.close()
    in_f.close()