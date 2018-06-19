##################################################
#  plink_recode_A_to_allelic_dosage_format.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/plink_recode_A_to_allelic_dosage_format.py
# 
#  This script will convert the genotypes from vcf format to allelic dosage format for eQTL analyses
#
#  Author: Ian McDowell
#
##################################################

import pandas as pd
from sys import argv

inf = argv[1]
outf = argv[2]

outfh = open(outf, 'w')

with open(inf, 'r') as f:
    for i, line in enumerate(f):
        if i==0: #i.e. header
            outfh.write('SNP' + line[3:])
        elif i in (1,2,3,4,5):
            pass # skip IID, PAT, MAT, SEX, PHENOTYPE lines
        else:
            split_line = line.split('\t')
            line = '\t'.join([split_line[0].split('_')[0]] + split_line[1:])
            outfh.write(line)
        
outfh.close()        

# inf = 'imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf05_HWEp1E6_dbSNPIDs.chr22.raw'
# outf = 'imputed_genotypes/allelic_dosage/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf05_HWEp1E6_dbSNPIDs.chr22.allelic_dosage.txt'

# g = pd.read_csv(inf, sep=' ')
# g['SNP'] = g['FID']
# g = g.T
# g.columns = g.ix['SNP']
# g = g.drop(['SNP','FID','IID','PAT','MAT','SEX','PHENOTYPE'], axis=0)
# g.index = [SNP.split('_')[0] for SNP in g.index] 
# g.to_csv(outf, float_format='%.0f', sep='\t')