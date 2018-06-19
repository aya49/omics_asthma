##################################################
#  obtain_imputed_genotype.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/obtain_imputed_genotype.py
# 
#  This script takes the vcf files and returns the allelic dosage format file in continuous genotypes
#
#  Author: Brian Jo
#
##################################################

from sys import argv
import gzip

vcf_dir = argv[1]
output_dir = argv[2]
chr_num = argv[3]

# Take the imputed genotypes directly from vcf files
f = gzip.open(vcf_dir + "GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput_info04_maf05_HWEp1E6_dbSNPIDs.chr" + chr_num + ".vcf.gz")
line = f.readline()
out_f = open(output_dir + "genotypesChr" + chr_num + "_temp.txt", 'w')
    
while True:
    line = f.readline()
    if line[0:2] != '##':
        break
            
# We went through the annotations, now the header line:
listHead = str.split(line.strip(), "\t")
# First 9 entries are other annotations
for j in range(9, len(listHead)):
    listHead[j] = '.'.join(str.split(listHead[j], "-")[0:2])
        
out_f.write('\t'.join(listHead[9:len(listHead)]))
    
count = 0
for line in f.readlines():
    count = count+1
    genotype_list = []
    listGen = str.split(line.strip(), "\t")
    # Add in rs SNP ID
    genotype_list.append(listGen[2])
    for k in range(9, len(listGen)):
        imput = str.split(listGen[k], ":")
        # Fixed version - just take the decimal value
        genotype_list.append(imput[2])
        # # Confidence under 90% - take the decimal form
        # if imput[0] == './.':
        #     genotype_list.append(imput[2])
        # # Confidence over 90% - take 0, 1 or 2
        # else:
        #     genotype_list.append(str(int(listGen[k][0]) + int(listGen[k][2])))
    if count%10000 == 0:
        print count
    out_f.write('\n' + '\t'.join(genotype_list))
        
f.close()
out_f.close()