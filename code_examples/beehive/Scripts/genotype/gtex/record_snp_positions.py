##################################################
#  record_snp_positions.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/record_snp_positions.py
# 
#  This script saves the SNP positions in a separate sep of files for quick loading
#
#  Author: Brian Jo
#
##################################################

metadata_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/SNP_metadata/'
output_dir = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/SNP_positions_hg19/'

out_file_all = output_dir + 'SNP_positions_ChrAll.txt'
out_f_all = open(out_file_all, 'w')
out_f_all.write('rsID\tchr\tpos\n')
for i in range(22):
	in_file = metadata_dir + 'GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr' + str(i+1) + '.metadata.txt'
	in_f = open(in_file)
	out_file = output_dir + 'SNP_positions_Chr' + str(i+1) + '.txt'
	out_f = open(out_file, 'w')
	out_f.write('rsID\tchr\tpos\n')
	for j in range(7):
		in_f.readline()
	SNP_dict = {}
	for line in in_f.readlines():
		entry = str.split(line.strip(), '\t')
		if entry[4] not in SNP_dict:
			SNP_dict[entry[4]] = []
			out_f.write(entry[4] + '\t' + 'chr' + entry[0] + '\t' + entry[1] + '\n')
			out_f_all.write(entry[4] + '\t' + 'chr' + entry[0] + '\t' + entry[1] + '\n')
	out_f.close()

out_f_all.close()