# f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_75mer.wig')
# out_f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_75mer.wig_mod', 'w')

# for line in f.readlines():
# 	if '  AC' in line:
# 		new = line.replace('  AC', '')
# 		out_f.write(new)
# 	else:
# 		out_f.write(line)

# f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_75mer.sizes')
# out_f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_75mer.sizes_mod', 'w')

# for line in f.readlines():
# 	if '  AC' in line:
# 		new = line.replace('  AC', '')
# 		out_f.write(new)
# 	else:
# 		out_f.write(line)

f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_36mer.wig')
out_f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_36mer.wig_mod', 'w')

for line in f.readlines():
	if '  AC' in line:
		new = line.replace('  AC', '')
		out_f.write(new)
	else:
		out_f.write(line)

f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_36mer.sizes')
out_f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_36mer.sizes_mod', 'w')

for line in f.readlines():
	if '  AC' in line:
		new = line.replace('  AC', '')
		out_f.write(new)
	else:
		out_f.write(line)