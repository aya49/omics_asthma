f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_75mer.wig')
out_f = open('/tigress/BEE/RNAseq/Output/processing/mappability/hg38_75mer.wig_mod', 'w')

for line in f.readlines():
	if '  AC' in line:
		new = line.replace('  AC', '')
		out_f.write(new)
	else:
		out_f.write(line)
