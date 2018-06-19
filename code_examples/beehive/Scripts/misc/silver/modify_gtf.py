##################################################
#  modify_gtf.py
#
#  $proj/Scripts/misc/silver/modify_gtf.py
#
#  For modifying the gtf files to have consistent chr names with new genome fasta file
#
#  Author: Brian Jo
#
##################################################

in_file = '/tigress/BEE/gtex/external_sources/hg38_silver/gencode.v26.annotation.gtf'
out_file = '/tigress/BEE/gtex/external_sources/hg38_silver/gencode.v26.annotation.mod.gtf'

in_f = open(in_file, 'r')
out_f = open(out_file, 'w')

for i in range(5):
	line = in_f.readline()
	out_f.write(line)

for line in in_f.readlines():
	out_f.write(line[3:])

in_f.close()
out_f.close()