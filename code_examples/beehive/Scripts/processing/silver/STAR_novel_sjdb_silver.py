##################################################
#  STAR_novel_sjdb_silver.py
#
#  $proj/Scripts/processing/silver/STAR_novel_sjdb_silver.py
#
#  For the creation of novel splice junction database based on completed STAR 1-pass samples
#
#  Authors: Ian McDowell, Brian Jo
#
##################################################

import glob
import pandas as pd
from collections import defaultdict, OrderedDict
import pickle
import matplotlib

matplotlib.use('Agg')
font = {'size'   : 8}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt

from sys import argv

base = argv[1]
out_splice_dict = argv[2]
out_plot = argv[3]
out_sjdb = argv[4]

# example
# base = '/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/STAR_1pass_hg38/'
# out_splice_dict = '/tigress/BEE/gtex/data/silver/STAR_hg38_1pass/STAR_1pass_sjdb.p'
# out_plot = '/tigress/BEE/gtex/data/silver/STAR_hg38_1pass/Number_of_splice_sites_by_number_of_samples_hg19.png'
# out_sjdb = '/tigress/BEE/gtex/data/silver/STAR_hg38_1pass/SJ.out.tab.Pass1.one.percent.sjdb'

# find all novel splice junctions
finished_file = '/Log.final.out'
suffix = '/SJ.out.tab'

finished_files = glob.glob(base + '*' + finished_file)
SJ_DBs = [s.replace(finished_file, suffix) for s in finished_files]

# remove base and suffix to get list of samples
samples = [s.replace(base, '').replace(finished_file, '') for s in finished_files]

# update the list of samples used to create novel splice junction database
star_1pass_samples = '/tigress/BEE/RNAseq/Scripts/processing/silver/checkpoints/STAR_1pass_hg38_samples.txt'
with open(star_1pass_samples, 'w') as f:
    f.write('\n'.join(samples) + '\n')

strand_dict = {}
strand_dict[0],strand_dict[1],strand_dict[2] = '.','+','-'

# SJ.out.tab format:
    
# Column 1: chromosome
# Column 2: first base of the intron (1-based)
# Column 3: last base of the intron (1-based)
# Column 4: strand
# Column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
# Column 6: reserved
# Column 7: number of uniquely mapping reads crossing the junction
# Column 8: reserved
# Column 9: maximum left/right overhang 

splice_dict = defaultdict(dict)
columns = ['chrom','start','end','strand','intron_motif','reserved','num_uniq_mapped_reads','reserved','left/right_overhang']
for i, (sample, SJ_DB) in enumerate(zip(samples, SJ_DBs)):
    print(i, sample)
    SJ_DB = pd.read_csv(SJ_DB, delim_whitespace=True, header=None, names = columns)
    #SJ_DB = SJ_DB[SJ_DB.chrom.isin(canonical_chroms)]
    # filter out non-canonical splice junctions
    SJ_DB = SJ_DB[SJ_DB.intron_motif != 0]
    # filter out splice juntions without support of greater than 5 reads
    SJ_DB = SJ_DB[SJ_DB.num_uniq_mapped_reads > 5]
    SJ_DB['strand'] = [strand_dict[strand] for strand in SJ_DB['strand']]
    for chrom, start, end, strand in zip(SJ_DB.chrom, SJ_DB.start, SJ_DB.end, SJ_DB.strand):
        if (chrom,start,end,strand) in splice_dict:
            splice_dict[(chrom,start,end,strand)] += 1
        else:
            splice_dict[(chrom,start,end,strand)] = 1

pickle.dump(splice_dict, open(out_splice_dict, 'wb'))
# splice_dict = pickle.load(open('/tigress/BEE/gtex/data/auxiliary/STAR_hg19/STAR_1pass_sjdb.p', 'rb'))

across_libraries_dict = {}
final_cutoff = int(len(samples) * 0.01)
# final_cutoff = 1
# print('Number of samples splice site required to be present, final =',final_cutoff)

cutoffs = list(range(1,10)) + list(range(10,1010,10)) + [final_cutoff]
cutoffs.sort()

for cutoff in cutoffs:
    # print 'Number of samples splice site required to be present =', cutoff
    num = sum([1 if v >= cutoff else 0 for v in splice_dict.values()])
    # print 'Number of splice sites at sample number cutoff =', num
    across_libraries_dict[cutoff] = num

od = OrderedDict(sorted(across_libraries_dict.items()))

# x = list(across_libraries_dict.keys())
# x.sort()
# y = list(across_libraries_dict.values())
# y.sort(reverse=True)

plt.figure()
plt.plot(list(across_libraries_dict.keys()), list(across_libraries_dict.values()), "o")
# plt.axvline(10, c='red')
# plt.axvline(20, c='red')
# plt.axvline(50, c='red')
plt.axvline(final_cutoff, c='red')
# plt.annotate('10', (10,across_libraries_dict[10]))
# plt.annotate('20', (20,across_libraries_dict[20]))
# plt.annotate('50', (50,across_libraries_dict[50]))
plt.annotate('%s'%(final_cutoff), (final_cutoff,across_libraries_dict[final_cutoff]))
plt.xlabel('Minimum number of samples in which novel splice site is present')
plt.ylabel('Number of novel splice sites')
plt.savefig(out_plot)

# It looks like splice site presence across 1% of samples is a good compromise between conservative and permissive
splice_sites = [k for k,v in splice_dict.items() if v >= final_cutoff ]
# canonical_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
# splice_sites = [splice_site for splice_site in splice_sites if splice_site[0] in canonical_chroms]
# splice_sites = sorted(splice_sites, key=lambda s: int(s[1]))
# splice_sites = sorted(splice_sites, key=lambda s: int(s[0][3:]))
splice_sites = ['\t'.join([str(x) for x in splice_site]) for splice_site in splice_sites]

out = open(out_sjdb,'w')
out.write('\n'.join(splice_sites) + '\n')
out.close()
