##################################################
#  v8_ld_pruned_subset_genotypes.py
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/silver/v8_ld_pruned_subset_genotypes.py
# 
#  This script saves a separate file of genotypes that are not LD-pruned
#
#  Author: Brian Jo
#
##################################################

import subprocess as sp
import os

vcf_directory = '/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/vcf/'
out_directory = '/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/ld_prune/'

os.chdir(out_directory)

for i in [str(x+1) for x in range(22)] + ['X']:
  # MAF 01
  input_f = vcf_directory + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + i + '_vcf_MAF_01.txt'
  sp.run(['/tigress/BEE/tools/plink', '--vcf', input_f, '--indep', '50', '5', '2'])
  sp.run(['mv', 'plink.log', out_directory+'chr'+i+'_MAF_01.log'])
  sp.run(['mv', 'plink.prune.in', out_directory+'chr'+i+'_MAF_01.in'])
  sp.run(['rm', 'plink.nosex'])
  sp.run(['rm', 'plink.prune.out'])
  # MAF 05
  input_f = vcf_directory + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + i + '_vcf_MAF_05.txt'
  sp.run(['/tigress/BEE/tools/plink', '--vcf', input_f, '--indep', '50', '5', '2'])
  sp.run(['mv', 'plink.log', out_directory+'chr'+i+'_MAF_05.log'])
  sp.run(['mv', 'plink.prune.in', out_directory+'chr'+i+'_MAF_05.in'])
  sp.run(['rm', 'plink.nosex'])
  sp.run(['rm', 'plink.prune.out'])

for i in [str(x+1) for x in range(22)] + ['X']:
  # MAF 01
  var_list = out_directory+'chr'+i+'_MAF_01.in'
  variants = [line.strip() for line in open(var_list).readlines()]
  input_f = vcf_directory + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + i + '_vcf_MAF_01.txt'
  output_f = out_directory + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + i + '_vcf_MAF_01_ld_pruned.txt'
  f = open(input_f)
  out_f = open(output_f, 'w')
  header = f.readline()
  out_f.write(header)
  for line in f.readlines():
    ID = line.split('\t')[2]
    if ID in variants:
      out_f.write(line)
  out_f.close()
  # MAF 05
  var_list = out_directory+'chr'+i+'_MAF_05.in'
  variants = [line.strip() for line in open(var_list).readlines()]
  input_f = vcf_directory + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + i + '_vcf_MAF_05.txt'
  output_f = out_directory + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + i + '_vcf_MAF_05_ld_pruned.txt'
  f = open(input_f)
  out_f = open(output_f, 'w')
  header = f.readline()
  out_f.write(header)
  for line in f.readlines():
    ID = line.split('\t')[2]
    if ID in variants:
      out_f.write(line)
  out_f.close()
