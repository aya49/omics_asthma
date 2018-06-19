##################################################
#  quantify_kallisto_silver.py
#
#  $proj/Scripts/processing/silver/quantify_kallisto_silver.py
#
#  Create expression matrix for kallisto
#
#  Authors: Brian Jo
#
##################################################
import pandas as pd
import pickle
from collections import OrderedDict
from sys import argv

tissue = argv[1]
# tissue = 'ovary'

tables_dir = '/tigress/BEE/RNAseq/Data/Resources/gtex/tables/'
quant_dir = '/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/kallisto_hg38/'
out_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/silver/kallisto/'

sample_table = pd.read_csv(tables_dir + 'sample_table.txt', sep='\t')
subject_table = pd.read_csv(tables_dir + 'subject_table.txt', sep='\t')

# which subjects/samples have genotypes?
subjects_with_geno = list(subject_table[subject_table['genotype_avail'] == True]['submitted_subject_id_s'])
subjects_without_geno = list(subject_table[subject_table['genotype_avail'] == False]['submitted_subject_id_s'])

samples_with_geno = sample_table[sample_table['submitted_subject_id_s'].isin(subjects_with_geno)]
samples_without_geno = sample_table[sample_table['submitted_subject_id_s'].isin(subjects_without_geno)]
isoform_suffix = 'abundance.tsv'

# create table for each tissue
print(tissue)
sample_list_with_geno = list(samples_with_geno[samples_with_geno['tissue_name'] == tissue]['Run_s'])
sample_list_without_geno = list(samples_without_geno[samples_without_geno['tissue_name'] == tissue]['Run_s'])
subject_list_with_geno = list(samples_with_geno[samples_with_geno['tissue_name'] == tissue]['submitted_subject_id_s'])
subject_list_without_geno = list(samples_without_geno[samples_without_geno['tissue_name'] == tissue]['submitted_subject_id_s'])

# Determine if all tissue samples have genotypes
if len(sample_list_without_geno) > 0:
    sample_list = [x for sublist in [sample_list_with_geno, sample_list_without_geno] for x in sublist]
    subject_list = [x for sublist in [subject_list_with_geno, subject_list_without_geno] for x in sublist]
else:
    sample_list = sample_list_with_geno
    subject_list = subject_list_with_geno

# Accumulate quantification across all samples in a tissue
isoform_expr_mat_TPM = pd.DataFrame()
isoform_expr_mat_counts = pd.DataFrame()
for i in range(len(sample_list)):
    print(str(i))
    sample = sample_list[i]
    subject = subject_list[i]
    # gene_table = pd.read_csv(quant_dir + sample + '/' + sample + gene_suffix, sep='\t')
    isoform_table = pd.read_csv(quant_dir + sample + '/' + isoform_suffix, sep='\t')
    
    isoform_expr_TPM = isoform_table[['tpm']]
    isoform_expr_counts = isoform_table[['est_counts']]
    # name expression column by subject
    isoform_expr_TPM.columns = [subject]
    isoform_expr_counts.columns = [subject]
    # build entire expression matrix incrementally by column (subject)
    isoform_expr_mat_TPM = pd.concat([isoform_expr_mat_TPM, isoform_expr_TPM], axis=1)
    isoform_expr_mat_counts = pd.concat([isoform_expr_mat_counts, isoform_expr_counts], axis=1)

# set row names
isoform_expr_mat_TPM.index = isoform_table['target_id'].tolist()
isoform_expr_mat_counts.index = isoform_table['target_id'].tolist()

# load transcript annotation
transcript_anno_file = '/tigress/BEE/RNAseq/Data/Resources/annotations/silver/gencode.v26.transcripts'
transcript_dict = pickle.load(open(transcript_anno_file, 'rb'))

gene_expr_mat_TPM_filter = pd.DataFrame()
isoform_expr_mat_TPM_filter = pd.DataFrame()
gene_expr_mat_counts_filter = pd.DataFrame()
isoform_expr_mat_counts_filter = pd.DataFrame()

temp_gene_expr_mat_TPM = pd.DataFrame()
temp_isoform_expr_mat_TPM = pd.DataFrame()
temp_gene_expr_mat_counts = pd.DataFrame()
temp_isoform_expr_mat_counts = pd.DataFrame()
all_chrs = [str(x) for x in range(1,23)] + ['X','Y','M']

# With ordered dict, transcripts come in gencode order (but it's not in order for kallisto files)
cur_gene = ''
for cur_chr in all_chrs:
    print(cur_chr)
    for k,v in transcript_dict[cur_chr].items():
        # is the transcript in the isoform table?
        if k not in isoform_expr_mat_TPM.index:
            continue
        # new gene
        if v['gene_id'] != cur_gene:
            if cur_gene != '':
                # sum across columns for gene counts and TPM
                temp_gene_expr_mat_TPM = temp_isoform_expr_mat_TPM.sum(axis = 0).to_frame()
                temp_gene_expr_mat_counts = temp_isoform_expr_mat_counts.sum(axis = 0).to_frame()
                temp_gene_expr_mat_TPM.columns = [v['gene_id']]
                temp_gene_expr_mat_counts.columns = [v['gene_id']]
                n_zeros = (temp_gene_expr_mat_TPM == 0).sum()
                # If there are at least 10 nonzero samples, write all transcripts in that gene
                if (n_zeros <= (isoform_expr_mat_TPM.shape[1] - 10))[0]:
                    temp_isoform_expr_mat_TPM['gene'] = cur_gene
                    temp_isoform_expr_mat_counts['gene'] = cur_gene
                    # add on to the current table
                    isoform_expr_mat_TPM_filter = isoform_expr_mat_TPM_filter.append(temp_isoform_expr_mat_TPM)
                    isoform_expr_mat_counts_filter = isoform_expr_mat_counts_filter.append(temp_isoform_expr_mat_counts)
                    gene_expr_mat_TPM_filter = gene_expr_mat_TPM_filter.append(temp_gene_expr_mat_TPM.transpose())
                    gene_expr_mat_counts_filter = gene_expr_mat_counts_filter.append(temp_gene_expr_mat_counts.transpose())
            cur_gene = v['gene_id']
            temp_gene_expr_mat_TPM = pd.DataFrame()
            temp_gene_expr_mat_counts = pd.DataFrame()
            temp_isoform_expr_mat_TPM = pd.DataFrame()
            temp_isoform_expr_mat_counts = pd.DataFrame()
        temp_isoform_expr_mat_TPM = temp_isoform_expr_mat_TPM.append(isoform_expr_mat_TPM.loc[[k]])
        temp_isoform_expr_mat_counts = temp_isoform_expr_mat_counts.append(isoform_expr_mat_counts.loc[[k]])

# Add the last row
temp_gene_expr_mat_TPM = temp_isoform_expr_mat_TPM.sum(axis = 0).to_frame()
temp_gene_expr_mat_counts = temp_isoform_expr_mat_counts.sum(axis = 0).to_frame()
temp_gene_expr_mat_TPM.columns = [v['gene_id']]
temp_gene_expr_mat_counts.columns = [v['gene_id']]
n_zeros = (temp_gene_expr_mat_TPM == 0).sum()
if (n_zeros <= (isoform_expr_mat_TPM.shape[1] - 10))[0]:
    temp_isoform_expr_mat_TPM['gene'] = cur_gene
    temp_isoform_expr_mat_counts['gene'] = cur_gene
    # add on to the current table
    isoform_expr_mat_TPM_filter = isoform_expr_mat_TPM_filter.append(temp_isoform_expr_mat_TPM)
    isoform_expr_mat_counts_filter = isoform_expr_mat_counts_filter.append(temp_isoform_expr_mat_counts)
    gene_expr_mat_TPM_filter = gene_expr_mat_TPM_filter.append(temp_gene_expr_mat_TPM.transpose())
    gene_expr_mat_counts_filter = gene_expr_mat_counts_filter.append(temp_gene_expr_mat_counts.transpose())

gene_expr_mat_TPM_filter[subject_list_with_geno].to_csv(out_dir + 'genes/' + tissue + '_kallisto_TPM_with_genotype.txt', sep='\t', float_format='%.4f')
gene_expr_mat_counts_filter[subject_list_with_geno].to_csv(out_dir + 'genes/' + tissue + '_kallisto_counts_with_genotype.txt', sep='\t', float_format='%.4f')
isoform_expr_mat_TPM_filter[subject_list_with_geno].to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_TPM_with_genotype.txt', sep='\t', float_format='%.4f')
isoform_expr_mat_counts_filter[subject_list_with_geno].to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_counts_with_genotype.txt', sep='\t', float_format='%.4f')

if len(sample_list_without_geno) > 0:
    gene_expr_mat_TPM_filter[subject_list_without_geno].to_csv(out_dir + 'genes/' + tissue + '_kallisto_TPM_without_genotype.txt', sep='\t', float_format='%.4f')
    gene_expr_mat_counts_filter[subject_list_without_geno].to_csv(out_dir + 'genes/' + tissue + '_kallisto_counts_without_genotype.txt', sep='\t', float_format='%.4f')
    isoform_expr_mat_TPM_filter[subject_list_without_geno].to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_TPM_without_genotype.txt', sep='\t', float_format='%.4f')
    isoform_expr_mat_counts_filter[subject_list_without_geno].to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_counts_without_genotype.txt', sep='\t', float_format='%.4f')
