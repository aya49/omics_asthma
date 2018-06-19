##################################################
#  quantify_RSEM_silver.py
#
#  $proj/Scripts/processing/silver/quantify_RSEM_silver.py
#
#  Create expression matrix for RSEM
#
#  Authors: Brian Jo
#
##################################################
import pandas as pd
from sys import argv

tables_dir = '/tigress/BEE/RNAseq/Data/Resources/gtex/tables/'
quant_dir = '/tigress/BEE/gtex/data/phenotype/expression/quantified_rna_seq_reads/silver/RSEM_hg38/'
out_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/silver/RSEM/'
gene_suffix = '.genes.results'
isoform_suffix = '.isoforms.results'

sample_table = pd.read_csv(tables_dir + 'sample_table.txt', sep='\t')
subject_table = pd.read_csv(tables_dir + 'subject_table.txt', sep='\t')
tissue_table = pd.read_csv(tables_dir + 'tissue_table.txt', sep='\t')

# which subjects/samples have genotypes?
subjects_with_geno = list(subject_table[subject_table['genotype_avail'] == True]['submitted_subject_id_s'])
subjects_without_geno = list(subject_table[subject_table['genotype_avail'] == False]['submitted_subject_id_s'])

samples_with_geno = sample_table[sample_table['submitted_subject_id_s'].isin(subjects_with_geno)]
samples_without_geno = sample_table[sample_table['submitted_subject_id_s'].isin(subjects_without_geno)]

# create table for each tissue
for tissue in tissue_table['tissue_name']:
    print(tissue)
    sample_list_with_geno = list(samples_with_geno[samples_with_geno['tissue_name'] == tissue]['Run_s'])
    sample_list_without_geno = list(samples_without_geno[samples_without_geno['tissue_name'] == tissue]['Run_s'])
    subject_list_with_geno = list(samples_with_geno[samples_with_geno['tissue_name'] == tissue]['submitted_subject_id_s'])
    subject_list_without_geno = list(samples_without_geno[samples_without_geno['tissue_name'] == tissue]['submitted_subject_id_s'])

    gene_expr_mat_TPM = pd.DataFrame()
    gene_expr_mat_counts = pd.DataFrame()
    isoform_expr_mat_TPM = pd.DataFrame()
    isoform_expr_mat_counts = pd.DataFrame()
    for i in range(len(sample_list_with_geno)):
        print(str(i))
        sample = sample_list_with_geno[i]
        subject = subject_list_with_geno[i]
        gene_table = pd.read_csv(quant_dir + sample + '/' + sample + gene_suffix, sep='\t')
        isoform_table = pd.read_csv(quant_dir + sample + '/' + sample + isoform_suffix, sep='\t')
        
        gene_expr_TPM = gene_table[['TPM']]
        gene_expr_counts = gene_table[['expected_count']]
        # name expression column by subject
        gene_expr_TPM.columns = [subject]
        gene_expr_counts.columns = [subject]
        # build entire expression matrix incrementally by column (subject)
        gene_expr_mat_TPM = pd.concat([gene_expr_mat_TPM, gene_expr_TPM], axis=1)
        gene_expr_mat_counts = pd.concat([gene_expr_mat_counts, gene_expr_counts], axis=1)

        isoform_expr_TPM = isoform_table[['TPM']]
        isoform_expr_counts = isoform_table[['expected_count']]
        # name expression column by subject
        isoform_expr_TPM.columns = [subject]
        isoform_expr_counts.columns = [subject]
        # build entire expression matrix incrementally by column (subject)
        isoform_expr_mat_TPM = pd.concat([isoform_expr_mat_TPM, isoform_expr_TPM], axis=1)
        isoform_expr_mat_counts = pd.concat([isoform_expr_mat_counts, isoform_expr_counts], axis=1)

        # program may use a string like '.' to signify 0, in which case pandas will turn it into NA
        # expr_mat_TPM_uncertain = expr_mat_TPM_uncertain.dropna(axis=0)
        # expr_mat_TPM_certain = expr_mat_TPM_certain.dropna(axis=0)
        # expr_mat_counts_uncertain = expr_mat_counts_uncertain.dropna(axis=0)
        # expr_mat_counts_certain = expr_mat_counts_certain.dropna(axis=0)

    # write files
    gene_expr_mat_TPM.to_csv(out_dir + 'genes/' + tissue + '_RSEM_TPM_with_genotype.txt', sep='\t')
    gene_expr_mat_counts.to_csv(out_dir + 'genes/' + tissue + '_RSEM_counts_with_genotype.txt', sep='\t')
    isoform_expr_mat_TPM.to_csv(out_dir + 'isoforms/' + tissue + '_RSEM_TPM_with_genotype.txt', sep='\t')
    isoform_expr_mat_counts.to_csv(out_dir + 'isoforms/' + tissue + '_RSEM_counts_with_genotype.txt', sep='\t')

    # If there are samples without genotype:
    gene_expr_mat_TPM = pd.DataFrame()
    gene_expr_mat_counts = pd.DataFrame()
    isoform_expr_mat_TPM = pd.DataFrame()
    isoform_expr_mat_counts = pd.DataFrame()
    for i in range(len(sample_list_without_geno)):
        print(str(i))
        sample = sample_list_without_geno[i]
        subject = subject_list_without_geno[i]
        gene_table = pd.read_csv(quant_dir + sample + '/' + sample + gene_suffix, sep='\t')
        isoform_table = pd.read_csv(quant_dir + sample + '/' + sample + isoform_suffix, sep='\t')
        
        gene_expr_TPM = gene_table[['TPM']]
        gene_expr_counts = gene_table[['expected_count']]
        # name expression column by subject
        gene_expr_TPM.columns = [subject]
        gene_expr_counts.columns = [subject]
        # build entire expression matrix incrementally by column (subject)
        gene_expr_mat_TPM = pd.concat([gene_expr_mat_TPM, gene_expr_TPM], axis=1)
        gene_expr_mat_counts = pd.concat([gene_expr_mat_counts, gene_expr_counts], axis=1)

        isoform_expr_TPM = isoform_table[['TPM']]
        isoform_expr_counts = isoform_table[['expected_count']]
        # name expression column by subject
        isoform_expr_TPM.columns = [subject]
        isoform_expr_counts.columns = [subject]
        # build entire expression matrix incrementally by column (subject)
        isoform_expr_mat_TPM = pd.concat([isoform_expr_mat_TPM, isoform_expr_TPM], axis=1)
        isoform_expr_mat_counts = pd.concat([isoform_expr_mat_counts, isoform_expr_counts], axis=1)

        # program may use a string like '.' to signify 0, in which case pandas will turn it into NA
        # expr_mat_TPM_uncertain = expr_mat_TPM_uncertain.dropna(axis=0)
        # expr_mat_TPM_certain = expr_mat_TPM_certain.dropna(axis=0)
        # expr_mat_counts_uncertain = expr_mat_counts_uncertain.dropna(axis=0)
        # expr_mat_counts_certain = expr_mat_counts_certain.dropna(axis=0)

    # write files
    gene_expr_mat_TPM.to_csv(out_dir + 'genes/' + tissue + '_RSEM_TPM_without_genotype.txt', sep='\t')
    gene_expr_mat_counts.to_csv(out_dir + 'genes/' + tissue + '_RSEM_counts_without_genotype.txt', sep='\t')
    isoform_expr_mat_TPM.to_csv(out_dir + 'isoforms/' + tissue + '_RSEM_TPM_without_genotype.txt', sep='\t')
    isoform_expr_mat_counts.to_csv(out_dir + 'isoforms/' + tissue + '_RSEM_counts_without_genotype.txt', sep='\t')

