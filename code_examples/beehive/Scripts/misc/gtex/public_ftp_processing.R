##################################################
#  public_ftp_processing.R
#
#  $proj/Scripts/misc/gtex/public_ftp_processing.R
# 
#  
#
#  Author: Brian Jo
#
##################################################

# prepare gene set
expr_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/'
suffix = '_v6p_consortium_autosomes_normalized.txt'

list_files = list.files(path=expr_dir, pattern=suffix)

for (file in list_files) {
	tissue_name = strsplit(file, '_')[[1]][1]
	print(tissue_name)
	expression_matrix = read.csv(file = paste0(expr_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	rownames(expression_matrix) = expression_matrix$gene_id
	expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]
	mappability_cutoff = 0.8
	mappability_list = read.table(paste0('/tigress/BEE/RNAseq/Data/Resources/annotations/avg_mappability_Exon_UTR.txt'), stringsAsFactors=F)
	rownames(mappability_list) = mappability_list$V1
	# Arbitrary 0.8 cutoff - can be modified for a different threshold
	mappability_list = mappability_list[(mappability_list$V2>mappability_cutoff),]
	# Filter out genes with low mappability
	expression_matrix = expression_matrix[rownames(expression_matrix) %in% rownames(mappability_list),]

	gene_list = data.frame(genes = rownames(expression_matrix))
	write.table(gene_list, file = paste0('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/gene_set_ldacc/', tissue_name, '_gene_list.txt'), quote=F, row.names=F, col.names=F)
}


# convert .RData to .txt files
in_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/eqtl_list_pval_thresh_1e-5/'
out_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/main_consortium_trans/eqtl_list_pval_1e-5/'

list_files = list.files(path=in_dir)
suffix = '_trans_list_by_pvalue.txt'

for (file in list_files) {
	load(paste0(in_dir, file))
	tissue_name = strsplit(file, '_')[[1]][1]
	print(tissue_name)
	write.table(cumulative_me_trans_total, file = paste0(out_dir, tissue_name, suffix), quote=F, row.names=F, sep='\t')
	#write.table(cumulative_me_trans_total, file = gzfile(paste0(out_dir, tissue_name, suffix, '.gz'), 'w'), quote=F, row.names=F, sep='\t')
}


in_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/intrachromosomal/eqtl_list/'
out_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/main_consortium_intrachromosomal/eqtl_list_pval_1e-5/'

list_files = list.files(path=in_dir)
suffix = '_trans_list_by_pvalue.txt'

for (file in list_files) {
	load(paste0(in_dir, file))
	tissue_name = strsplit(file, '_')[[1]][1]
	print(tissue_name)
	write.table(cumulative_me_trans_total, file = paste0(out_dir, tissue_name, suffix), quote=F, row.names=F, sep='\t')
	# write.table(cumulative_me_trans_total, file = gzfile(paste0(out_dir, tissue_name, suffix, '.gz'), 'w'), quote=F, row.names=F, sep='\t')
}

in_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/intrachromosomal/eqtl_list/'
out_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/main_consortium_intrachromosomal/eqtl_list_pval_1e-5/'

list_files = list.files(path=in_dir)
suffix = '_trans_list_by_pvalue.txt'

for (file in list_files) {
	load(paste0(in_dir, file))
	tissue_name = strsplit(file, '_')[[1]][1]
	print(tissue_name)
	write.table(cumulative_me_trans_total, file = paste0(out_dir, tissue_name, suffix), quote=F, row.names=F, sep='\t')
	# write.table(cumulative_me_trans_total, file = gzfile(paste0(out_dir, tissue_name, suffix, '.gz'), 'w'), quote=F, row.names=F, sep='\t')
}

in_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/subset_runs/eqtl_list/'
out_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/main_consortium_subsetting_runs/eqtl_list_pval_1e-5/'

for (subsetting in c('cis_subset', 'gwas_subset', 'ld_subset')) {
	list_files = list.files(path=in_dir, pattern=subsetting)
	suffix = '_trans_list_by_pvalue.txt'

	for (file in list_files) {
		load(paste0(in_dir, file))
		tissue_name = strsplit(file, '_')[[1]][1]
		print(tissue_name)

		if (subsetting == 'cis_subset') {
			cumulative_me_trans_total = cis_subset_cumulative_me_trans_total
		} else if (subsetting == 'gwas_subset') {
			cumulative_me_trans_total = gwas_subset_cumulative_me_trans_total
		} else if (subsetting == 'ld_subset') {
			cumulative_me_trans_total = ld_subset_cumulative_me_trans_total
		}

		write.table(cumulative_me_trans_total, file = paste0(out_dir, tissue_name, '_', subsetting, suffix), quote=F, row.names=F, sep='\t')
		# write.table(cumulative_me_trans_total, file = gzfile(paste0(out_dir, tissue_name, suffix, '.gz'), 'w'), quote=F, row.names=F, sep='\t')
	}
}
