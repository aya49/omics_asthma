##################################################
#  prepare_gwas_subset.R
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/prepare_gwas_subset.R
# 
#  This script saves a separate genotype file for the gwas-subsetting run, by taking the best cis SNP, in a tissue-specific fashion
#
#  Author: Brian Jo
#
##################################################

args <-commandArgs(TRUE)
# tissue = 'uterus'
tissue = args[1]

gwas_df = read.csv('/tigress/BEE/eQTLs/Data/References/Regulatory/gwas_catalog_v1.0.1-associations_e84_r2016-06-12.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)
#dim(gwas_df)
#25226    36
SNP_list = unlist(sapply(gwas_df$SNPS, function(x) {strsplit(x, ', ')[[1]]}))
SNP_list = unique(SNP_list)

SNP_info_dfs = vector('list', 22)
#SNP_info_dfs_rand = vector('list', 22)
SNP_info_dfs_rand = vector('list', 242)
SNP_info_dfs_rand_no_tss = vector('list', 242)
genotype_mats = vector('list', 22)
#genotype_mats_rand = vector('list', 22)
genotype_mats_rand = vector('list', 242)
genotype_mats_rand_no_tss = vector('list', 242)


# Import tissue-specific subject list:
expr_dir = '/tigress/BEE/eQTLs/Data/Expression/Expression_matrices/GTEx/hg19/filtered/normalized_gtex_gcp_v6p/'
file_suffix = '_nonoverlapping_certain_biotype_normalized_GTEx_gct_RNA-SeQCv1.1.8_gene_rpkm.txt'
subject_list = colnames(read.table(paste(expr_dir, tissue, file_suffix, sep=''), nrows = 1, header =TRUE, sep ='\t'))

# Also import locations of exons:
annotation_gtf = read.table('/tigress/BEE/eQTLs/Data/References/Annotations/gencode.v19.annotation.gtf', skip = 5, header = FALSE, sep ='\t')
annotation_gtf = annotation_gtf[annotation_gtf$V3 == 'exon', c('V1','V2','V3','V4','V5')]

gene_pos = read.table('/tigress/BEE/eQTLs/Data/References/GTEx/gene_certain_biotype_list.bed', header=FALSE, stringsAsFactors=FALSE, sep='\t')
colnames(gene_pos)=c('chr','start','end','gene','strand')

for (i in c(1:22)) {
	print(i)
	annotation_gtf_part = annotation_gtf[annotation_gtf$V1 == paste('chr',i,sep=''),]
	metadata_file = paste('/tigress/BEE/eQTLs/Data/Genotype/GTEx/SNP_metadata/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr',i,'.metadata.expanded.txt',sep='')
	snp_metadata = read.table(metadata_file, header=TRUE, stringsAsFactors=FALSE, sep='\t')
	snp_metadata = snp_metadata[snp_metadata$maf05_FILTER == 1,]
	snp_metadata = snp_metadata[!(snp_metadata$IN_REPEAT),]
	snp_metadata = snp_metadata[!(duplicated(snp_metadata$dbSNPID)),]

	# Also read in the genotypes:
	genotype_in_file = paste('/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr',i,'_Final.txt',sep='')
	genotype_matrix = read.table(genotype_in_file, header=TRUE, stringsAsFactors=FALSE, sep='\t')

	# First process the GWAS snps - which variants to we have information for?
	SNP_info_df_part = snp_metadata[snp_metadata$dbSNPID %in% SNP_list,c('CHROM', 'POS', 'dbSNPID', 'maf05_FILTER', 'A1_FREQ', 'IN_REPEAT', 'CIS.GENE', 'IN.GENE')]
	rownames(SNP_info_df_part) = SNP_info_df_part$dbSNPID
	SNP_info_df_part = SNP_info_df_part[SNP_info_df_part$dbSNPID %in% rownames(genotype_matrix),]

	# Tissue_specific SNP MAF
	genotype_mat_part = genotype_matrix[SNP_info_df_part$dbSNPID,subject_list]
	MAF = rowSums(genotype_mat_part)/dim(genotype_mat_part)[2]/2
	genotype_mat_part = genotype_mat_part[MAF>0.05 & MAF<0.95,]
	SNP_info_df_part = SNP_info_df_part[MAF>0.05 & MAF<0.95,]

	gene_part = gene_pos[gene_pos$chr == i,]
	SNP_info_df_part$tss_dist = sapply(SNP_info_df_part$POS, function(x) { min(min(abs(gene_part$start[gene_part$strand == '+']-x)),min(abs(gene_part$end[gene_part$strand == '-']-x))) })
	SNP_info_df_part$MAF = MAF[MAF>0.05 & MAF<0.95]
	SNP_info_df_part$MAF = sapply(SNP_info_df_part$MAF, function(x) {min(x, 1-x)})

	SNP_info_dfs[[i]] = SNP_info_df_part
	genotype_mats[[i]] = genotype_mat_part

	# Now process the random selection:
	# We have two sets - one with tss_dist constraint, one without. We also import exon positions to make sure the SNPs are non-coding
	for (n in c(1:11)) {
		ind = (i-1)*11+n
		print(ind)
		# First select a subset from snp_metadata to match gwas tss distance distribution
		SNP_info_dfs_rand_part = snp_metadata[,c('CHROM', 'POS', 'dbSNPID', 'maf05_FILTER', 'A1_FREQ', 'IN_REPEAT', 'CIS.GENE', 'IN.GENE')]
		# Just in case: remove all snps that overlap
		SNP_info_dfs_rand_part = SNP_info_dfs_rand_part[!(SNP_info_dfs_rand_part$dbSNPID %in% SNP_info_df_part),]

		# Sample extra for the constraints - let's say 20X
		set.seed(ind)
		r_inds = sample(c(1:dim(SNP_info_dfs_rand_part)[1]), dim(SNP_info_df_part)[1]*20, replace = FALSE)
		SNP_sample_df = SNP_info_dfs_rand_part[r_inds,]
		rownames(SNP_sample_df) = SNP_sample_df$dbSNPID
		SNP_sample_df = SNP_sample_df[SNP_sample_df$dbSNPID %in% rownames(genotype_matrix),]

		# MAF filter
		genotype_mat_rand_part = genotype_matrix[SNP_sample_df$dbSNPID,subject_list]
		MAF = rowSums(genotype_mat_rand_part)/dim(genotype_mat_rand_part)[2]/2
		SNP_sample_df = SNP_sample_df[MAF>0.05 & MAF<0.95,]

		SNP_sample_df$tss_dist = sapply(SNP_sample_df$POS, function(x) { min(min(abs(gene_part$start[gene_part$strand == '+']-x)),min(abs(gene_part$end[gene_part$strand == '-']-x))) })
		SNP_sample_df$MAF = MAF[MAF>0.05 & MAF<0.95]
		SNP_sample_df$MAF = sapply(SNP_sample_df$MAF, function(x) {min(x, 1-x)})

		# Now we select a subset from snp_metadata to match MAF distribution
		temp = sapply(SNP_sample_df$MAF, function(x) {sum(x >= c(1:9)/20)})
		temp2 = sapply(SNP_info_df_part$MAF, function(x) {sum(x >= c(1:9)/20)})

		# Make sure that we always have enough SNPs so that the code doesn't error out
		multiplier = min(5, min(floor(sapply(c(1:9), function(x) { sum(temp==x)/sum(temp2==x) }))))
		rr_inds = c()
		for (k in c(1:9)) {
			# sample extra for the second-pass MAF filter
			set.seed(100*k + 10*n + i)
			rr_inds = c(rr_inds, sample(which(temp==k), sum(temp2==k)*multiplier, replace = FALSE))
		}

		# Make sure the snps are non-coding:
		SNP_sample_df$in_exon = sapply(SNP_sample_df$POS, function(x) { max(colSums(rbind(x >= annotation_gtf_part$V4, annotation_gtf_part$V5 >= x))) })
		SNP_sample_df$in_exon = SNP_sample_df$in_exon-1
		SNP_sample_df = SNP_sample_df[SNP_sample_df$in_exon == 0,]

		# Selection without tss_dist
		set.seed(10*n + i)
		r_inds = sample(c(1:dim(SNP_sample_df)[1]), dim(SNP_info_df_part)[1], replace = FALSE)
		SNP_info_dfs_rand_no_tss[[ind]] = SNP_sample_df[r_inds,]
		genotype_mats_rand_no_tss[[ind]] = genotype_matrix[SNP_info_dfs_rand_no_tss[[ind]]$dbSNPID,subject_list]

		# Selection with MAF filter

		temp = sapply(log10(SNP_sample_df$tss_dist), function(x) {sum(x > c(8:12)/2)})
		temp2 = sapply(log10(SNP_info_df_part$tss_dist), function(x) {sum(x > c(8:12)/2)})
		rr_inds = c()
		for (k in c(0:min(max(temp),max(temp2)))) {
			set.seed(1000*k + 10*n + i)
			rr_inds = c(rr_inds, sample(which(temp==k), sum(temp2==k), replace = FALSE))
		}
		SNP_info_dfs_rand[[ind]] = SNP_sample_df[rr_inds,]
		genotype_mats_rand[[ind]] = genotype_matrix[SNP_info_dfs_rand[[ind]]$dbSNPID,subject_list]
	}
}

genotype_mat_gwas = do.call("rbind", genotype_mats)
# this should be 0
print(sum(is.na(genotype_mat_gwas[,1])))

output_dir1 = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous/GWAS_snp/'
write.table(genotype_mat_gwas, file=paste(output_dir1, 'GTEx_genotypes_maf05_continuous_not_in_repeat_GWAS_', tissue, '.txt', sep=''), quote = FALSE, sep = "\t")

for (i in c(1:11)) {
	genotype_mat_rand = do.call("rbind", genotype_mats_rand[c(0:21)*11+i])
	genotype_mat_rand_no_tss = do.call("rbind", genotype_mats_rand_no_tss[c(0:21)*11+i])

	# this should also be 0
	print(sum(is.na(genotype_mat_rand[,1])))
	print(sum(is.na(genotype_mat_rand_no_tss[,1])))

	output_dir2 = '/tigress/BEE/eQTLs/Data/Genotype/GTEx/imputed_genotypes/allelic_dosage/continuous_random/GWAS_snp/'
	write.table(genotype_mat_rand, file=paste(output_dir2, 'GTEx_genotypes_maf05_continuous_not_in_repeat_GWAS_', tissue, '_', i, '.txt', sep=''), quote = FALSE, sep = "\t")
	write.table(genotype_mat_rand_no_tss, file=paste(output_dir2, 'GTEx_genotypes_maf05_continuous_not_in_repeat_GWAS_', tissue, '_', i, '_no_tss_filter.txt', sep=''), quote = FALSE, sep = "\t")
}

# SNP_info = do.call("rbind", SNP_info_dfs)
# SNP_info_rand = do.call("rbind", SNP_info_dfs_rand_part)

# save(SNP_info, SNP_info_rand, file='/home/bj5/temp_artery.RData')
