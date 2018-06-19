##################################################
#  trans_matrix_eqtl_subset_runs_final.R
#
#  $proj/Scripts/eqtls/trans/gtex/consortium/trans_matrix_eqtl_subset_runs_final.R
# 
#  This version performs the subsetting runs (cis-eQTLs, GWAS variants, and LD) for the consortium manuscript.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# Example
# args = c(1:7)
# args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/wholeblood_v6p_consortium_autosomes_normalized.txt'
# args[2] = '8'
# args[3] = '2.5e8'
# args[4] = 'wholeblood'
# args[5] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg19/GTEx_phs000424_v6p/normalized/covariates/'
# args[6] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/subset_runs/wholeblood/wholeblood_nonverlapping_certain_autosomes_normalized_MatrixEQTL'
# # Attach PEER number at the end
# args[7] = '1'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
geno_option = 'continuous'
part_number = as.numeric(args[2])
cis_dist = as.numeric(args[3])
tissue_name = args[4]
cov_dir = args[5]
out_file = args[6]
num_split = as.numeric(args[7])
chr_number = ceiling(part_number/num_split)

# Take all SNPs for the trans runs
# genotype_file_name = paste(args[8], 'Genotype/GTEx/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_', disc_or_cont , '_Chr', chr_number, '_Final.txt', sep="")
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/', geno_option, '/GTEx_genotypes_maf05_continuous_not_in_repeat_Chr', chr_number, '_Final.txt', sep="")
SNP_position_file = paste0(proj_dir, '/Data/Genotype/gtex/SNP_positions_hg19/SNP_positions_Chr', chr_number, '.txt', sep="")
gene_position_file = paste0(proj_dir, '/Data/Expression/gene_metadata_hg19/gene_metadata_chrAll.txt')

# Need to check that the subjects read in from expression matrices match the subjects in genotype data
library(dplyr)

expression_matrix = read.csv(file = expression_file_location, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,2:dim(expression_matrix)[2]]

# Load in the mappability list:
mappability_cutoff = 0.8
mappability_list = read.table(paste0(proj_dir, '/Data/Resources/annotations/avg_mappability_Exon_UTR.txt'), stringsAsFactors=F)
rownames(mappability_list) = mappability_list$V1
# Arbitrary 0.8 cutoff - can be modified for a different threshold
mappability_list = mappability_list[(mappability_list$V2>mappability_cutoff),]
# Filter out genes with low mappability
expression_matrix = expression_matrix[rownames(expression_matrix) %in% rownames(mappability_list),]

genotype_matrix = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

# Only for all-by-all runs - divide the chromosomes into ten because the jobs because there are a lot of tests
snp_set_size = nrow(genotype_matrix)

# Also obtain genes and snps that are in the expression and genotype files, respectively.
gene_positions = read.table(file = gene_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gene_positions = gene_positions[,c('gene_id', 'chr', 'start', 'end')]
gene_positions = mutate(gene_positions, chr = paste("chr", chr, sep = ""))
gene_positions = select(gene_positions, c(gene_id, chr, start, end))

snp_positions = read.table(SNP_position_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(snp_positions) = snp_positions$rsID
snp_positions = snp_positions[rownames(snp_positions) %in% rownames(genotype_matrix),]
genotype_matrix = genotype_matrix[rownames(genotype_matrix) %in% rownames(snp_positions),]

# Make sure the list of genes and their ordering are the same
# We are also filtering out non-autosomal genes
expression_matrix = expression_matrix[rownames(expression_matrix) %in% gene_positions$gene_id,]
gene_positions = gene_positions[gene_positions$gene_id %in% rownames(expression_matrix),]
rownames(gene_positions) = gene_positions$gene_id
gene_positions = gene_positions[rownames(expression_matrix),]
# gene_positions = gene_positions[sapply(rownames(expression_matrix), function(x) {match(x, gene_positions$gene_id)}),]
# rownames(gene_positions) = gene_positions$gene_id

# Incorporate a tissue-specific MAF filter so that we are not testing variants that actually have a low MAF in our test set
MAF = rowSums(genotype_matrix)/dim(genotype_matrix)[2]/2
MAF = sapply(MAF, function(x) {min(x, 1-x)})
genotype_matrix = genotype_matrix[MAF>0.05,]
snp_positions = snp_positions[MAF>0.05,]
MAF = MAF[MAF>0.05]

# Load in the covariates
suffix = '_Analysis.covariates.txt'
covars = read.csv(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

print('finished loading data files')

# Load the SNP metadata file
metadata_file = paste('/tigress/BEE/RNAseq/Data/Genotype/gtex/SNP_metadata/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_allchr_genot_imput.chr',chr_number,'.metadata.txt',sep='')
snp_metadata = read.table(metadata_file, header=TRUE, stringsAsFactors=FALSE, sep='\t')
snp_metadata = snp_metadata[!duplicated(snp_metadata$dbSNPID),]
rownames(snp_metadata) = snp_metadata$dbSNPID
snp_metadata = snp_metadata[snp_metadata$dbSNPID %in% rownames(genotype_matrix),]
snp_metadata = snp_metadata[rownames(genotype_matrix),]
# record MAF
snp_metadata$MAF = MAF[snp_metadata$dbSNPID]
# record TSS
gene_pos_subset = gene_positions[gene_positions$chr == paste0('chr', snp_metadata$CHROM[1]),]
snp_metadata$TSS = sapply(snp_metadata$POS, function(x) {min(abs(gene_pos_subset$start - x))})

print('finished loading SNP metadata')

# Import MatrixEQTL
library(MatrixEQTL)
source(paste0(proj_dir, '/Scripts/eqtls/MatrixEQTL_wrapper.R'))

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1e-5
gene_set_size = dim(expression_matrix)[1]

# Prepare subsetting runs
# First prepare the cis subsetting:
cis_subset_dir = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/subset_SNPs/cis-best-subset/')

if (file.exists(paste0(cis_subset_dir, tissue_name, '_part', part_number, '.txt'))) {
	cis_best_snp = read.table(paste0(cis_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', header=T, stringsAsFactors=F)
} else {
	fastqtl_dir = '/tigress/BEE/RNAseq/Data/Resources/gtex/dbGaP/GTEx_phs000424/v6p_fastQTL_FOR_QC_ONLY/'
	suffix = '_Analysis.v6p.FOR_QC_ONLY.snpgenes.txt.gz'
	cis_eqtl_files = list.files(path = fastqtl_dir, pattern = suffix)
	ind = which(as.character(sapply(cis_eqtl_files, function(x) {tolower(gsub(" ", "", gsub("[[:punct:]]", "", strsplit(x, suffix)[[1]][1])))})) == tissue_name)

	snpgene_df = read.table(gzfile(paste(fastqtl_dir, cis_eqtl_files[ind], sep='')), header=TRUE, stringsAsFactors=FALSE, sep='\t')
	snpgene_df = snpgene_df[(snpgene_df$variant_id %in% snp_metadata$ID),]

	unique_genes = unique(snpgene_df$gene_id)
	cis_best_snp = vector('list', length(unique_genes))

	for (i in c(1:length(unique_genes))) {
		gene = unique_genes[i]
		snpgene_subset = snpgene_df[(which(snpgene_df$gene_id == gene)),]
		cis_best_snp[[i]] = snpgene_subset[snpgene_subset$pval_nominal == min(snpgene_subset$pval_nominal),]
	}

	cis_best_snp = do.call("rbind", cis_best_snp)
	rownames(snp_metadata) = snp_metadata$ID
	cis_best_snp$dbSNPID = snp_metadata[cis_best_snp$variant_id,]$dbSNPID
	cis_best_snp = cis_best_snp[(!duplicated(cis_best_snp$dbSNPID)),]
	rownames(snp_metadata) = snp_metadata$dbSNPID
	
	cis_best_snp$MAF = MAF[cis_best_snp$dbSNPID]
	cis_best_snp$position = as.numeric(sapply(cis_best_snp$variant_id, function(x) {strsplit(x, '_')[[1]][2]}))

	write.table(cis_best_snp, file = paste0(cis_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', quote=F, row.names=F)
}

print('cis-subsetting tests')

genotype_subset = genotype_matrix[cis_best_snp$dbSNPID,]
snp_position_subset = snp_positions[cis_best_snp$dbSNPID,]
cis_subset_size = nrow(cis_best_snp)

me_cis = MatrixEQTL_wrapper(genotype_subset, expression_matrix, snp_position_subset, gene_positions, pvThresh = pvOutputThreshold, pvThresh_cis = pvOutputThreshold, cis_dist = cis_dist, covariates = covars)

# get the random-matched subset for cis-variants:
cis_subset_dir = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/subset_SNPs/cis-best-subset-random/')

if (file.exists(paste0(cis_subset_dir, tissue_name, '_part', part_number, '.txt'))) {
	cis_best_snp_random = read.table(paste0(cis_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', header=T, stringsAsFactors=F)
} else {
	# selection scheme - for each test SNP, select random SNP that is matched in TSS and MAF, while not being too close to any of the other SNPs
	cis_best_snp_random = vector('list', nrow(cis_best_snp))
	for (i in c(1:nrow(cis_best_snp))) {
		# try to match - let's say 1000 bp and 0.01 in MAF
		tss_match = which(abs(snp_metadata$TSS - abs(cis_best_snp$tss_distance[i])) < 1000)
		maf_match = which(abs(snp_metadata$MAF - cis_best_snp$MAF[i]) < 0.01)
		both_match = intersect(tss_match, maf_match)
		# if too stringent, relax a little bit
		if (length(both_match) == 0) {
			tss_match = which(abs(snp_metadata$TSS - abs(cis_best_snp$tss_distance[i])) < 5000)
			maf_match = which(abs(snp_metadata$MAF - cis_best_snp$MAF[i]) < 0.03)
			both_match = intersect(tss_match, maf_match)
		}
		# if still too stringent, skip
		if (length(both_match) == 0) {next}
		# not in LD with any of the snps being tested - let's say 10000 bp
		snp_metadata_subset = snp_metadata[both_match,]
		snp_metadata_subset = snp_metadata_subset[sapply(snp_metadata_subset$POS, function(x) {sum(abs(cis_best_snp$position - x) < 10000) == 0}),]
		# if still too stringent, skip
		if (nrow(snp_metadata_subset) == 0) {next}
		# now randomly select
		set.seed = 42
		cis_best_snp_random[[i]] = snp_metadata_subset[sample(nrow(snp_metadata_subset), 1),]
	}

	cis_best_snp_random = do.call("rbind", cis_best_snp_random[!duplicated(sapply(cis_best_snp_random, function(x) {x$dbSNPID}))])
	write.table(cis_best_snp_random, file = paste0(cis_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', quote=F, row.names=F)
}

print('cis-subsetting tests, matched random')

genotype_subset = genotype_matrix[cis_best_snp_random$dbSNPID,]
snp_position_subset = snp_positions[cis_best_snp_random$dbSNPID,]
cis_subset_size_random = nrow(cis_best_snp_random)

me_cis_random = MatrixEQTL_wrapper(genotype_subset, expression_matrix, snp_position_subset, gene_positions, pvThresh = pvOutputThreshold, pvThresh_cis = pvOutputThreshold, cis_dist = cis_dist, covariates = covars)

# Then prepare the GWAS subsetting:
gwas_subset_dir = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/subset_SNPs/gwas-subset/')

if (file.exists(paste0(gwas_subset_dir, tissue_name, '_part', part_number, '.txt'))) {
	gwas_df = read.table(paste0(gwas_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', header=F, stringsAsFactors=F)
} else {
	gwas_catalogue = read.csv('/tigress/BEE/RNAseq/Data/Resources/annotations/gwas_catalog_v1.0.1-associations_e84_r2016-06-12.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)
	gwas_catalogue = gwas_catalogue[(gwas_catalogue$SNPS %in% rownames(genotype_matrix)),]
	gwas_snps = vector('list', nrow(gwas_catalogue))
	for (i in c(1:nrow(gwas_catalogue))) {
		snp = gwas_catalogue$SNPS[i]
		if (grepl('chr', snp)) {
			# key = paste(gwas_df$CHR_ID[i], gwas_df$CHR_POS[i], sep='_')
			# ind = which(sapply(snp_metadata$ID, function(x) {grepl(key, x)}))
			# if (length(ind) == 1) {
			# 	gwas_snps[[i]] = snp_metadata$dbSNPID[ind]
			# }
		} else {
			gwas_snps[[i]] = strsplit(snp, ', ')[[1]]
		}
	}
	gwas_snps = unlist(gwas_snps, recursive=T)
	gwas_snps = gwas_snps[gwas_snps %in% rownames(genotype_matrix)]
	gwas_snps = unique(gwas_snps)
	gwas_df = snp_metadata[gwas_snps,]
	write.table(gwas_df, file = paste0(gwas_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', quote=F, row.names=F)
}

print('gwas-subsetting tests')

genotype_subset = genotype_matrix[gwas_df$dbSNPID,]
snp_position_subset = snp_positions[gwas_df$dbSNPID,]
gwas_subset_size = nrow(gwas_df)

me_gwas = MatrixEQTL_wrapper(genotype_subset, expression_matrix, snp_position_subset, gene_positions, pvThresh = pvOutputThreshold, pvThresh_cis = pvOutputThreshold, cis_dist = cis_dist, covariates = covars)

# Also get a set of matched random variants
gwas_subset_dir = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/subset_SNPs/gwas-subset-random/')

if (file.exists(paste0(gwas_subset_dir, tissue_name, '_part', part_number, '.txt'))) {
	gwas_df_random = read.table(paste0(gwas_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', header=F, stringsAsFactors=F)
} else {
	# selection scheme - for each test SNP, select random SNP that is matched in TSS and MAF, while not being too close to any of the other SNPs
	gwas_df_random = vector('list', nrow(gwas_df))
	for (i in c(1:nrow(gwas_df))) {
		# try to match - let's say 1000 bp and 0.01 in MAF
		tss_match = which(abs(snp_metadata$TSS - abs(gwas_df$TSS[i])) < 1000)
		maf_match = which(abs(snp_metadata$MAF - gwas_df$MAF[i]) < 0.01)
		both_match = intersect(tss_match, maf_match)
		# if too stringent, relax a little bit
		if (length(both_match) == 0) {
			tss_match = which(abs(snp_metadata$TSS - abs(gwas_df$TSS[i])) < 5000)
			maf_match = which(abs(snp_metadata$MAF - gwas_df$MAF[i]) < 0.03)
			both_match = intersect(tss_match, maf_match)
		}
		# if still too stringent, skip
		if (length(both_match) == 0) {next}
		# not in LD with any of the snps being tested - let's say 10000 bp
		snp_metadata_subset = snp_metadata[both_match,]
		snp_metadata_subset = snp_metadata_subset[sapply(snp_metadata_subset$POS, function(x) {sum(abs(gwas_df$POS - x) < 10000) == 0}),]
		# if still too stringent, skip
		if (nrow(snp_metadata_subset) == 0) {next}
		# now randomly select
		set.seed = 42
		gwas_df_random[[i]] = snp_metadata_subset[sample(nrow(snp_metadata_subset), 1),]
	}

	gwas_df_random = do.call("rbind", gwas_df_random[!duplicated(sapply(gwas_df_random, function(x) {x$dbSNPID}))])
	write.table(gwas_df_random, file = paste0(gwas_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', quote=F, row.names=F)
}

print('gwas-subsetting tests, matched random')

genotype_subset = genotype_matrix[gwas_df_random$dbSNPID,]
snp_position_subset = snp_positions[gwas_df_random$dbSNPID,]
gwas_subset_size_random = nrow(gwas_df_random)

me_gwas_random = MatrixEQTL_wrapper(genotype_subset, expression_matrix, snp_position_subset, gene_positions, pvThresh = pvOutputThreshold, pvThresh_cis = pvOutputThreshold, cis_dist = cis_dist, covariates = covars)

# # Then prepare the LD subsetting:
# ld_subset_dir = paste0(proj_dir, '/Data/Genotype/gtex/imputed_genotypes/allelic_dosage/continuous/subset_SNPs/ld-subset/')
# if (file.exists(paste0(ld_subset_dir, tissue_name, '_part', part_number, '.txt'))) {
# 	ld_subset_df = read.table(paste0(ld_subset_dir, tissue_name, '_part', part_number, '.txt'), sep='\t', header=F, stringsAsFactors=F)
# 	colnames(ld_subset_df) = 'snps'
# } else {
# 	ld_subset_list = read.table(paste0(ld_subset_dir, 'master/chr', chr_number, '.txt'), header=F, stringsAsFactors=FALSE)
# 	ld_subset_list = ld_subset_list[ld_subset_list$V1 %in% rownames(genotype_matrix),'V1']
# 	ld_subset_df = data.frame(snps = ld_subset_list)
# 	write.table(ld_subset_df, file = paste0(ld_subset_dir, tissue_name, '_part', part_number, '.txt'), quote=F, row.names=F, col.names=F)
# }

# genotype_subset = genotype_matrix[ld_subset_df$snps,]
# snp_position_subset = snp_positions[ld_subset_df$snps,]
# ld_subset_size = nrow(ld_subset_df)

# me_ld = MatrixEQTL_wrapper(genotype_subset, expression_matrix, snp_position_subset, gene_positions, pvThresh = pvOutputThreshold, pvThresh_cis = pvOutputThreshold, cis_dist = cis_dist, covariates = covars)

# save(cis_subset_size, me_cis, gwas_subset_size, me_gwas, ld_subset_size, me_ld, file=paste0(out_file, '_part', sprintf("%03d", part_number), '.RData') )

save(cis_subset_size, cis_subset_size_random, me_cis, me_cis_random, gwas_subset_size, gwas_subset_size_random, me_gwas, me_gwas_random, file=paste0(out_file, '_part', sprintf("%03d", part_number), '.RData') )
