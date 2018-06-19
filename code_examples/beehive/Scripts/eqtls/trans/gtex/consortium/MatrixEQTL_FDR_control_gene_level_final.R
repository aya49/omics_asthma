eGene_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/gene_level_FDR/'
main_trans_list = read.table('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/main_consortium_trans/trans-eQTLs_FDR-0.1.txt', header=T, stringsAsFactors=F) 

tissues = list.files(path=paste0(eGene_dir, 'tissues/'))

cumulative_table = main_trans_list[,c('tissue', 'gene', 'pvalue', 'FDR', 'snps')]
cumulative_table$identifier = paste0(cumulative_table$tissue, cumulative_table$gene)
cumulative_table = cumulative_table[!(duplicated(cumulative_table$identifier)),]

post_crossmapping_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/all-by-all-PEER-increments/eqtl_list_post_crossmapping/'


BH_FDR = function(pvalues, ntests) {
    # Assumes that pvalues are ordered
    # ordered_pvals_list = pvalues[order(pvalues)]
    # alpha_hat = sapply(c(1:length(ordered_pvals_list)), function(x) {ordered_pvals_list[x] * ntests / x})

    alpha_hat = sapply(c(1:length(pvalues)), function(x) {pvalues[x] * ntests / x})
    alpha_hat[(alpha_hat >= 1.0)] = 1.0

    # Perform Benjamini-Hochberg
    alpha_hat_copy = alpha_hat
    cur_index = 1
    max_index = max(which(alpha_hat_copy<1))
    while (cur_index <= max_index) {
        temp = alpha_hat_copy[cur_index:length(alpha_hat_copy)]
        min_index = which.min(temp)
        min_alpha = min(temp)
        alpha_hat_copy[c(cur_index:(cur_index+min_index-1))] = min_alpha
        cur_index = cur_index + min_index
    }
    return(alpha_hat_copy)
}

cumulative_table$BonHack_BH = -1
addendum = list()
for (tissue in tissues) {
	tissue_name = strsplit(tissue, '\\.')[[1]][1]
	print(tissue_name)

	min_table = read.table(paste0(eGene_dir, '/tissues/', tissue), sep='\t', stringsAsFactors=F, header=T)
	post_crossmap_files = list.files(post_crossmapping_dir, pattern = tissue_name)
	post_crossmap_table = read.table(paste0(post_crossmapping_dir, post_crossmap_files[length(post_crossmap_files)-1]), sep='\t', stringsAsFactors=F, header=T)
	post_crossmap_table = post_crossmap_table[!(duplicated(post_crossmap_table$gene)),]

	if (nrow(post_crossmap_table) > 0) {
		post_crossmap_table$BonHack_BH = BH_FDR(post_crossmap_table$pvalue * 1e6, nrow(min_table))
		rownames(post_crossmap_table) = post_crossmap_table$gene

		if (sum(cumulative_table$tissue == tissue_name) > 0) {
			cumulative_table[cumulative_table$tissue == tissue_name,]$BonHack_BH = post_crossmap_table[(cumulative_table[cumulative_table$tissue == tissue_name,]$gene),]$BonHack_BH
		}
	}



	if (sum(post_crossmap_table$BonHack_BH <= 0.1) > 0) {
		for (j in c(1:sum(post_crossmap_table$BonHack_BH <= 0.1))) {
			identifier = paste0(tissue_name, post_crossmap_table$gene)[j]
			if (!(identifier %in% cumulative_table$identifier)) {
				addendum[[identifier]] = post_crossmap_table[j,c('snps', 'gene', 'pvalue', 'BonHack_BH')]
			}
		}
	}
}

addendum = do.call('rbind', addendum)
cumulative_table$eGene = cumulative_table$BonHack_BH <= 0.1
write.table(cumulative_table, file = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/gene_level_FDR/eGene_FDR_0.1.txt', row.names=F, sep=',', quote=F)
write.table(addendum, file = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/gene_level_FDR/eGene_FDR_0.1_addendum.txt', sep=',', quote=F)



# tissue_table = read.table('/tigress/BEE/RNAseq/Data/Resources/gtex/tables/tissue_table.txt', sep='\t')

# cumulative_table = main_trans_list[,c('tissue', 'gene', 'pvalue', 'FDR')]
# cumulative_table$identifier = paste0(cumulative_table$tissue, cumulative_table$gene)
# cumulative_table = cumulative_table[!(duplicated(cumulative_table$identifier)),]

# BonHack_table = list()
# for (tissue in tissues) {
# 	tissue_name = strsplit(tissue, '\\.')[[1]][1]
# 	print(tissue_name)

# 	ind = which(rownames(tissue_table) == tissue_name)
# 	table = read.table(paste0(eGene_dir, '/tissues/', tissue), sep='\t', stringsAsFactors=F, header=T)

# 	if (sum(table$qval_Uniform_Bonferroni <= 0.1) > 0) {
# 		BonHack_table[[tissue]] = table[c(1:sum(table$qval_Uniform_Bonferroni <= 0.1)),c('gene', 'pvalue', 'qval_Uniform_Bonferroni')]
# 		BonHack_table[[tissue]]$tissue = tissue_name
# 	}
# }

# BonHack_table = do.call('rbind', BonHack_table)

# BonHack_table$identifier = paste0(BonHack_table$tissue, BonHack_table$gene)

# cumulative_table$found_in_BonHack_qval_0.1 = sapply(cumulative_table$identifier, function(x) {x %in% BonHack_table$identifier})
# BonHack_table$found_in_eQTL_FDR_0.1 = sapply(BonHack_table$identifier, function(x) {x %in% cumulative_table$identifier})





# eGene_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/gene_level_FDR/'
# main_trans_list = read.table('/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/public_ftp/main_consortium_trans/trans-eQTLs_FDR-0.1.txt', header=T, stringsAsFactors=F) 

# main_trans_list$qvalue_BonHack = -1
# prev_tissue = ''

# for (i in c(1:nrow(main_trans_list))) {
# 	tissue = main_trans_list$tissue[i]
# 	if (prev_tissue != tissue) {
# 		prev_tissue = tissue
# 		table = read.table(paste0(eGene_dir, '/tissues/', tissue, '.txt'), sep='\t', stringsAsFactors=F, header=T)
# 	}
# 	ind = which(table$gene == main_trans_list$gene[i])
# 	main_trans_list$qvalue_BonHack[i] = table$qval_Uniform_Bonferroni[ind]
# }