##################################################
#  v8_MatrixEQTL_FDR_control.R
#
#  $proj/Scripts/eqtls/trans/gtex/v8_consortium/v8_MatrixEQTL_FDR_control.R
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# FDR Control of p-values from the null distribution
args = c(1:2)
args[1] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/v8/all-by-all/'
args[2] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/v8/all-by-all/'

summary_dir = args[1]
write_dir = args[2]

tissue_list = c('Adipose_Subcutaneous','Cells_EBV-transformed_lymphocytes','Skin_Sun_Exposed_Lower_leg','Testis','Thyroid','Whole_Blood')

library(dplyr)

trans_eqtl_list_5 = list()
trans_eqtl_list_10 = list()

trans_egene_list_5 = list()
trans_egene_list_10 = list()

for (tissue in tissue_list) {
    print(tissue)
    list_files = list.files(path = paste0(summary_dir, tissue, '/summary/', sep=''))
    print(length(list_files))

    # Required data to aggregate
    cumulative_me_trans_total = vector('list', length(list_files))
    all_hist_vals_trans_total = rep(0,100)
    all_n_tests_trans_total = 0
    snps_tested_total = 0

    i = 1
    for (f in list_files) {
        load(paste0(summary_dir, tissue, '/summary/', f))
        
        cumulative_me_trans_total[[i]] = me_trans
        all_hist_vals_trans_total = all_hist_vals_trans_total + hist_counts
        all_n_tests_trans_total = all_n_tests_trans_total + total_ntests
        snps_tested_total = snps_tested_total + snps_tested
        i = i+1
    }

    cumulative_me_trans_total = do.call("rbind", cumulative_me_trans_total)
    print(dim(cumulative_me_trans_total))

    cumulative_me_trans_total = cumulative_me_trans_total[order(cumulative_me_trans_total$pvalue),]
    # BH procedure for SNP-gene FDR:
    cumulative_me_trans_total$FDR = p.adjust(cumulative_me_trans_total$pvalue, method = 'BH', n = all_n_tests_trans_total)

    # Save all pairs with pval < 1e-5
    if(!file.exists(paste0(write_dir, 'eqtl_list_pval_1e-5/', tissue, '.RData'))) {
        save(cumulative_me_trans_total, file = paste0(write_dir, 'eqtl_list_pval_1e-5/', tissue, '.RData'))
    }

    # For pair FDR:
    n_eqtls_5 = sum(cumulative_me_trans_total$FDR <= 0.05)
    n_eqtls_10 = sum(cumulative_me_trans_total$FDR <= 0.10)
    print(n_eqtls_5)
    print(n_eqtls_10)

    print(length(unique(cumulative_me_trans_total$gene[1:n_eqtls_5])))
    print(length(unique(cumulative_me_trans_total$gene[1:n_eqtls_10])))
    print(length(unique(cumulative_me_trans_total$snps[1:n_eqtls_5])))
    print(length(unique(cumulative_me_trans_total$snps[1:n_eqtls_10])))

    trans_eqtl_list_5[[tissue]] = cumulative_me_trans_total[1:n_eqtls_5,]
    trans_eqtl_list_5[[tissue]]$tissue = tissue
    trans_eqtl_list_10[[tissue]] = cumulative_me_trans_total[1:n_eqtls_10,]
    trans_eqtl_list_10[[tissue]]$tissue = tissue

    # Get the extreme values for eGene FDR:
    cumulative_egene_df = cumulative_me_trans_total[(!duplicated(cumulative_me_trans_total$gene)),]
    cumulative_egene_df$BonHack_FDR = p.adjust((cumulative_egene_df$pvalue*1000000), method = 'BH', n = gene_tested)

    n_egene_5 = sum(cumulative_egene_df$BonHack_FDR <= 0.05)
    n_egene_10 = sum(cumulative_egene_df$BonHack_FDR <= 0.10)

    trans_egene_list_5[[tissue]] = cumulative_egene_df[1:n_egene_5,]
    trans_egene_list_5[[tissue]]$tissue = tissue
    trans_egene_list_10[[tissue]] = cumulative_egene_df[1:n_egene_10,]
    trans_egene_list_10[[tissue]]$tissue = tissue
}

trans_eqtls_5 = do.call("rbind", trans_eqtl_list_5)
trans_eqtls_10 = do.call("rbind", trans_eqtl_list_10)

trans_egene_5 = do.call("rbind", trans_egene_list_5)
trans_egene_10 = do.call("rbind", trans_egene_list_10)

write.table(trans_eqtls_5, file = paste0(write_dir, 'v8_trans_eqtls_FDR_0.05.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(trans_eqtls_10, file = paste0(write_dir, 'v8_trans_eqtls_FDR_0.10.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(trans_egene_5, file = paste0(write_dir, 'v8_trans_egenes_FDR_0.05.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(trans_egene_10, file = paste0(write_dir, 'v8_trans_egenes_FDR_0.10.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
