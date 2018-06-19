##################################################
#  compile_gene_list.R
#
#  $proj/Scripts/misc/gtex/compile_gene_list.R
# 
#  Compile a list of genes represented in the current output
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# FDR Control of p-values from the null distribution
# original_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/'
# write_dir = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/gene_level_FDR/genes_represented/'

original_dir = paste0(proj_dir, args[1])
write_dir = paste0(proj_dir, args[2])

library(dplyr)
tissue_list = as.character(read.table(paste0(original_dir, 'tissue_list.txt'))$V1)

for (tissue in tissue_list) {
    print(tissue)
    suffix = '_part440.RData'
    list_files = list.files(path = paste0(original_dir, tissue, '/'), pattern = suffix)
    max_PEER = max(as.numeric(sapply(sapply(list_files, function(x) {strsplit(x, 'PEER')[[1]][2]}), function(y) {strsplit(y, '_')[[1]][1]})))
    list_max_PEER_files = list.files(path = paste0(original_dir, tissue, '/'), pattern = paste0('PEER', max_PEER))
    
    gene_list = list()
    i = 1
    for (item in list_max_PEER_files) {
        load(paste(original_dir, tissue, '/', item, sep=""))
        gene_list[[i]] = me$trans$eqtls[(!duplicated(me$trans$eqtls$gene)),c('gene', 'statistic', 'pvalue')]
        i = i+1
    }

    unique_gene_list = do.call('rbind', gene_list)
    unique_gene_list = unique_gene_list[order(unique_gene_list$pvalue),]
    unique_gene_list = unique_gene_list[(!duplicated(unique_gene_list$gene)),]
    write.table(unique_gene_list, file = paste0(write_dir, tissue, '_genes_present.txt'), quote = F, row.names = F)
    print(sum(unique_gene_list$pvalue < 5e-9))
}
