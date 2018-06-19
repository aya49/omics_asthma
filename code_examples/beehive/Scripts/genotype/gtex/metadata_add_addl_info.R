##################################################
#  metadata_add_addl_info.R
#
#  /tigress/BEE/RNAseq/Scripts/genotype/gtex/metadata_add_addl_info.R
# 
#  This script creates the metadata files that also have the information about the cis genes and repeat regions
#
#  Author: Brian Jo
#
##################################################

setwd('/tigress/BEE/eQTLs/Data/Genotype/GTEx/SNP_metadata')

# Load in RepeatMasker annotations
load('/tigress/BEE/eQTLs/Data/References/Regulatory/RepeatMasker/repeat_masker.RData')

# which chromsome are we looking at?
args <-commandArgs(TRUE)

i = args[1]
in_file = args[2]
out_file = args[3]

SNP_metadata = read.table(in_file, skip = 6, header=TRUE, stringsAsFactors=FALSE)

# cut down repeat_masker
repeat_masker = repeat_masker[repeat_masker$V5 == paste('chr',i,sep=''),]
total = dim(repeat_masker)[1]

ptm = proc.time()
repeat_element_ind = sapply(SNP_metadata$POS, function(x) {
	ind = 0
	if ((sum(repeat_masker$V6 < x) + sum(repeat_masker$V7 > x)) > total) {
		# Probably safe with just taking the sum, but just to be safe
		ind = which(repeat_masker$V6 < x)[sum(repeat_masker$V6 < x)]
	}
	return(ind)
})
print(proc.time() - ptm)

SNP_metadata[,'IN_REPEAT'] = repeat_element_ind>0
SNP_metadata[,'REPEAT_ELEMENT'] = NA
SNP_metadata[,'REPEAT_FAMILY'] = NA

SNP_metadata[,'REPEAT_ELEMENT'][which(SNP_metadata$IN_REPEAT)] = sapply(which(SNP_metadata$IN_REPEAT), function(x) { as.character(repeat_masker$V10[repeat_element_ind[x]]) })
SNP_metadata[,'REPEAT_FAMILY'][which(SNP_metadata$IN_REPEAT)] = sapply(which(SNP_metadata$IN_REPEAT), function(x) { as.character(repeat_masker$V11[repeat_element_ind[x]]) })

# Now that we have added in repeat element information, now let's add in cis-gene information
# caveat: this just takes in all protein-coding annotations. It doesn't take into account overlapping nature of the genes

gene_list = read.table('/tigress/BEE/eQTLs/Data/References/GTEx/gene_names_positions.txt', header=TRUE, stringsAsFactors=FALSE)
gene_list = gene_list[(gene_list$chr == i),]
gene_list = gene_list[(gene_list$gene_type == 'protein_coding'),]

# For now, let's say cis- distance is 150kb
cis_dist = 150000
cis_gene_list = sapply(SNP_metadata$POS, function(x) {
	list = NA
	ind1 = x > (gene_list$start - cis_dist)
	ind2 = x < (gene_list$end + cis_dist)
	if (length(which(ind1 & ind2)) > 0) {
		list = paste(gene_list$gene_name[which(ind1 & ind2)], collapse=',')
	}
	return(list)
})

in_gene_list = sapply(SNP_metadata$POS, function(x) {
	list = NA
	ind1 = x > gene_list$start
	ind2 = x < gene_list$end
	if (length(which(ind1 & ind2)) > 0) {
		list = paste(gene_list$gene_name[which(ind1 & ind2)], collapse=',')
	}
	return(list)
})

SNP_metadata[,'CIS-GENE'] = cis_gene_list
SNP_metadata[,'IN-GENE'] = in_gene_list

# Now we can write out the updated metadata files:
write.table(SNP_metadata, file=out_file, quote=FALSE, sep='\t', na = "NA", row.names = FALSE)