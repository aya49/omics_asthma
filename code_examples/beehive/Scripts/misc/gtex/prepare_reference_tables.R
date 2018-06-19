##################################################
#  prepare_reference_tables.R
#
#  $proj/Scripts/misc/gtex/prepare_reference_tables.R
# 
#  This script prepares the three tables with necessary information (sample-level, subject-level and tissue-level) for downstream processing.
#
#  Author: Brian Jo
#
##################################################

# There are two main annotation files that we want to use for this task:

proj_dir = Sys.getenv('proj')

sample_attr = read.csv(paste0(proj_dir, '/Data/Resources/gtex/covariates/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt'), header=T, stringsAsFactors=F, sep='\t')
srr_runtable = read.csv(paste0(proj_dir, '/Data/Resources/gtex/covariates/SRARunTable.txt'), header=T, stringsAsFactors=F, sep='\t')

sample_attr_use = sample_attr[sample_attr$SMAFRZE == 'USE ME',]
# dim(sample_attr_use)
# [1] 9109   75
sample_attr_use = sample_attr_use[sample_attr_use$ANALYTE_TYPE == 'RNA:Total RNA',]
# dim(sample_attr_use)
# [1] 8555   75

sample_table = sample_attr_use[,c('SAMPID', 'SMTSD')]
match_ind = as.numeric(sapply(sample_table$SAMPID, function(x) {which(srr_runtable$Sample_Name_s == x)}))

sample_table = cbind(sample_table, srr_runtable[match_ind, c('Run_s', 'histological_type_s', 'sex_s', 'submitted_subject_id_s')])
sample_table = sample_table[,c('Run_s', 'SAMPID', 'SMTSD', 'histological_type_s', 'submitted_subject_id_s', 'sex_s')]
# There are four entries where there is no matching to SRR ID
sample_table = sample_table[!is.na(sample_table$Run_s),]
# dim(sample_table)
# [1] 8551    6

# minor fixes:
which(sample_table$histological_type_s == "<not provided>")
sample_table$histological_type_s[7379] = 'Stomach'
sample_table$histological_type_s[7443] = 'Skin'

# convert tissue names and histological types

sample_table$tissue_name = sapply(sample_table$SMTSD, function(x) {gsub("[[:punct:]]", "", x)})
sample_table$tissue_name = sapply(sample_table$tissue_name, function(x) {tolower(gsub(" ", "", x))})

sample_table$histology = sapply(sample_table$histological_type_s, function(x) {tolower(gsub(" ", "", x))})

# which samples have genotypes?
genotype_subject_list = read.table(paste0(proj_dir, '/Data/Resources/gtex/genotype/subjects_with_genotypes.txt'), header=F, stringsAsFactors=F)

sample_table$genotype_avail = sapply(sample_table$submitted_subject_id_s, function(x) {x %in% genotype_subject_list$V1})
rownames(sample_table) = sample_table$Run_s
sample_table = sample_table[,c('Run_s', 'submitted_subject_id_s', 'genotype_avail', 'sex_s', 'SAMPID', 'tissue_name', 'histology', 'SMTSD', 'histological_type_s')]

# Write sample_level table
write.table(sample_table, sep='\t', file = paste0(proj_dir, '/Data/Resources/gtex/tables/sample_table.txt'), quote=F)

# Make subject_level table
subject_table = data.frame(submitted_subject_id_s = unique(sample_table$submitted_subject_id_s))
match_index = as.numeric(sapply(subject_table$submitted_subject_id_s, function(x) {which(sample_table$submitted_subject_id_s == x)[1]}))
subject_table$genotype_avail = sample_table$genotype_avail[match_index]
subject_table$sex_s = sample_table$sex_s[match_index]

subject_table$num_samples = as.numeric(sapply(subject_table$submitted_subject_id_s, function(x) {sum(sample_table$submitted_subject_id_s == x)[1]}))
rownames(subject_table) = subject_table$submitted_subject_id_s

# Write subject_level table
write.table(subject_table, sep='\t', file = paste0(proj_dir, '/Data/Resources/gtex/tables/subject_table.txt'), quote=F)

# Make tissue_level table
tissue_table = data.frame(tissue_name = unique(sample_table$tissue_name))
match_index = as.numeric(sapply(tissue_table$tissue_name, function(x) {which(sample_table$tissue_name == x)[1]}))

tissue_table$histology = sample_table$histology[match_index]
tissue_table$SMTSD = sample_table$SMTSD[match_index]
tissue_table$histological_type_s = sample_table$histological_type_s[match_index]

tissue_table$num_samples = as.numeric(sapply(tissue_table$tissue_name, function(x) {sum(sample_table$tissue_name == x)[1]}))

tissue_table$num_samples_with_geno = as.numeric(sapply(tissue_table$tissue_name, function(x) {sum(sample_table[sample_table$tissue_name == x,'genotype_avail'])}))

rownames(tissue_table) = tissue_table$tissue_name
tissue_table = tissue_table[order(-tissue_table$num_samples),]

# Write tissue_level table
write.table(tissue_table, sep='\t', file = paste0(proj_dir, '/Data/Resources/gtex/tables/tissue_table.txt'), quote=F)

# > sum(tissue_table$num_samples)
# [1] 8551
# > sum(tissue_table$num_samples_with_geno)
# [1] 7329
