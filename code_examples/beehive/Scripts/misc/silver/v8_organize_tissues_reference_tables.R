##################################################
#  v8_organize_tissues_reference_tables.R
#
#  $proj/Scripts/misc/silver/v8_organize_tissues_reference_tables.R
# 
#  This script prepares the three tables with necessary information (sample-level, subject-level and tissue-level) for downstream processing.
#
#  Author: Brian Jo
#
##################################################

# There are two main annotation files that we want to use for this task:

proj_dir = Sys.getenv('proj')

subject_attr = read.csv(paste0(proj_dir, '/Data/Resources/gtex/covariates/v8/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt'), header=T, stringsAsFactors=F, sep='\t')
sample_attr = read.csv(paste0(proj_dir, '/Data/Resources/gtex/covariates/v8/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'), header=T, stringsAsFactors=F, sep='\t')
# R takes in the first column name as 'X'
rownames(sample_attr) = sample_attr$X

# dim(sample_attr)
# [1] 17395    74

# srr_runtable = read.csv(paste0(proj_dir, '/Data/Resources/gtex/covariates/SRARunTable.txt'), header=T, stringsAsFactors=F, sep='\t')
# sample_attr_use = sample_attr[sample_attr$SMAFRZE == 'USE ME',]
# sample_attr_use = sample_attr_use[sample_attr$ANALYTE_TYPE == 'RNA:Total RNA',]

# R takes in the first column name as 'X'
sample_table = sample_attr[,c('X', 'SMTSD')]
colnames(sample_table)[1] = 'SAMPID'
sample_table$SUBJID = sapply(sample_table$SAMPID, function(x) {paste(strsplit(x, '-')[[1]][1], strsplit(x, '-')[[1]][2], sep = '-')})
match_ind = as.numeric(sapply(sample_table$SUBJID, function(x) {which(subject_attr$SUBJID == x)}))
sample_table = cbind(sample_table, subject_attr[match_ind, c('SEX', 'RACE', 'AGE', 'BMI')])

# convert tissue names - update for v8: Now we just remove punctuations, and replace whitespace with underscore
sample_table$tissue_name = sapply(sample_table$SMTSD, function(x) {gsub("[[:punct:]]", "", x)})
# sample_table$tissue_name = sapply(sample_table$tissue_name, function(x) {tolower(gsub(" ", "", x))})
sample_table$tissue_name = sapply(sample_table$tissue_name, function(x) {gsub(" ", "_", (gsub("  ", " ", x)))})

# sample_table$histology = sapply(sample_table$histological_type_s, function(x) {tolower(gsub(" ", "", x))})

# All samples have genotypes in v8
# genotype_subject_list = read.table(paste0(proj_dir, '/Data/Resources/gtex/genotype/subjects_with_genotypes.txt'), header=F, stringsAsFactors=F)

# Write sample_level table
write.table(sample_table, sep='\t', file = paste0(proj_dir, '/Data/Resources/gtex/tables/v8/sample_table.txt'), quote=F, row.names=F)

# Make subject_level table
# subject_table = data.frame(submitted_subject_id_s = unique(sample_table$submitted_subject_id_s))
# match_index = as.numeric(sapply(subject_table$submitted_subject_id_s, function(x) {which(sample_table$submitted_subject_id_s == x)[1]}))
# subject_table$genotype_avail = sample_table$genotype_avail[match_index]
# subject_table$sex_s = sample_table$sex_s[match_index]

subject_table = subject_attr[,c('SUBJID', 'SEX', 'AGE', 'RACE', 'ETHNCTY', 'BMI')]
subject_table$num_samples = as.numeric(sapply(subject_table$SUBJID, function(x) {sum(sample_table$SUBJID == x)[1]}))
rownames(subject_table) = subject_table$SUBJID

# Write subject_level table
write.table(subject_table, sep='\t', file = paste0(proj_dir, '/Data/Resources/gtex/tables/v8/subject_table.txt'), quote=F, row.names=F)

# Make tissue_level table
tissue_table = data.frame(tissue_name = unique(sample_table$tissue_name), stringsAsFactors=F)
match_index = as.numeric(sapply(tissue_table$tissue_name, function(x) {which(sample_table$tissue_name == x)[1]}))

tissue_table$SMTS = sample_attr$SMTS[match_index]
tissue_table$SMUBRTRM = sample_attr$SMUBRTRM[match_index]

# tissue_table$histology = sample_table$histology[match_index]
# tissue_table$SMTSD = sample_table$SMTSD[match_index]
# tissue_table$histological_type_s = sample_table$histological_type_s[match_index]

tissue_table$num_samples = as.numeric(sapply(tissue_table$tissue_name, function(x) {sum(sample_table$tissue_name == x)[1]}))

# tissue_table$num_samples_with_geno = as.numeric(sapply(tissue_table$tissue_name, function(x) {sum(sample_table[sample_table$tissue_name == x,'genotype_avail'])}))
# rownames(tissue_table) = tissue_table$tissue_name
tissue_table = tissue_table[order(-tissue_table$num_samples),]

# Write tissue_level table
write.table(tissue_table, sep='\t', file = paste0(proj_dir, '/Data/Resources/gtex/tables/v8/tissue_table.txt'), quote=F, row.names=F)

