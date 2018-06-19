##################################################
#  v8_sfamix_invert_matrix.R
#
#  $proj/Scripts/processing/covariates_exploratory/v8_sfamix_invert_matrix.R
#
#  SFAmix takes in the inverted matrix - this script inverts the input matrix
#
#  Authors: Brian Jo
#
##################################################

args <-commandArgs(TRUE)
exp_dir = args[1]
exp_mat = args[2]

z = read.table(paste0(exp_dir, exp_mat), header=T, sep='\t', row.names=1)
inv_z = t(z)
write.table(inv_z, paste0('/scratch/gpfs/bj5/', exp_mat), row.names=F, col.names=F, sep='\t')
