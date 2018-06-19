##################################################
#  v8_peer_factor_calculation.R
#
#  $proj/Scripts/processing/covariates_exploratory/v8_peer_factor_calculation.R
#
#  PEER factor calculation for GTEx v8 data, with genotype PCs
#
#  Authors: Brian Jo
#
##################################################

args <-commandArgs(TRUE)
exp_dir = args[1]
exp_mat = args[2]

