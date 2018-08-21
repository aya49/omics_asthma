## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720



## logistics
root = "~/projects/asthma"
no_cores = 5#detectCores()-3
commandArgs <- function(...) list(root,F,c("qdapTools","GGally","ROCR","rafalib","pROC","CellCODE"))  # devtools::install_github("mchikina/CellCODE")
flipper_col = "flipper_calc"
class_col = "response"
controls = "ER"
experiments = "DR"
source(paste0(root, "/code/_prep.R"))


## options
class_col = "response"

pthres = .01

overwrite = F

height = 400
width = 500

good_na = .75 #proportion of na more than this, then delete the column in matrix
good_col = 3

ncomp = 2


for (feat_type in diab) {

} #feat_type
