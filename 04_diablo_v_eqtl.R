## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720



## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/visualizationFunctions.R"))
libr(append(pkgs(),c("qdapTools")))


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
