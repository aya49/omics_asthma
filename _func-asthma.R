## list of commonly used packages for this pipeline + options
## aya43@sfu.ca
## created 20180509
## last modified 20180816


## HOW TO USE: source this after sourcing _func.R


## packages ---------------------------------------------
pkgs <- function() 
  c("gwascat", # interface to the [NHGRI's](http://www.genome.gov/) database of gwas
    "entropy",
    "edgeR", "limma",
    "data.table", "Matrix", "gdata", #read xls
    "foreach", "doMC",
    "stringr", "pracma",
    "ggplot2", "reshape2","plyr", "dplyr", "tidyr", 
    "devtools")

## input: list of package names to load
## output: none; load/install package
## note: customized for asthma pipeline
libr <- function(pkgs) {
  # install devtools to install packages from github
  if (!"devtools"%in%installed.packages())
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  
  if ("amritr"%in%pkgs & !"amritr"%in%installed.packages())
    devtools::install_github("singha53/amritr")
  
  if ("sear"%in%pkgs & !"sear"%in%installed.packages())
    devtools::install_github('cashoes/sear')

  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) 
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}


## options -----------------------------------------
options(stringsAsFactors=F)
options(.parallel=T)
id_col = "id"
date_col = "date"
class_col = "response"
controls = "ER"
experiments = "DR"
flipper_col = "flipper_calc"
good_col = 3 #each gene must have > good_col samples of each genotype; else delete
good_col_nap = .5 #each gene must have > good_col_na x 100% samples with data
good_catecont = .75
good_count = 10 #each gene must have >10 abundence in more than half the samples; else delete
good_na = .75 #proportion of na more than this, then delete the column in matrix



## load things --------------------------------------

## load meta
if (file.exists(paste0(meta_file_dir,".Rdata")))
  meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))


## load feature names
if (length(list.files(feat_dir, pattern=".Rdata", full.names=F))>0) {
  feat_types_annot = gsub(".Rdata","",list.files(feat_dir, pattern=".Rdata", full.names=F))
  feat_types_annot = feat_types_annot[order(file.size(paste0(feat_dir,"/",feat_types_annot)))]
  feat_types_annot = feat_types_annot[!grepl(".raw",feat_types_annot)]
  feat_types_times = sapply(strsplit(feat_types_annot,"[.]"), function(x) x[2])
  feat_types_names = sapply(strsplit(feat_types_annot,"[.]"), function(x) x[1])
  
  feat_types_annots = list()
  for (ftai in 1:(length(feat_types_annot)-1)) {
    for (ftaj in (ftai+1):length(feat_types_annot)) {
      fta = feat_types_annot[c(ftai,ftaj)]
      ftt = feat_types_times[c(ftai,ftaj)]
      ftn = feat_types_names[c(ftai,ftaj)]
      if (all(grepl("rna",ftn)) | all(grepl("cell",ftn)) | ftn[1] == ftn[2]) next()
      if (all(!is.na(ftt)) & ftt[1]!=ftt[2]) next()
      if (grepl("rna",ftn[1]) | grepl("genotype",ftn[2])) fta = c(fta[2],fta[1])
      feat_types_annots[[1+length(feat_types_annots)]] = fta
    }
  }
}


## features and indices
inds_paths = list.files(meta_dir,pattern="_id_",full.names=T)
inds_temp = str_split(gsub(".Rdata","",fileNames(inds_paths)),"_")
inds_types = sapply(inds_temp, function(x) x[3])
inds_names = sapply(inds_temp, function(x) str_split(x[1],"[-]")[[1]][2])
inds_filecol = sapply(inds_temp, function(x) str_split(x[1],"[-]")[[1]][1])

uif = unique(inds_names[inds_filecol=="col"])
col_inds0 = lapply(uif, function(x) list(all=c("")))
names(col_inds0) = uif
file_inds = list(all=c(""))


## load file indices
inds_paths = list.files(meta_dir,pattern="_id_",full.names=T)
inds_temp = str_split(gsub(".Rdata","",fileNames(inds_paths)),"_")
inds_types = sapply(inds_temp, function(x) x[3])
inds_names = sapply(inds_temp, function(x) str_split(x[1],"[-]")[[1]][2])
inds_filecol = sapply(inds_temp, function(x) str_split(x[1],"[-]")[[1]][1])

uif = unique(inds_names[inds_filecol=="col"])
col_inds0 = lapply(uif, function(x) list(all=c("")))
names(col_inds0) = uif
file_inds = list(all=c(""))

# col inds separate by feature type, file inds don't
for (i in 1:length(inds_paths)) {
  if (inds_filecol[i]=="col") {
    col_inds0[[inds_names[i]]][[inds_types[i]]] = get(load(inds_paths[i]))
  } else if (inds_filecol[i]=="file") {
    file_inds[[inds_types[i]]] = get(load(inds_paths[i]))
  }
}
if (exists("class_col") & exists("controls") & exists("flipper_col"))
  file_inds$flipperdr = meta_file0$id[meta_file0[,class_col]==controls | meta_file0[,flipper_col]] #flippers = DR, no actual DR's used
