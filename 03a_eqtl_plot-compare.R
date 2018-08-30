## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
source(paste0(root, "/code/visualizationFunctions.R"))
libr(append(pkgs(),c("qdapTools")))

no_cores = 8#detectCores()-3
registerDoMC(no_cores)


## options
overwrite = T
writecsv = T #write results as csv on top of Rdata?

categorical = T # is class column categorical?
interested_cols = c("response")
# interested_cols = c("age","bmi","sex","centre","batch","race","response") 
# interested_cont_cols = ""

useModels = c("modelLINEAR", "modelANOVA", "modelLINEAR_CROSS")

pthres = .01
max_plots = 150 # max number of plots to make

# plotting size
width = 500
height = 300



## plot eQTL ---------------------------------------------

start = Sys.time()

pval_paths = list.files(gwas_dir, full.names=T, pattern=".Rdata")
eqtl_paths = list.files(eqtl_dir, full.names=T, pattern=".Rdata"); eqtl_paths = eqtl_paths[grepl(".post-",eqtl_paths)]
eqtl_names = fileNames(eqtl_paths)
feat_type_sets = lapply(str_split(eqtl_names,"_"), function(x) c(str_split(x,"[-]")[[1]][1], str_split(x,"[-]")[[1]][2]))


# a = foreach (ei = 1:length(eqtl_paths)) %dopar% { try({
a = llply(1:length(eqtl_paths), function(ei) { try({
  # a = foreach (ei = 1:length(eqtl_paths)) %dopar% {
  # a = foreach (ei = 1:length(eqtl_paths)) %dopar% {
  eqtl_path = eqtl_paths[ei]
  eqtl_name = eqtl_names[ei]
  feat_type_set = feat_type_sets[[ei]]
  if (!grepl(".post",feat_type_set[2])) next()
  
  file_ind = unlist(read.csv(gsub(".Rdata","_id.csv",eqtl_path)))
  file_ind_n = gsub("[-]|X","",str_extract(eqtl_name,"-[a-zA-Z]+X"))
  
  # load eqtl / matrix
  mecisl = NULL
  for (met in c("pre","post","diff")) { try({ 
    fn = gsub(feat_type_set[2], gsub("post",met,feat_type_set[2]), eqtl_path)
    if (!file.exists(fn)) next
    mecis = get(load(fn)) 
    mecis$snps = as.character(mecis$snps)
    mecis$gene = as.character(mecis$gene)
    mecis$id_temp = paste(mecis$snps, mecis$gene)
    mecisl[[met]] = mecis
  }) }
  if (length(mecisl)<2) next()
  break()
  
  for (metcomp in list(c("pre","diff"),c("pre","post"))) {
    idint = Reduce(intersect,lapply(metcomp, function(xi) mecisl[[xi]]$id_temp))
    if (length(idint)==0) next()
    mecisll = lapply(metcomp, function(xi) mecisl[[xi]][match(idint,mecisl[[xi]]$id_temp),])
    names(mecisll) = metcomp
    mecislt = Reduce(merge,mecisll)
    write.csv(mecislt, file=gsub(".Rdata",paste0("_",paste(metcomp,collapse="VS"),".csv"),
                                 gsub(feat_type_set[2], gsub("post",met,feat_type_set[2]), eqtl_path)))
    
  }
}) }, .progress="tk", .parallel=T)
time_output(start)










