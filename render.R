## input: rmd file paths
## output: rendered rmd's
## aya43@sfu.ca


root = "~/projects/asthma" # root directory, used for _dirs.R
rmd_paths = list.files(paste0(root,"/src"),full.names=T,pattern=".Rmd$")
rmd_paths = rmd_paths[!grepl("_.Rmd",rmd_paths)]

if (!"rmarkdown" %in% rownames(installed.packages())) install.packages("rmarkdown", verbose=F)
for (rp in rmd_paths) rmarkdown::render(rp)
