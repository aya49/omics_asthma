## input: rmd file paths
## output: rendered rmd's
## aya43@sfu.ca


if (!"rmarkdown" %in% rownames(installed.packages())) install.packages("rmarkdown", verbose=F)

rmd_paths = list.files(paste0(root,"/code"),full.names=T,pattern=".Rmd$")
rmd_paths = rmd_paths[!grepl("_.Rmd",rmd_paths)]

for (rp in rmd_paths) rmarkdown::render(rp)
