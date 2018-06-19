#install CRAN packages, if missing
packages <- c("shiny", "shinyjs", "data.table", "dplyr", "tidyr", "ggplot2", "knitr", "markdown", "stringr","DT","seqminer")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)  
} else { print("All required CRAN packages installed")}

#install Bioconductor packages if missing
source("https://bioconductor.org/biocLite.R")
bioc <- c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db")
if (length(setdiff(bioc, rownames(installed.packages()))) > 0) {
  biocLite(setdiff(bioc, rownames(installed.packages())))  
} else { print("All required Bioconductor packages installed")}

library(shiny)  
runGitHub("LocusExplorer", "oncogenetics", launch.browser = TRUE)

library(shiny)  
runApp(launch.browser = TRUE)