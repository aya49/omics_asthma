## input: GWAS & EQTL
## output: visualize GWAS & EQTL with coMET
## aya43@sfu.ca
## created 20180630



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result")


# asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")

stat_dir = paste0(result_dir,"/stat")
eqtl_dir = paste0(stat_dir,"/eqtl")
gwas_dir = paste0(stat_dir,"/gwas")


## output directory



## libraries
# source("https://bioconductor.org/biocLite.R")
# install.packages('shiny', repos='http://cran.rstudio.com/')
# install.packages('rmarkdown', repos='http://cran.rstudio.com/')
source("code/_func.R")
library("psych")
library("Gviz")
library("corrplot")
library("coMET")

libr("biomaRt")
libr("rtracklayer")

libr("data.table")
libr("MatrixEQTL")
libr("foreach")
libr("doMC")
libr("stringr")
libr("qdapTools") # make dummy variables
libr("Matrix")

## 3 main functions:
#
# 1. comet.web: pre-customized function that allows us to visualise quickly
# EWAS (or results/enrichr omic-WAS) results, annotation tracks, and correlations 
# between features. This version is installed in the Shiny web-service. 
# Currently, it is formated only to visualise human data.
# 
# 2. comet: generic function that allows us to 
# visualise quickly EWAS results, annotation tracks, and correlations 
# between features. Users can visualise more personalised annotation tracks 
# and give multiple extra EWAS/omic-WAS results to plot.
# 
# 3. comet.list: additional function that allows us to extract the values of 
# correlations, the pvalues, and estimates and confidence intervals 
# for all datapoints that surpass a particular threshold.

## five types of files that can be given by the user to produce the plot:
# Warning: These are mandatory and has to be in tabular format with a header
# 1. Info: defined in the option mydata.file
# 
# 2. Correlation: defined in the option cormatrix.file
# 
# 3. Extra info: are defined in the option mydata.file.large
# 
# 4. Annotation info file is defined in the option biofeat.user.file
#    exists only in the function comet.web
#    and should inform also the format to visualise data with the options
#    biofeat.user.type and biofeat.user.type.plot
# 
# 5. Configuration file contains the values of these options 
#    instead of defining these by command line.
#    Warning: Each line in the file is one option. The name of
#    the option is in capital letters and is separated by its value by "="
#    If there are multiple values such as for the option list.tracks
#    or the options for additional data, you need to separated them by a "comma"

# Info mydata.file: name, chr, (region/region_asso: start, end) (site/site_asso: pos), p, (site_asso/region_asso: sign/ass effect direction)

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
infofile <- file.path(extdata, "cyp1b1_infofile.txt")
data_info <- read.csv(infofile, header = TRUE, sep = "\t", quote = "")
# TargetID: "cg22248750", CHR, MAPINFO: 38294160, Pval: 0.027

infoexp <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
data_infoexp <- read.csv(infoexp, header = TRUE, sep = "\t", quote = "")
# TargetID: "ENSG00000138061.7_38294652_38298453", CHR, MAPINFO.START/STOP: 38294160, Pval: 0.027, BETA: +/-

# Correlation cormatrix.file: 
# 
# 1. cormatrix: pre-computed correlation matrix provided by the user; Dimension
# of matrix : CpG_number X CpG_number. Need to put the CpG sites/regions
# in the ascending order of positions and to have a header with the name of CpG
# sites/regions;
# 
# 2. raw: Raw data format. Correlations of these can be computed by one of 3
# methods Spearman, Pearson, Kendall (option cormatrix.method
# ). Dimension of matrix : sample_size X CpG_number. Need to have a header with the name
# of CpG sites/regions ;
# 
# 3.raw_rev: Raw data format. Correlations of these can be computed by one of 3
# methods Spearman, Pearson, Kendall (option cormatrix.method). Dimension
# of matrix : CpG_number X sample_size. Need to have the row names of CpG
# sites/regions and a header with the name of samples ;

# extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
corfile <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
data_cor <- read.csv(corfile, header = TRUE, sep = "\t", quote = "")
# "cg19753864" .. vs "356" ...

# config file for plot

# extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver_Grch38.txt")
data_config <- read.csv(configfile, quote = "", sep="\t", header=FALSE)


## plot via the coMET website (http://epigen.kcl.ac.uk/comet).

# extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
myinfofile <- file.path(extdata, "cyp1b1_infofile_Grch38.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region_Grch38.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver_Grch38.txt")
comet.web(config.file=configfile, mydata.file=myinfofile,
          cormatrix.file=mycorrelation ,mydata.large.file=myexpressfile,
          print.image=FALSE,verbose=FALSE)



## plot via coMET

# pvalue plot, annotation tracks, and correlation matrix
# annotation tracks by Gviz, track viewer or ggbio


# extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
# mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")

chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg38"
strand <- "*"

#??
# BROWSER.SESSION="UCSC"
# mySession <- browserSession(BROWSER.SESSION)
# genome(mySession) <- gen

genetrack <- genes_ENSEMBL(gen,chrom,start,end,showId=TRUE)
# snptrack <- snpBiomart_ENSEMBL(gen,chrom, start, end,
#                                dataset="ENSEMBL_MART_SNP",showId=FALSE)
# cpgIstrack <- cpgIslands_UCSC(gen,chrom,start,end)
# prombedFilePath <- file.path(extdata, "/RoadMap/regions_prom_E063.bed")
# promRMtrackE063 <- DNaseI_RoadMap(gen,chrom,start, end, prombedFilePath,
#                                  featureDisplay='promotor', type_stacking="squish")
# bedFilePath <- file.path(extdata, "RoadMap/E063_15_coreMarks_mnemonics.bed")
# chromHMM_RoadMapAllE063 <- chromHMM_RoadMap(gen,chrom,start, end,
#                                             bedFilePath, featureDisplay = "all", colorcase='roadmap15')

#Data no more available in UCSC (from September 2015)
# iscatrack <- ISCA_UCSC(gen,chrom,start,end,mySession, table="iscaPathogenic")
# listgviz <- list(genetrack,snptrack,iscatrack)
# listgviz <- list(genetrack,snptrack,cpgIstrack,promRMtrackE063,chromHMM_RoadMapAllE063)

matrix.dnamethylation <- read.delim(myinfofile, header=TRUE, sep="\t", as.is=TRUE,
                                    blank.lines.skip = TRUE, fill=TRUE)
matrix.expression <- read.delim(myexpressfile, header=TRUE, sep="\t", as.is=TRUE,
                                blank.lines.skip = TRUE, fill=TRUE)
cormatrix.data.raw <- read.delim(mycorrelation, sep="\t", header=TRUE, as.is=TRUE,
                                 blank.lines.skip = TRUE, fill=TRUE)

listmatrix.expression <- list(matrix.expression)
listcormatrix.data.raw <- list(cormatrix.data.raw)

comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
      cormatrix.file=mycorrelation,  cormatrix.type="listfile",
      mydata.large.file=myexpressfile, mydata.large.type="listfile",
      # tracks.gviz=listgviz, 
      verbose=F, print.image=F)



# if(interactive()){
#   cat("interactive")
#   genetrack <-genesENSEMBL(gen,chrom,start,end,showId=TRUE)
#   # snptrack <- snpBiomart(chrom, start, end,
#   #                        dataset="hsapiens_snp_som",showId=FALSE)
#   # strutrack <- structureBiomart(chrom, start, end,
#   #                               strand, dataset="hsapiens_structvar_som")
#   clinVariant<-ClinVarMainTrack(gen,chrom,start,end)
#   clinCNV<-ClinVarCnvTrack(gen,chrom,start,end)
#   gwastrack <-GWASTrack(gen,chrom,start,end)
#   geneRtrack <-GeneReviewsTrack(gen,chrom,start,end)
#   # listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
#   #                  clinCNV,gwastrack,geneRtrack)
#   comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
#         cormatrix.file=mycorrelation, cormatrix.type="listfile",
#         mydata.large.file=myexpressfile, mydata.large.type="listfile",
#         tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE,disp.pvalueplot=FALSE)
# } else {
#   cat("Non interactive")
#   data(geneENSEMBLtrack)
#   data(snpBiomarttrack)
#   data(ISCAtrack)
#   data(strucBiomarttrack)
#   data(ClinVarCnvTrack)
#   data(clinVarMaintrack)
#   data(GWASTrack)
#   data(GeneReviewTrack)
#   # listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
#   #                  clinCNV,gwastrack,geneRtrack)
#   comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
#         cormatrix.file=mycorrelation, cormatrix.type="listfile",
#         mydata.large.file=myexpressfile,  mydata.large.type="listfile",
#         # tracks.gviz=listgviz, 
#         verbose=T, print.image=T,disp.pvalueplot=T)
# }






