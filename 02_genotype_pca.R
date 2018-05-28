## input: features & meta_file
## output: p vaue features testing significance of ER DR correlation
## aya43@sfu.ca
## created 20180524
## last modified 20180524



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result/genotype")


asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col",asthma)

feat_dir = paste0(result_dir,"/feat")
feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype",asthma)


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)


## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
library(data.table)
library(MatrixEQTL)
library(RDRToolbox)
library(entropy)
library(foreach)
library(doMC)
library(stringr)
library(Matrix)
source(paste0(root,"/codes/_func.R"))



## options
no_cores = 15#detectCores()-3
registerDoMC(no_cores)

writecsv = T #write results as csv on top of Rdata?

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis

id_col = "well"
class_col = "response"
categorical = T # is class column categorical?
interested_cols = c("sex","centre","batch","race","response") #"centre"
interested_cont_cols = ""

dofullPCA = F #do PCA for all cell popoulations not just first level
pca_k = 5 #number of pca pc to plot

doISO = T #do ISO feature reduction
iso_k = 3 # number of neighbours
mds_type = c("iso", "mds")

doTsne = T #do Tsne feature reduction
theta=.5 #parameter for Tsne

doHC = F #do hierarchical clustering













start = Sys.time()


meta_file = get(load(paste0(meta_file_dir,".Rdata")))
meta_col = as.data.frame(get(load(paste0(meta_col_dir,".Rdata"))))
m = get(load(paste0(feat_genotype_dir,".Rdata")))
if (sum(colnames(m)%in%meta_file[,id_col])>=ncol(m)) m = t(m)

good_col_ind = apply(m,2,function(x) sum(is.na(x))<good_col)
m = m[,good_col_ind]; m[is.na(m)] = -1
meta_col = meta_col[good_col_ind,]


## get interested columns
meta_file_factor = as.data.frame(sapply(meta_file, function(x) as.numeric(factor(x, ordered=T))))

interested_col_ind = which(colnames(meta_file)%in%interested_cols)
uniquecols = apply(meta_file_factor, 2, function(x) nrow(meta_file_factor)==length(unique(x)))
meta_file_attonly = meta_file_factor[,!uniquecols]



## pca analysis ------------------------------------------------
if (doISO) { cat("iso; ")
  fit <- Isomap((m),k=iso_k)
}
cat("pca; ")
pc <- prcomp((m))

#pca scatterplot
pngname = paste0(gsub("feat/","stat/",feat_genotype_dir),"_pca-iso.png")
png(filename=pngname, width=length(interested_col_ind)*400, height=(1+doISO+pca_k)*400)
layout(matrix(c(rep(1,length(interested_col_ind)), 2:(((2*pca_k)+doISO)*length(interested_col_ind)+1)),ncol=length(interested_col_ind),byrow=T))
par(cex=1)

plot(pc$sdev^2/sum(pc$sdev^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
for (i in 1:pca_k) {
  for (col in interested_col_ind) {
    colname = colnames(meta_file)[col]
    coloursm = meta_file_factor[,col]
    plot(pc$x[,i], pc$x[,i+1], col = coloursm, main = paste0("PCA ", colname), xlab = paste0("PC_",i), ylab = paste0("PC_",i+1))
    points(0, 0, pch = 3, cex = 4, lwd = 4)
  }
  for (col in interested_col_ind) {
    colname = colnames(meta_file)[col]
    coloursm = meta_file_factor[,col]
    
    attribute = meta_file_factor[,col]
    attributen = meta_file[,col]
    attributenames = sort(unique(attributen))
    
    cor = cor(attribute, pc$x[,i])
    
    if (colname%in%interested_cont_cols) {
      plot(attribute, pc$x[,i], col=coloursm, main=paste0("PCA ", colname," Pearson Cor = ", cor), xlab = colname, ylab = paste0("PC_",i))
    } else {
      xy_boxplot = lapply(attributenames, function(x) pc$x[attributen==x,i])
      boxplot(xy_boxplot, lwd = 1, outline=F, ylim=c(min(pc$x[,i]),max(pc$x[,i])),
              main = paste0("PCA ", colname," Pearson Cor (ok if 2 var values) = ", round(cor,3)),
              xaxt = 'n', xlab=colname, ylab = paste0("PC_",i)) #,xaxt ='n', xlab=testcol
      axis(1, at=1:length(attributenames), labels=attributenames, las=2)
      points(jitter(attribute, factor=1), pc$x[,i], col = coloursm)
    }
  }
}
if (doISO) {
  for (col in interested_col_ind) {
    coloursm <- meta_file_factor[,col]
    plot(fit$dim2, col = coloursm, main = paste0("ISO ", colname), las=2)
    points(0, 0, pch = 3, cex = 4, lwd = 4)
  }
}

graphics.off()


time_output(start)



