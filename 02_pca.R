## input: features & meta_file
## output: plot PCA plot to fin dconfounders
## aya43@sfu.ca
## created 20180524
## last modified 20180524

## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result") #"genotype", "RNAseq_genes", "RNAseq_isoforms"


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
pca_dir = paste0(stat_dir,"/pca_plot"); dir.create(pca_dir, showWarnings=F)


## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("data.table")
libr("MatrixEQTL")
libr("RDRToolbox")
libr("entropy")
libr("foreach")
libr("doMC")
libr("stringr")
libr("Matrix")



## options
no_cores = 15#detectCores()-3
registerDoMC(no_cores)

writecsv = T #write results as csv on top of Rdata?

good_col = 3 #each SNP must have <good_col NA or -1; else delete from analysis
cont_col = 15 #if a column has more than this unique values, it's continuous
scale_cont = T #if a column is continuous, scale it


id_col = "id"
class_col = "response"
# categorical = T # is class column categorical?
interested_cols = c("sex","genotype_centre","genotype_batch", "allergen", 
                    "race","response",
                    "age","bmi","weight", "height") 
interested_cont_cols = c("age","bmi","weight", "height")

# cid_col = "id"

dofullPCA = F #do PCA for all cell popoulations not just first level
pca_k = 5 #number of pca pc to plot

doISO = T #do ISO feature reduction
iso_k = 3 # number of neighbours
mds_type = c("iso", "mds")

doTsne = T #do Tsne feature reduction
theta=.5 #parameter for Tsne

doHC = F #do hierarchical clustering

#plot
wdth = 400
ht = 600












## features and indices
feat_types = list.files(feat_dir,full.names=F,pattern=".Rdata")
feat_types = feat_types[!grepl("raw",feat_types)]
feat_types = gsub(".Rdata","",feat_types)
feat_temp = str_split(feat_types,"[.]")
feat_names = sapply(feat_temp, function(x) x[1])
feat_times = sapply(feat_temp, function(x) x[2])

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




start = Sys.time()


meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
# meta_col0 = as.data.frame(get(load(paste0(meta_col_dir,".Rdata"))))

for (feat_type in feat_types) {
  cat("\n", feat_type)
  
  m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
  if (sum(colnames(m0)%in%meta_file0[,id_col])==ncol(m0)) m0 = t(m0)
  
  col_inds_i = sapply(names(col_inds0), function(x) grepl(x,feat_type))
  if (sum(col_inds_i)==0) {
    col_inds = list(all=c(""))
  } else {
    col_inds = col_inds0[[col_inds_i]]
  }
  
  
  for (file_ind_n in names(file_inds)) {
    file_ind = file_inds[[file_ind_n]]
    if (file_ind_n=="all") file_ind = rownames(m0)
    # file_ind_n = paste0("-",file_ind_n)
    
    for (col_ind_n in names(col_inds)) {
      col_ind = col_inds[[col_ind_n]]
      if (col_ind_n=="all") col_ind = colnames(m0)
      # col_ind_n = paste0(".",col_ind_n)

      m = m0[rownames(m0)%in%file_ind, colnames(m0)%in%col_ind]
      m = delna(m)
      if (any(dim(m)==0) | sum(!is.na(m))<sum(is.na(m))) next()
      m_cont_col = apply(m,2,function(x) length(unique(x))>cont_col)
      if (scale_cont) m[,m_cont_col] = scale(as.numeric(as.matrix(m[,m_cont_col])))
      m[is.na(m)] = -1
      
      cat(" (",file_ind_n, " x ",col_ind_n,"; ",nrow(m)," x ",ncol(m),") ")
      
      
      meta_file = meta_file0[match(rownames(m),meta_file0[,id_col]),]
      # meta_col = meta_col0[match(colnames(m),meta_col0[,cid_col]),]
      
      
      ## get interested columns
      meta_file_factor = as.data.frame(sapply(meta_file, function(x) as.numeric(factor(x, ordered=T))))
      
      interested_col_ind = which(colnames(meta_file)%in%interested_cols)
      uniquecols = apply(meta_file_factor, 2, function(x) nrow(meta_file_factor)==length(unique(x)))
      meta_file_attonly = meta_file_factor[,!uniquecols]
      
      
      
      ## pca analysis ------------------------------------------------
      cat("pca; ")
      pc <- prcomp(as.matrix(m))
      
      #pca scatterplot
      pngname = paste0(pca_dir,"/",feat_type,"-",file_ind_n,"X",col_ind_n,"_pca-iso.png")
      png(filename=pngname, width=length(interested_col_ind)*wdth, height=(1+doISO+pca_k)*ht)
      layout(matrix(c(rep(1,length(interested_col_ind)), 2:(((2*pca_k)+doISO+1)*length(interested_col_ind)+1)),ncol=length(interested_col_ind),byrow=T))
      par(cex=1)
      
      plot(pc$sdev^2/sum(pc$sdev^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
      for (i in 1:(pca_k-1)) {
        for (col in interested_col_ind) {
          colname = colnames(meta_file)[col]
          coloursm = meta_file_factor[,col]
          plot(pc$x[,i+1], pc$x[,i], col = coloursm, main = paste0("PCA ", colname), xlab = paste0("PC_",i+1), ylab = paste0("PC_",i))
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
            boxplot(pc$x[,i]~, lwd = 1, outline=F, ylim=c(min(pc$x[,i]),max(pc$x[,i])),
                    main = paste0("PCA ", colname," Pearson Cor (ok if 2 var values) = ", round(cor,3)),
                    xaxt = 'n', xlab=colname, ylab = paste0("PC_",i)) #,xaxt ='n', xlab=testcol
            axis(1, at=1:length(attributenames), labels=attributenames, las=2)
            points(jitter(attribute, factor=1), pc$x[,i], col = coloursm)
          }
        }
      }
      if (doISO) {
        fit <- Isomap(as.matrix(m),k=iso_k,verbose=F,plotResiduals=T)
        for (col in interested_col_ind) {
          coloursm <- meta_file_factor[,col]
          plot(fit$dim2, col = coloursm, main = paste0("ISO ", colname), las=2)
          points(0, 0, pch = 3, cex = 4, lwd = 4)
        }
      }
      graphics.off()
      
    } #col_ind
  } #row_ind
} #feat_type

time_output(start)



