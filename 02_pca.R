## input: features & meta_file
## output: plot PCA plot to fin dconfounders
## aya43@sfu.ca
## created 20180524
## last modified 20180524

## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/src/_dirs.R"))
source(paste0(root, "/src/_func.R"))
source(paste0(root, "/src/_func-asthma.R"))
source(paste0(root, "/src/_func-classifiers.R"))
libr(append(pkgs(),c("RDRToolbox")))

# no_cores = 1#detectCores()-2 #number of cores to use for parallel processing
# registerDoMC(no_cores)


## options
writecsv = T #write results as csv on top of Rdata?

cont_col = 15 #if a column has more than this unique values, it's continuous
scale_cont = T #if a column is continuous, scale it


interested_cols = c("sex","dna_centre","dna_batch", "allergen", 
                    "race","response",
                    "age","bmi","weight", "height") 
interested_cont_cols = c("age","bmi","weight", "height")

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
feat_types = feat_types_annot


## start
start = Sys.time()

for (feat_type in feat_types) {
  cat("\n", feat_type)
  
  # load m0 using: meta_file0, feat_dir, feat_type, col_inds0 --> m0, col_inds
  source(paste0(root, "/src/_func-asthma_m0-load.R")) 
  
  class_coli = 1
  for (file_ind_n in names(file_inds)) {
    for (col_ind_n in names(col_inds)) {
      
      # prep m, meta_file, meta_col: m0 class_coli, class_cols, file_ind_n, file_inds
      source(paste0(root, "/src/_func-asthma_m-trim.R")); if (nextm) next()
      
      
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
            plot(attribute, pc$x[,i], main=paste0("PCA ", colname," Pearson Cor = ", cor), 
                 xlab = colname, ylab = paste0("PC_",i)) #, col=coloursm
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



