## input: features & meta_file
## output: plot PCA plot to fin dconfounders
## aya43@sfu.ca
## created 20180524
## last modified 20180524



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result/genotype") #"genotype", "RNAseq_genes", "RNAseq_isoforms"


# asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)


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

id_col = "fileName"
class_col = "response"
categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""

# cid_col = "id"

dofullPCA = F #do PCA for all cell popoulations not just first level
pca_k = 5 #number of pca pc to plot

doISO = T #do ISO feature reduction
iso_k = 3 # number of neighbours
mds_type = c("iso", "mds")

doTsne = T #do Tsne feature reduction
theta=.5 #parameter for Tsne

doHC = F #do hierarchical clustering













## features and indices
feat_types = list.files(feat_dir,full.names=F,pattern=".Rdata")
feat_types = gsub(".Rdata","",feat_types)

col_inds_paths = list.files(meta_dir,pattern="col_id_",full.names=T)
col_inds_names = sapply(str_split(gsub(".Rdata","",col_inds_paths),"_"), function(x) x[length(x)])
col_inds = lapply(col_inds_paths, function(x) get(load(x)))
names(col_inds) = col_inds_names
col_inds = append(list(all=""), col_inds)

file_inds_paths = list.files(meta_dir,pattern="file_id_",full.names=T)
file_inds_names = sapply(str_split(gsub(".Rdata","",file_inds_paths),"_"), function(x) x[length(x)])
file_inds = lapply(file_inds_paths, function(x) get(load(x)))
names(file_inds) = file_inds_names
file_inds = append(list(all=""), file_inds)




start = Sys.time()


meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
# meta_col0 = as.data.frame(get(load(paste0(meta_col_dir,".Rdata"))))

for (feat_type in feat_types) {
  m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
  if (sum(colnames(m0)%in%meta_file0[,id_col])==ncol(m0)) m0 = t(m0)
  
  for (file_ind_n in names(file_inds)) {
    file_ind = file_inds[[file_ind_n]]
    if (file_ind_n=="all") file_ind = rownames(m0)
    # file_ind_n = paste0("-",file_ind_n)
    
    for (col_ind_n in names(col_inds)) {
      col_ind = col_inds[[col_ind_n]]
      if (col_ind_n=="all") col_ind = colnames(m0)
      # col_ind_n = paste0(".",col_ind_n)
      
      m = m0[file_ind[file_ind%in%rownames(m0)], col_ind[col_ind%in%colnames(m0)]]
      mna = apply(m,2,function(x) sum(is.na(x)))
      print(table(mna[mna>0]))
      m = m[,mna<=good_col]
      m[is.na(m)] = -1
      
      meta_file = meta_file0[match(rownames(m),meta_file0[,id_col]),]
      # meta_col = meta_col0[match(colnames(m),meta_col0[,cid_col]),]
      
      
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
      pc <- prcomp(m)
      
      #pca scatterplot
      pngname = paste0(stat_dir,"/",feat_type,"-",file_ind_n,"X",col_ind_n,"_pca-iso.png")
      png(filename=pngname, width=length(interested_col_ind)*400, height=(1+doISO+pca_k)*400)
      layout(matrix(c(rep(1,length(interested_col_ind)), 2:(((2*pca_k)+doISO)*length(interested_col_ind)+1)),ncol=length(interested_col_ind),byrow=T))
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
      
    } #col_ind
  } #row_ind
} #feat_type

time_output(start)



