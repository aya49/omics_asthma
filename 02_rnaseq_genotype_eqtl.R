## input: features & meta_file
## output: EQTL
## aya43@sfu.ca
## created 20180614



## root directory
root = "~/projects/asthma"
setwd(root)

result0_dir = paste0(root, "/result")
result_dir = list.dirs(result0_dir, recursive=F)


# asthma = "asthma" # "asthma" if only test asthma related SNP; else ""


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_genotype_dir = paste0(feat_dir,"/snp-file-genotype")


## output directory
stat_dir = paste0(result0_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
eqtl_cis_dir = paste0(stat_dir,"/eqtlcis")



## libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affy","derfinder"))
source("code/_func.R")
libr("rtracklayer") #ucsc
libr("TxDb.Hsapiens.UCSC.hg19.knownGene")
libr("org.Hs.eg.db") #org.* annotation packages; can forge own and interact with using library("AnnotationDbi")
libr("data.table")
libr("MatrixEQTL")
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
cid_col = "id"
class_col = "response"
categorical = T # is class column categorical?
interested_cols = c("age","bmi","sex","centre","batch","race","response") 
interested_cont_cols = ""



# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");

## Location of the package with the data files.
# base.dir = find.package('MatrixEQTL');
# base.dir = '.';




gt_inds = grep("genotyp",result_dir)
rs_inds = grep("RNAseq",result_dir)

for (gt_ind in gt_inds) {
  for (rs_ind in rs_inds) {
    
    # load data
    gt_meta_file0 = get(load(paste0(meta_file_dir[gt_ind],".Rdata")))
    gt_meta_col0 = get(load(paste0(meta_col_dir[gt_ind],".Rdata")))
    gt_m0_names = list.files(feat_dir[gt_ind],full.names=F)
    gt_m0_name = gt_m0_names[which.min(nchar(gt_m0_names))]
    gt_m0 = as.matrix(get(load(paste0(feat_dir[gt_ind],"/",gt_m0_name))))
    if (colnames(gt_m0)[1]%in%gt_meta_file0[,id_col]) gt_m0 = t(gt_m0)
    
    rs_meta_file0 = get(load(paste0(meta_file_dir[rs_ind],".Rdata")))
    rs_meta_col0 = get(load(paste0(meta_col_dir[rs_ind],".Rdata")))
    rs_m0_names = list.files(feat_dir[rs_ind],full.names=F)
    rs_m0_name = rs_m0_names[which.min(nchar(rs_m0_names))]
    rs_m0 = as.matrix(get(load(paste0(feat_dir[rs_ind],"/",rs_m0_name))))
    if (colnames(rs_m0)[1]%in%rs_meta_file0[,id_col]) rs_m0 = t(rs_m0)

    # prepare data as SlicedData for input into MatrixEQTL()
    gt_sd = SlicedData$new()
    gt_sd$CreateFromMatrix(t(gt_m0))
    
    rs_sd = SlicedData$new()
    rs_sd$CreateFromMatrix(t(rs_m0))
    # # Slice it in pieces of 2 rows
    # sd$ResliceCombined(sliceSize = 2L)  
    # # Show the number of slices (equivalent function calls)
    # sd$nSlices()
    # length(sd)
    # # Is it all in one slice? (No)
    # sd$IsCombined()
    # # Show the column names (equivalent function calls)
    # sd$columnNames
    # colnames(sd)
    # # Show all row names (equivalent function calls)
    # sd$GetAllRowNames()
    # rownames(sd)
    # # Print the second slice
    # print(sd[[2]])
    # # Reorder and subset columns
    # sd$ColumnSubsample( c(1,3,4) )  
    # # Reorder and subset rows
    # sd$RowReorder( c(3,1) )  
    # # Find the row with name "row1" (it is second in the first slice)
    # sd$FindRow("row1")
    # # Show the detail of the object (one slice again)
    # show(sd)
    
    # prepare SNP and gene positions
    gt_pos = data.frame(snpid=gt_meta_col0$dbSNP, chr=paste0("chr",gt_meta_col0$chromosome), pos=gt_meta_col0$pos_phys)
    # levels(gt_pos$chr) = paste("chr", c(1:22, "X", "Y", "M"), sep="") #convert to bioconductor format
    rs_pos = data.frame(geneid=, chr=, left=, right=)
    cgene_col = ifelse("gene"%in%colnames(rs_meta_col0),"gene","id")
    
    ## gene look up
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene #transcriptDb; behind the scenes, everythign is SQLite
    # group individual ranges into groups: transcriptsBy, cdsBy, exonsBy
    tx.by.gene = transcriptsBy(txdb, "gene") #list names are Entrez gene ID's
    # gene names; Entrez Gene ID
    # tx.by.gene # GRangesList object
    # columns(org.Hs.eg.db)
    # keys: APOE gene
    select(org.Hs.eg.db, keys="APOE", 
           columns=c("ENTREZID", "SYMBOL", "GENENAME"), 
           keytype="SYMBOL") #keytypes()
    # look up Gene ID
    tx.by.gene["348"]
    
  }
}




# Gene expression file name
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");

# # Output file name
# output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvalthres_cis = 2e-2;
pvalthres_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;



## Run the analysis
# snpspos = read.table(snps_location_file_name, header=T, stringsAsFactors=T);
# genepos = read.table(gene_location_file_name, header=T, stringsAsFactors=T);

me = Matrix_eQTL_main(
  snps = gt_sd, 
  gene = rs_sd, 
  # cvrt = cvrt, # covariate SlicedData
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvalthres_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = paste0(eqtl_cis_dir,"_",gsub(".Rdata","",gt_m0_name),".",gsub(".Rdata","",rs_m0_name),".Rdata"), #output_file_name_cis=tempfile()
  pvOutputThreshold.cis = pvalthres_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the histogram of local and distant p-values

plot(me)












