## input: meta_col & indices
## output: string function-ed gene etc maps with snps (pretty-fied info about snps)
## aya43@sfu.ca
## created 20180524



## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
libr(append(pkgs(),c("qqman","biomaRt")))



## options
options(stringsAsFactors=F)

no_cores = 5#detectCores()-3
registerDoMC(no_cores)

overwrite = F
writecsv = T #write results as csv on top of Rdata?

pthres = .05
padjust = p.adjust.methods


cid_col = "probe"
cinclud_col = c("dbSNP","pos_phys","chromosome","Associated Gene","EBI MAPPED GENE(S)","ClinVar OMIM Gene", "ClinVar OMIM Description")
key_word = "asthma"


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
meta_col0 = as.data.frame(get(load(paste0(meta_col_gt_dir,".Rdata"))))



key_row = apply(meta_col0,1,function(x) any(grepl(key_word,x,ignore.case=T))) & !duplicated(meta_col0$dbSNP)
key_col = apply(meta_col0,2,function(x) any(grepl(key_word,x,ignore.case=T)))

meta_col = meta_col0[key_row,c(cinclud_col,colnames(meta_col0)[key_col])]
meta_col = meta_col[,apply(meta_col,2,function(x) sum(!duplicated(x))>1)]


associated_gene0 = meta_col[,"Associated Gene"]
associated_gene = str_split(associated_gene0,"///") #separate genes
names(associated_gene) = meta_col$dbSNP
associated_gene1 = lapply(associated_gene, function(ag) {
  genes = (str_split(ag,"//")) #separate gene attributes
  gene_names = sapply(genes,function(x) x[5]) #get gene names
  genes = genes[order(gene_names)] #sort by gene names
  genes = lapply(genes, trimws) #trim leading and trailing whitespace
  
  #merge genes with the same name; different transcript_accession
  new_genes = list()
  dupgene = duplicated(gene_names)
  newgenei = 0
  dupgenei = 0
  while (dupgenei<length(gene_names)) {
    dupgenei = dupgenei+1
    if (!dupgene[dupgenei]) {
      newgenei = newgenei+1
      new_genes[[newgenei]] = genes[[dupgenei]]
    } else {
      new_genes[[newgenei]][1] = 
        paste(new_genes[[newgenei]][1], 
              genes[[dupgenei]][1], collapse=", ")
    }
  }
  return(new_genes)
})
associated_gene1_len = sapply(associated_gene1, length)

meta_col_snp = as.data.frame( cbind(rep(names(associated_gene1_len),associated_gene1_len),
                                    Reduce("rbind",unlist(associated_gene1,recursive=F))) )

colnames(meta_col_snp) = c("dbSNP","transcript_accession", "SNP-gene_relation", 
                           "dist_0-if-in-gene", "unigene_cluster", "gene_name", 
                           "NCBI_id", "GenBank_descr")

meta_col2 = as.data.frame( 
  cbind(meta_col_snp[,c("transcript_accession", "SNP-gene_relation", 
                       "dist_0-if-in-gene", "gene_name","GenBank_descr")],
        meta_col[match(as.vector(meta_col_snp[,"dbSNP"]),meta_col$dbSNP),]) )
meta_col2 = meta_col2[order(meta_col2$chromosome,meta_col2$pos_phys),!colnames(meta_col2)%in%"Associated Gene"]
meta_col2 = meta_col2[,c("dbSNP","chromosome","pos_phys","gene_name",
                         "SNP-gene_relation","dist_0-if-in-gene","GenBank_descr","transcript_accession",
                         "EBI MAPPED GENE(S)","EBI DISEASE/TRAIT","EBI MAPPED_TRAIT")]

write.csv(meta_col2,file=meta_col_pretty_dir, row.names=F)
time_output(start)


