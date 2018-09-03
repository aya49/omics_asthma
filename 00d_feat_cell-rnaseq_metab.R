## input: rnaseq data from rnaseq_analysis-meta
## output: well formatted data (cell counts inferred from rnaseq data) and meta col files
## aya43@sfu.ca
## created 20180720


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/src/_dirs.R"))
source(paste0(root, "/src/_func.R"))
source(paste0(root, "/src/_func-asthma.R"))
libr(append(pkgs(),c("mixOmics","CellCODE","biomaRt")))


## import clinical datasets
demo = readRDS(meta_file_rnaseqa_dir)
rownames(demo) = paste(demo$concealedID, demo$Time, sep=".")

## load RNA-Seq datasets
load(rnaseqa_data_dir)
load(rnaseqa_datanorm_dir)


## import clinical datasets
demo = readRDS(meta_file_rnaseqa_dir)
rownames(demo) = paste(demo$concealedID, demo$Time, sep=".")



## METAB ######################################

# Cohort 1 (14 subjects)
# these variables aren't used...

asthma_templt = read.csv(asthma_templt_dir, header=F, stringsAsFactors=F)

cohort1lod = as.numeric(as.matrix(asthma_templt[5, -c(1:23)]))
names(cohort1lod) = as.character(as.matrix(asthma_templt[2, -c(1:23)]))

cohort1Dat = t(sapply(6:nrow(asthma_templt), function(i) as.numeric(as.matrix(asthma_templt[i, -c(1:23)])) ))
rownames(cohort1Dat) = paste(asthma_templt[-c(1:5), 4], "cohort1", sep = " ")
colnames(cohort1Dat) = as.character(as.matrix(asthma_templt[2, -c(1:23)]))
all(names(cohort1lod)==colnames(cohort1Dat))

cohort1pathway = factor(as.character(as.matrix(asthma_templt[3, -c(1:23)])))
names(cohort1pathway) = as.character(as.matrix(asthma_templt[2, -c(1:23)]))
all(names(cohort1lod)==colnames(cohort1pathway))


# Cohort 2 (14 subjects)
# these variables aren't used...

## import dataset
asthma_discov = read.csv(asthma_discov_dir, header=F)

cohort2lod = as.numeric(as.matrix(asthma_discov[16, -c(1:22)]))
names(cohort2lod) = as.character(as.matrix(asthma_discov[2, -c(1:22)]))

cohort2Dat = t(sapply(17:nrow(asthma_discov), function(i) as.numeric(as.matrix(asthma_discov[i, -c(1:22)])) ))
rownames(cohort2Dat) = paste(asthma_discov[-c(1:16), 4], "cohort2", sep = " ")
colnames(cohort2Dat) = as.character(as.matrix(asthma_discov[2, -c(1:22)]))
all(names(cohort2lod)==colnames(cohort2Dat))

cohort2pathway = factor(as.character(as.matrix(asthma_discov[4, -c(1:22)])))
names(cohort2pathway) = as.character(as.matrix(asthma_discov[2, -c(1:22)]))
all(names(cohort2lod)==colnames(cohort2pathway))

## remove quality controls, standards, pbs and blank
cohort2Dat.new = cohort2Dat[
  -unlist(lapply(c("Blank", "KIT1", "PBS", "Standard"),
                 function(i) agrep(i, rownames(cohort2Dat)) )), ]


# Combine both plate and extract pre-challenge samples

cohortDat = rbind(cohort1Dat[c("L_ST_003  Pre cohort1", "L_ST_018  Pre cohort1",
                               "L_ST_026  Pre cohort1", "L_ST_018  Pre (repeat) cohort1",
                               "L_ST_026  Pre (repeat) cohort1"), ],
                  cohort2Dat.new)
rownames(cohortDat) = gsub("^ *|(?<= ) | *$", "", rownames(cohortDat), perl = T)

subjNames = unlist(lapply(strsplit(rownames(cohortDat), " "), function(i) i[1]))
time = unlist(lapply(strsplit(rownames(cohortDat), " "), function(i) i[2]))
time[time == "0h"] = "Pre"
time[time == "3h"] = "Post"

## filter metabolites with mean expression less than the LOD; scale
cohortDatFiltered = cohortDat#[, colMeans(cohortDat) > cohort2lod]
cohortDatFilterednorm = scale(cohortDatFiltered)

col = paste(subjNames, time)
col[col %in% names(table(col)[table(col) == 1])] = "results/enrichr" # unique name time
col2 = as.numeric(factor(col))

## heatmap
pdf(metab_hist_dir, width = 9, height = 9)
cim(cor(t(cohortDatFilterednorm), method = "spearman"),
    # row.names = F,
    # col.names = F,
    row.sideColors = color.mixo(col2),
    col.sideColors = color.mixo(col2),
    margins = c(5, 15))
legend(x = 1.5, y = 62, levels(factor(col)),
       col = color.mixo(1:6),
       pch = 19, bty = "n")
dev.off()

## Take average of replicates (transposed)
expDat0 = as.data.frame(cohortDatFilterednorm)
expDat0$Rep = paste(subjNames, time)

# expDat = expDat0 %>% group_by(Rep) %>% do(x=colMeans(.[,-ncol(.)]))
# exp = do.call(cbind, expDat$x)
# colnames(exp) = expDat$Rep

expDat1 = expDat0[,!colnames(expDat0)%in%"Rep"]
reptemp = unique(gsub(" cohort[0-9]","",rownames(expDat0)))
exp = t(Reduce("rbind",lapply(reptemp, function(x) 
  colMeans(expDat1[expDat0$Rep==x,,drop=F]))))
colnames(exp) = reptemp


hcExp0 = exp[, agrep("HLI", colnames(exp))]
colnames(hcExp0) = c("HLI-Normal 1 3h (12PM)", "HLI-Normal 1 0h (9AM)",
                     "HLI-Normal 2 3h (12PM)", "HLI-Normal 2 0h (9AM)",
                     "HLI-Normal 3 3h (12PM)", "HLI-Normal 3 0h (9AM)",
                     "HLI-Normal 4 3h (12PM)", "HLI-Normal 4 0h (9AM)",
                     "HLI-Normal 5 3h (12PM)", "HLI-Normal 5 0h (9AM)",
                     "HLI-Normal 6 3h (12PM)", "HLI-Normal 6 0h (9AM)")
hcDemo = demo[!is.na(match(as.character(demo$UniqueID), colnames(hcExp0))), ]
hcExp = hcExp0[, as.character(hcDemo$UniqueID)]
all(colnames(hcExp) == as.character(hcDemo$UniqueID))

asthmaDemo = demo[paste(as.character(demo$MST_LST_numbers), as.character(demo$Time), sep = " ") %in% gsub(" cohort[0-9]+","",colnames(exp)), ]
asthmaExp = exp[, paste(as.character(asthmaDemo$MST_LST_numbers), as.character(asthmaDemo$Time), sep = " ")]
colnames(asthmaExp) = rownames(asthmaDemo)

# save
metab = t(asthmaExp)
metab_pre = metab[grepl("Pre",asthmaDemo$Time),]
rownames(metab_pre) = as.character(asthmaDemo$NAME[grepl("Pre",asthmaDemo$Time)])
metab_pre = rbind(colMeans(metab_pre[rownames(metab_pre)%in%"L_ST_026",]),
                  metab_pre[!rownames(metab_pre)%in%"L_ST_026",])
rownames(metab_pre)[1] = "L_ST_026"
metab_post = metab[grepl("Post",asthmaDemo$Time),]
rownames(metab_post) = as.character(asthmaDemo$NAME[grepl("Post",asthmaDemo$Time)])
metab_post = rbind(colMeans(metab_post[rownames(metab_post)%in%"L_ST_026",]),
                   metab_post[!rownames(metab_post)%in%"L_ST_026",])
rownames(metab_post)[1] = "L_ST_026"

rownames(metab_pre)[grepl("WRF",rownames(metab_pre))] =
  rownames(metab_post)[grepl("WRF",rownames(metab_post))] = "WRF"

metab_diff = metab_post-metab_pre

save(metab_pre, file=paste0(feat_metab_dir,".pre.Rdata"))
write.csv(metab_pre, file=paste0(feat_metab_dir,".pre.csv"))
save(metab_post, file=paste0(feat_metab_dir,".post.Rdata"))
write.csv(metab_post, file=paste0(feat_metab_dir,".post.csv"))
save(metab_diff, file=paste0(feat_metab_dir,".diff.Rdata"))
write.csv(metab_diff, file=paste0(feat_metab_dir,".diff.csv"))



## pre&post-challenge data RNASEQ #################
for (time in c("Pre","Post")) {
  
  asthmaDemoPre = asthmaDemo[asthmaDemo$Time == time, ]
  asthmaExpPre = t(asthmaExp[, rownames(asthmaDemoPre)])
  
  ## use ensembl datasets
  # remove ERCC controls
  ensemblDat = normalizelibSum(starEnsemblExp[-(60156:nrow(starEnsemblExp)), intersect(rownames(asthmaExpPre), colnames(genDats$starEnsemblExp))])
  dim(ensemblDat)
  starDemo = asthmaDemo[colnames(ensemblDat), ]
  
  ## get markers genes: DMAP, iris, lm22, pancancer genes #############################
  data(DMAP) #cellCODE
  dmapTag=tagData(DMAP, 0.5, max=15, ref=NULL, ref.mean=F)
  data(IRIS) #cellCODE
  irisTag=tagData(IRIS, 1, max=15, ref=NULL, ref.mean=F);
  lm22 = read.delim(cell_gene2_dir, row.names = 1) ####MISSING
  lm22Tag=tagData(lm22, 0.5, max=15, ref=NULL, ref.mean=F)
  pancancer0 = read.csv(cell_pc_immune_dir)
  dim(dmapTag); dim(irisTag); dim(lm22Tag);
  
  dmapList = apply(dmapTag, 2, function(i) names(i)[which(i != 0)])
  irisList = apply(irisTag, 2, function(i) names(i)[which(i != 0)])
  lm22List = apply(lm22Tag, 2, function(i) names(i)[which(i != 0)])
  pancancerList = split(as.character(pancancer0[pancancer0$Cell.Type != "", c("HUGO.Name")]), as.character(pancancer0[pancancer0$Cell.Type != "", c("Cell.Type")]))
  
  dmapDat = data.frame(Gene = unlist(dmapList), Cell.Type = paste("dmap",rep(names(dmapList), unlist(lapply(dmapList, length))), sep = "_"))
  irisDat = data.frame(Gene = unlist(irisList), Cell.Type = paste("iris",rep(names(irisList), unlist(lapply(irisList, length))), sep = "_"))
  lm22Dat = data.frame(Gene = unlist(lm22List), Cell.Type = paste("lm22",rep(names(lm22List), unlist(lapply(lm22List, length))), sep = "_"))
  pancancerDat = data.frame(Gene = unlist(pancancerList), Cell.Type = paste("pancancer",rep(names(pancancerList), unlist(lapply(pancancerList, length))), sep = "_"))
  
  cellSpecificList0 = rbind(dmapDat, irisDat, pancancerDat, lm22Dat)
  cellSpecificList = split(cellSpecificList0$Gene, cellSpecificList0$Cell.Type)
  
  save(cellSpecificList, file=cell_gene_dir)
  
  ids = unique(unlist(cellSpecificList))
  Dat = do.call(cbind, lapply(cellSpecificList, function(i){
    x = rep(0, length(ids))
    x[ids %in% i] = 1
    x
  }))
  rownames(Dat) = ids
  
  ## annoate gene symbols to ensembl ids
  
  ensemIDs = 
    getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
          filters = 'hgnc_symbol',
          values = ids,
          mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                         host="grch37.ensembl.org", 
                         path="/biomart/martservice", 
                         dataset="hsapiens_gene_ensembl"))
  
  cc = unique(ensemIDs$hgnc_symbol)
  Dat2 = Dat[cc, ]
  rowNames = unlist(lapply(strsplit(rownames(ensemblDat), "\\."), function(i) i[1]))
  ccDat = unlist(lapply(cc, function(x){
    ensemblID = rownames(ensemblDat)[rowNames %in% subset(ensemIDs, hgnc_symbol == x)$ensembl_gene_id]
    if(length(ensemblID) == 0) ensemblID = NA
    ensemblID
  }))
  Dat3 = Dat2[!is.na(ccDat), ]
  rownames(Dat3) = ccDat[!is.na(ccDat)]
  Dat4 = Dat3[, colSums(Dat3) > 3]
  
  # correlation btwn ENSG x samples and ENSG x marker genes affected cell types (binary matrix)
  SPVs=getAllSPVs(ensemblDat, grp=starDemo$calculated_Response, Dat4, method="mixed", plot=F, mix.par = 0.3)
  rownames(SPVs) = rownames(starDemo)
  
  pdf(paste0(cell_type_fig_dir,time,".pdf"), width = 9, height = 9)
  cim(cor(SPVs), margins = c(15, 15))
  dev.off()
  
  rownames(SPVs) = as.character(starDemo$NAME)
  rownames(SPVs)[grepl("WRF",rownames(SPVs))] = "WRF"
  
  SPVs = exp(SPVs)
  
  save(SPVs, file=paste0(feat_cellseqgenes_dir,".",tolower(time),".Rdata"))
  write.csv(SPVs, file=paste0(feat_cellseqgenes_dir,".",tolower(time),".csv"))
}
a = get(load(paste0(feat_cellseqgenes_dir,".post.Rdata"))) - 
  get(load(paste0(feat_cellseqgenes_dir,".pre.Rdata")))
save(a, file=paste0(feat_cellseqgenes_dir,".diff.Rdata"))
write.csv(a, file=paste0(feat_cellseqgenes_dir,".diff.csv"))

#genes that went into making cells
celldiff = names(rowSums(Dat4)>0)
save(celldiff, file=paste0(meta_col_dir,"rnaseqstarensgenes_id_cell.Rdata"))
