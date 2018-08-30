## input: meta + rnaseq
## output: diablo
## aya43@sfu.ca
## created 20180720



## root directory
root = "~/projects/asthma"
setwd(root)

result_dir = paste0(root, "/result")


## input directory
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = passte0(meta_dir,"/file")
meta_col_dir = paste0(meta_dir,"/col")

feat_dir = paste0(result_dir,"/feat")
# feat_dna_dir = paste0(feat_dir,"/snp-file-dna")


## output directory
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
eqtl_dir = paste0(stat_dir,"/eqtl"); dir.create(eqtl_dir, showWarnings=F)
gwas_dir = paste0(stat_dir,"/gwas"); dir.create(gwas_dir,showWarnings=F)


















## import clinical datasets
demo = readRDS("~/projects/asthma/data/RNAelements/data/demo/allsitesDemo.rds")
rownames(demo) = paste(demo$concealedID, demo$Time, sep=".")

## load RNA-Seq datasets
load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_rawData.RDATA")
load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_normalized.RDATA")

asthma_templt_dir   = paste0(root, "/data/RNAelements/data/asthmaTestplate.Jan16_2015.csv")
asthma_discov_dir = paste0(root, "/data/RNAelements/data/2015-02-10_Conc_asthmaDiscoveryPlateFeb10.2015.csv")

## import clinical datasets
demo = readRDS("~/projects/asthma/data/RNAelements/data/demo/allsitesDemo.rds")
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

cohort2lod = as.numeric(as.matrix(asthma_discov[16, -c(1:23)]))
names(cohort2lod) = as.character(as.matrix(asthma_discov[2, -c(1:23)]))

cohort2Dat = t(sapply(17:nrow(asthma_discov), function(i) as.numeric(as.matrix(asthma_discov[i, -c(1:23)])) ))
rownames(cohort2Dat) = paste(asthma_discov[-c(1:16), 4], "cohort2", sep = " ")
colnames(cohort2Dat) = as.character(as.matrix(asthma_discov[2, -c(1:23)]))
all(names(cohort2lod)==colnames(cohort2Dat))

cohort2pathway = factor(as.character(as.matrix(asthma_discov[4, -c(1:23)])))
names(cohort2pathway) = as.character(as.matrix(asthma_discov[2, -c(1:23)]))
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
pdf(paste0(stat_dir, "/hist_metab.pdf"), width = 9, height = 9)
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
expDat = expDat0 %>% group_by(Rep) %>% do(x=colMeans(.[,-ncol(.)]))
exp = do.call(cbind, expDat$x)
colnames(exp) = expDat$Rep

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

asthmaDemo = demo[paste(as.character(demo$MST_LST_numbers), as.character(demo$Time), sep = " ") %in% colnames(exp), ]
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

save(metab_pre, file=paste0(feat_dir,"/metab.pre.Rdata"))
write.csv(metab_pre, file=paste0(feat_dir,"/metab.pre.csv"))
save(metab_post, file=paste0(feat_dir,"/metab.post.Rdata"))
write.csv(metab_post, file=paste0(feat_dir,"/metab.post.csv"))
save(metab_diff, file=paste0(feat_dir,"/metab.diff.Rdata"))
write.csv(metab_diff, file=paste0(feat_dir,"/metab.diff.csv"))



## pre-challenge data RNASEQ #################
asthmaDemoPre = asthmaDemo[asthmaDemo$Time == "Pre", ]
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
lm22 = read.delim(paste0(root, "/data/LM22.txt"), row.names = 1) ####MISSING
lm22Tag=tagData(lm22, 0.5, max=15, ref=NULL, ref.mean=F)
pancancer0 = read.csv(paste0(root, "/data/Cells_nCounter_Human_PanCancer_Immune_Profiling_Panel_Gene_List.csv"))
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

save(cellSpecificList, file=paste0(result_dir,"/cellgene_dmap.iris.lm22.pc.Rdata"))

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

pdf(paste0(root, "/data/RNAelements/figs/cellTypes.pdf"), width = 9, height = 9)
cim(cor(SPVs), margins = c(15, 15))
dev.off()

## pre-challenge
Y = factor(starDemo$calculated_Response, levels = c("ER", "DR"))
names(Y) = rownames(starDemo)
X = list(Cells = SPVs, Genes = t(genDats$starEnsemblExp[is.na(match(rownames(genDats$starEnsemblExp), rownames(Dat4))), names(Y)]), 
         Metabolites = asthmaExpPre[rownames(starDemo), ])
dim(X[[1]]); dim(X[[2]]); dim(X[[3]]); 

#ncomp = rep(2, length(X))
ncomp = 2
design = matrix(c(0, 0, 0,
                  0, 0, 0,
                  0, 0, 0), nrow = 3, ncol = 3)
keepX = list(Cells = rep(2, ncomp[1]), Genes = rep(10, ncomp[1]), Metabolites = rep(5, ncomp[1]))
## block.splsda = horizontal integration PLS-DA model with a specified number of components per block either by Y or by its position indY in the list of blocks X
# result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
#                       keepX = keepX, design = design,
#                       mode = "regression", bias = T)
result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                      keepX = keepX, design = design,
                      mode = "regression") # PLS regression ("regression"), PLS canonical analysis ("canonical"), redundancy analysis ("invariant") and the classical PLS algorithm ("classic")
# get variables that have >0 weights
feat1 = lapply(result$loadings, function(x) apply(x, 2, function(i) names(i)[which(i != 0)]))

pdf(paste0(root, "/data/RNAelements/figs/SamplePlot_rnaseq_dmap-iris-lm22-pc_v_celltype.pdf"), width = 9, height = 5)
# plotDiablo3(result, ncomp = 1, groupOrder = c("DR", "ER"))
plotDiablo(result, ncomp = 1, groupOrder = c("DR", "ER"))
dev.off()

pdf(paste0(root, "/data/RNAelements/figs/circosPlot_rnaseq_dmap-iris-lm22-pc_v_celltype.pdf"), width = 7)
# circosPlot_diabloModif(result, corThreshold = 0.5, cex.label = 0.5, showIntraLinks = F)
circosPlot(result, showIntraLinks = F, cutoff=.5)#, corThreshold = 0.5, cex.label = 0.5)
dev.off()

pdf(paste0(root, "/data/RNAelements/figs/heatmap_rnaseq_dmap-iris-lm22-pc_v_celltype.pdf"), width = 7)
# heatmap_diablo(result, margins = c(2, 12))
cimDiablo(result, margins = c(2, 12))
dev.off()


## get genes associated with splsda (sparse pls-da) loadings
hk.known = 
  getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', "description"),
        filters = 'ensembl_gene_id',
        values = unlist(lapply(strsplit(as.character(feat1$Genes), "\\."), function(i) i[1])),
        mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                       host="grch37.ensembl.org", 
                       path="/biomart/martservice", 
                       dataset="hsapiens_gene_ensembl"))
hk.known$hgnc_symbol
hk.known$description
write.csv(hk.known, paste0(root, "/data/RNAelements/figs/multiOmicBiomarkerPanelGenes.csv"))

# performance of splsda
cv = perf(result, validation = "loo")
cv$WeightedPredict.error.rate
# cv$error.rate

biomarkerGenescv = unlist(lapply(strsplit(names(unlist(cv$features$stable$Genes)), "\\."), function(i) i[2]))
write.csv(biomarkerGenescv, paste0(root, "/data/RNAelements/figs/bimarkerPanleGenes.csv"))

cv$predict$nrep1$Cells[[2]]

## plot area under the curve
predictScores = Reduce("+", lapply(cv$predict$nrep1, function(i) i[[2]]))/length(X) #average prediction
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T, direction = "<")   ## had the wrong direction!!
roc.res4 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res4$Specificity = 100 - as.numeric(roc.res4$Specificity)
roc.res4$Sensitivity = as.numeric(roc.res4$Sensitivity)
roc.score$auc

## Cells
predictScores = cv$predict$nrep1$Cells[[2]]
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T, direction = "<")
roc.res1 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res1$Specificity = 100 - as.numeric(roc.res1$Specificity)
roc.res1$Sensitivity = as.numeric(roc.res1$Sensitivity)
roc.score$auc

## Genes
predictScores = cv$predict$nrep1$Genes[[2]]
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T,
                direction = "<")
roc.res2 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res2$Specificity = 100 - as.numeric(roc.res2$Specificity)
roc.res2$Sensitivity = as.numeric(roc.res2$Sensitivity)
roc.score$auc

## Metabolites
predictScores = cv$predict$nrep1$Metabolites[[2]]
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T,
                direction = "<")
roc.res3 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res3$Specificity = 100 - as.numeric(roc.res3$Specificity)
roc.res3$Sensitivity = as.numeric(roc.res3$Sensitivity)
roc.score$auc

## plot AUC curves
pdf(paste0(root, "/data/RNAelements/figs/aucs.pdf"), height = 6, width = 6)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 2))
plot(roc.res1$Sensitivity ~ roc.res1$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", col = "#1F78B4", main = "Cells")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 74.5%", col = "#1F78B4")
plot(roc.res2$Sensitivity ~ roc.res2$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", col = "#33A02C", main = "Genes")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 59.2%", col = "#33A02C")
plot(roc.res3$Sensitivity ~ roc.res3$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", col = "#E31A1C", main = "Metabolites")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 88.6%", col = "#E31A1C")
plot(roc.res4$Sensitivity ~ roc.res4$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", main = "Combined")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 68.6%")

plot(X[[1]][, "dmap_HSC3"]~ X[[3]][, "xLeucine"])
abline(lm(X[[1]][, "dmap_HSC3"]~ X[[3]][, "xLeucine"]))
dev.off()



## performance on each dataset separately; do plsda again #############
# Cells
result.cc = splsda(X = SPVs, Y = Y, keepX = c(2, 2))
cv.cc = perf(result.cc, validation = "loo")
predictScores = cv.cc$predict$`comp 2`
roc.score = roc(response = as.character(Y), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T)
roc.res3 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res3$Specificity = 100 - as.numeric(roc.res3$Specificity)
roc.res3$Sensitivity = as.numeric(roc.res3$Sensitivity)
roc.score$auc

Genes = t(genDats$starEnsemblExp[is.na(match(rownames(genDats$starEnsemblExp), rownames(Dat4))), names(Y)])
# Genes
result.gen = splsda(X = Genes, Y = Y, keepX = c(10, 10))
cv.gen = perf(result.gen, validation = "loo")
predictScores = cv.gen$predict$`comp 2`
roc.score = roc(response = as.character(Y), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T)
roc.res3 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res3$Specificity = 100 - as.numeric(roc.res3$Specificity)
roc.res3$Sensitivity = as.numeric(roc.res3$Sensitivity)
roc.score$auc

met = asthmaExpPre[rownames(starDemo), ]
## Metabolomics
result.met = splsda(X = met, Y = Y, keepX = c(5, 5))
cv.met = perf(result.met, validation = "loo")
predictScores = cv.met$predict$`comp 2`
roc.score = roc(response = as.character(Y), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T)
roc.res3 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res3$Specificity = 100 - as.numeric(roc.res3$Specificity)
roc.res3$Sensitivity = as.numeric(roc.res3$Sensitivity)
roc.score$auc


result = block.splsda(X = X, Y = Y, ncomp = ncomp,
                      keepX = keepX, design = design,
                      mode = "regression")#bias = T)

# definition of the keepX value to be tested for each block mRNA miRNA and protein
# names of test.keepX must match the names of 'data'
test.keepX = list(Cells = seq(5,20,5), Genes = seq(5,20,5), Metabolites = seq(5,20,5))

# the following may take some time to run, note that for through tuning
# nrepeat should be > 1
tune = tune.block.splsda(X = X, Y = Y,
                         ncomp = ncomp, test.keepX = test.keepX, design = design, nrepeat = 1,
                         cpus = 8, validation = "Mfold")

tune$choice.keepX.constraint # NULL as constraint = F per default
tune$choice.keepX

pdf(paste0(root,"/data/rnaelements/figs/tune_blockplsda.pdf"))
par(mfrow=c(2,1))
plot(tune$error.rate[, "comp1"][order(tune$error.rate[, "comp2"])] ~ factor(rownames(tune$error.rate)), col = 1, ylim = c(0.2, 1))
points(tune$error.rate[, "comp2"][order(tune$error.rate[, "comp2"])] ~ factor(rownames(tune$error.rate)), col = 2)
dev.off()

which.min(tune$error.rate[, "comp2"])

# use tuned keepX
keepX = tune$choice.keepX
result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                      keepX = keepX, design = design,
                      mode = "regression")#, bias = T)
feat1 = lapply(result$loadings, function(x) apply(x, 2, function(i) names(i)[which(i != 0)]))


cv = perf(result, validation = "Mfold")
cv$WeightedPredict.error.rate
#cv$error.rate

biomarkerGenescv = unlist(lapply(strsplit(names(unlist(cv$features$stable$Genes)), "\\."), function(i) i[2]))
#write.csv(biomarkerGenescv, paste0(root, "/data/RNAelements/figs/bimarkerPanleGenes.csv"))

cv$predict$nrep1$Cells[[2]]

## plot area under the curve
predictScores = Reduce("+", lapply(cv$predict$nrep1, function(i) i[[2]]))/length(X)
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T, direction = "<")   ## had the wrong direction!!
roc.res4 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res4$Specificity = 100 - as.numeric(roc.res4$Specificity)
roc.res4$Sensitivity = as.numeric(roc.res4$Sensitivity)
roc.score$auc

## Cells
predictScores = cv$predict$nrep1$Cells[[2]]
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T, direction = "<")
roc.res1 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res1$Specificity = 100 - as.numeric(roc.res1$Specificity)
roc.res1$Sensitivity = as.numeric(roc.res1$Sensitivity)
roc.score$auc

## Genes
predictScores = cv$predict$nrep1$Genes[[2]]
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T,
                direction = "<")
roc.res2 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res2$Specificity = 100 - as.numeric(roc.res2$Specificity)
roc.res2$Sensitivity = as.numeric(roc.res2$Sensitivity)
roc.score$auc

## Metabolites
predictScores = cv$predict$nrep1$Metabolites[[2]]
roc.score = roc(response = factor(as.character(Y), levels = c("ER", "DR")), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T,
                direction = "<")
roc.res3 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res3$Specificity = 100 - as.numeric(roc.res3$Specificity)
roc.res3$Sensitivity = as.numeric(roc.res3$Sensitivity)
roc.score$auc





























## post-challenge data RNASEQ #################

asthmaDemoPost = asthmaDemo[asthmaDemo$Time == "Post", ]
asthmaExpDiff = t(asthmaExp[, rownames(asthmaDemoPost)]) - t(asthmaExp[, rownames(asthmaDemoPre)])

## load RNA-Seq datasets
load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_rawData.RDATA")
load("~/projects/asthma/data/RNAseq/allRnaseqDatasets_normalized.RDATA")
## use ensembl datasets
# remove ERCC controls
ensemblDat0 = normalizelibSum(starEnsemblExp[-(60156:nrow(starEnsemblExp)), colnames(genDats$starEnsemblExp)])
ensemblDat = ensemblDat0[, !is.na(match(colnames(ensemblDat0), rownames(asthmaDemoPost)))]-ensemblDat0[, !is.na(match(colnames(ensemblDat0), rownames(asthmaDemoPre)))]
dim(ensemblDat)
starDemo = asthmaDemoPost[colnames(ensemblDat), ]
all(rownames(starDemo) == colnames(ensemblDat ))

# ## get markers genes
# ## DMAP
# data(DMAP)
# data(IRIS)
# 
# dmapTag=tagData(DMAP, 0.5, max=15, ref=NULL, ref.mean=F)
# irisTag=tagData(IRIS, 1, max=15, ref=NULL, ref.mean=F);
# lm22 = read.delim("~/projects/asthma/data/LM22.txt", row.names = 1)
# lm22Tag=tagData(lm22, 0.5, max=15, ref=NULL, ref.mean=F)
# pancancer0 = read.csv("~/projects/asthma/data/Cells_nCounter_Human_PanCancer_Immune_Profiling_Panel_Gene_List.csv")
# dim(dmapTag); dim(irisTag); dim(lm22Tag);
# 
# dmapList = apply(dmapTag, 2, function(i) names(i)[which(i != 0)])
# irisList = apply(irisTag, 2, function(i) names(i)[which(i != 0)])
# lm22List = apply(lm22Tag, 2, function(i) names(i)[which(i != 0)])
# pancancerList = split(pancancer0[pancancer0$Cell.Type != "", c("HUGO.Name")], pancancer0[pancancer0$Cell.Type != "", c("Cell.Type")])
# 
# dmapDat = data.frame(Gene = unlist(dmapList), Cell.Type = paste("dmap",rep(names(dmapList), unlist(lapply(dmapList, length))), sep = "_"))
# irisDat = data.frame(Gene = unlist(irisList), Cell.Type = paste("iris",rep(names(irisList), unlist(lapply(irisList, length))), sep = "_"))
# lm22Dat = data.frame(Gene = unlist(lm22List), Cell.Type = paste("lm22",rep(names(lm22List), unlist(lapply(lm22List, length))), sep = "_"))
# pancancerDat = data.frame(Gene = unlist(pancancerList), Cell.Type = paste("pancancer",rep(names(pancancerList), unlist(lapply(pancancerList, length))), sep = "_"))
# 
# cellSpecificList0 = rbind(dmapDat, irisDat, lm22Dat, pancancerDat)
# cellSpecificList = split(cellSpecificList0$Gene, cellSpecificList0$Cell.Type)
# 
# ids = unique(unlist(cellSpecificList))
# Dat = do.call(cbind, lapply(cellSpecificList, function(i){
#   x = rep(0, length(ids))
#   x[ids %in% i] = 1
#   x
# }))
# rownames(Dat) = ids
# 
# ## annoate gene symbols to ensembl ids
# 
# ensemIDs = 
#   getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
#         filters = 'hgnc_symbol',
#         values = ids,
#         mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
#                        host="grch37.ensembl.org", 
#                        path="/biomart/martservice", 
#                        dataset="hsapiens_gene_ensembl"))

# cc = unique(ensemIDs$hgnc_symbol)
# Dat2 = Dat[cc, ]
# rowNames = unlist(lapply(strsplit(rownames(ensemblDat), "\\."), function(i) i[1]))
# ccDat = unlist(lapply(cc, function(x){
#   ensemblID = rownames(ensemblDat)[rowNames %in% subset(ensemIDs, hgnc_symbol == x)$ensembl_gene_id]
#   if(length(ensemblID) == 0){
#     ensemblID = NA
#   }
#   ensemblID
# }))
# Dat3 = Dat2[!is.na(ccDat), ]
# rownames(Dat3) = ccDat[!is.na(ccDat)]
# Dat4 = Dat3[, colSums(Dat3) > 3]

SPVs=getAllSPVs(ensemblDat, grp=starDemo$calculated_Response, 
                Dat4, method="mixed",
                plot=F, mix.par = 0.3)
rownames(SPVs) = rownames(starDemo)

pdf(paste0(root, "/data/RNAelements/figs/cellTypes_post.pdf"), width = 9, height = 9)
cim(cor(SPVs), margins = c(15, 15))
dev.off()

## post-challenge
Y = factor(starDemo$calculated_Response, levels = c("ER", "DR"))
names(Y) = rownames(starDemo)
X = list(Cells = SPVs[names(Y), feat1$Cells], Genes = t(genDats$starEnsemblExp[feat1$Genes, names(Y)]), 
         Metabolites = asthmaExpDiff[names(Y), feat1$Metabolites])
dim(X[[1]]); dim(X[[2]]); dim(X[[3]]); 

ncomp = rep(2, length(X))
# design = matrix(c(0, 0, 0,
#                   0, 0, 0,
#                   0, 0, 0), nrow = 3, ncol = 3)
keepX = list(rep(2, ncomp[1]), rep(10, ncomp[1]), rep(5, ncomp[1]))
result = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                      keepX = keepX, design = design,
                      mode = "regression")#, bias = T)
feat2 = lapply(result$loadings, function(x)
  apply(x, 2, function(i) names(i)[which(i != 0)]))

pdf(paste0(root, "/data/RNAelements/figs/SamplePlot_post.pdf"), width = 10)
# plotDiablo2(result, ncomp = 1, groupOrder = c("DR", "ER"))
plotDiablo(result, ncomp = 1, groupOrder = c("DR", "ER"))
dev.off()

pdf(paste0(root, "/data/RNAelements/figs/circosPlot_post.pdf"), width = 7)
# circosPlot_diabloModif(result, corThreshold = 0.7)
circosPlot(result, cutoff=.7)#corThreshold = 0.7)
dev.off()

pdf(paste0(root, "/data/RNAelements/figs/heatmap_post.pdf"), width = 7)
# heatmap_diablo(result, margins = c(2, 12))
cimDiablo(result, margins = c(2, 12))
dev.off()

hk.known = 
  getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
        filters = 'ensembl_gene_id',
        values = unlist(lapply(strsplit(as.character(feat1$genes), "\\."), function(i) i[1])),
        mart =  useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        host="grch37.ensembl.org", 
                        path="/biomart/martservice", 
                        dataset="hsapiens_gene_ensembl"))$hgnc_symbol
write.csv(unlist(lapply(strsplit(as.character(feat1$genes), "\\."), function(i) i[1])), paste0(root, "/data/RNAelements/figs/bimarkerPanleGenes.csv"))
unique(hk.known)

cv = perf(result, validation = "loo")
cv$AveragePredict.error.rate

#biomarkerGenescv = unlist(lapply(strsplit(names(unlist(cv$features$stable$genes)), "\\."), function(i) i[2]))
#write.csv(biomarkerGenescv, paste0(root, "/data/RNAelements/figs/bimarkerPanleGenes.csv"))


## plot area under the curve
predictScores = Reduce("+", lapply(cv$predict, function(i) i[[2]]))/length(X)
roc.score = roc(response = as.character(Y), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T)
roc.res4 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res4$Specificity = 100 - as.numeric(roc.res4$Specificity)
roc.res4$Sensitivity = as.numeric(roc.res4$Sensitivity)
roc.score$auc

## Cells
predictScores = cv$predict$Cells[[2]]
roc.score = roc(response = as.character(Y), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T)
roc.res1 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res1$Specificity = 100 - as.numeric(roc.res1$Specificity)
roc.res1$Sensitivity = as.numeric(roc.res1$Sensitivity)
roc.score$auc

## Genes
predictScores = cv$predict$Genes[[2]]
roc.score = roc(response = as.character(Y), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T)
roc.res2 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res2$Specificity = 100 - as.numeric(roc.res2$Specificity)
roc.res2$Sensitivity = as.numeric(roc.res2$Sensitivity)
roc.score$auc

## Metabolites
predictScores = cv$predict$Metabolites[[2]]
roc.score = roc(response = as.character(Y), predictor = predictScores[, "DR"], plot = T, percent = T, na.rm =T)
roc.res3 = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res3$Specificity = 100 - as.numeric(roc.res3$Specificity)
roc.res3$Sensitivity = as.numeric(roc.res3$Sensitivity)
roc.score$auc

## plot AUC curves
pdf(paste0(root, "/data/RNAelements/figs/aucs_validation_post.pdf"), height = 3, width = 10)
par(mfrow = c(1, 4), mar = c(4, 4, 2, 2))
plot(roc.res1$Sensitivity ~ roc.res1$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", col = "#1F78B4", main = "Cells")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 48.2%", col = "#1F78B4")
plot(roc.res2$Sensitivity ~ roc.res2$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", col = "#33A02C", main = "Genes")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 77.7%", col = "#33A02C")
plot(roc.res3$Sensitivity ~ roc.res3$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", col = "#E31A1C", main = "Metabolites")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 64.7%", col = "#E31A1C")
plot(roc.res4$Sensitivity ~ roc.res4$Specificity, type = "o", pch = 19, 
     xlab = "100-Specificity", ylab = "Sensitivity", main = "Combined")
abline(a = 0, b = 1, col = "gray", lty = 2)
text(x = 80, y = 20, labels = "AUC = 67.8%")
dev.off()

