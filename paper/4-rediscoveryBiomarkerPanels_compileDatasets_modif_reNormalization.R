##############################################################
#
# 4-rediscoveryBiomarkerPanels_compileDatasets.R
# Date: February 26, 2016
#
#############################################################
WhereAmI <- "~/Dropbox/Asthma/biomarkerPanels/reAnalysis/"

## load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(psych)  # geometric.mean()
library(NormqPCR) # selectHKs
library(amritr)

source("~/Dropbox/Asthma/biomarkerPanels/code/discovery/functions.R")
load("~/Dropbox/Asthma/biomarkerPanels/data/discovery/rnaseq/allRnaseqDatasets_normalized.RDATA")

#----------------------------------------------
# 1) HBA attenutation (Discovery cohort)
#----------------------------------------------
filePaths <- c("~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_rediscovery/20160114_test0 max fovjan14-2016 _RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_rediscovery/20160116_set2 jan 16-2016 max fov_RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_rediscovery/20160116_set 3 3-1 to 3-12 jan16 2016 max fov_RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_rediscovery/20160116_set 4 4-1 to 4-12 maxfov jan 16-2016_RCC")
setName <- paste("Set", 1:length(filePaths))
date <- c("2016-01-14","2016-01-16","2016-01-16", "2016-01-16")
cartridgeLotNum <- rep("31020635468_exp02.2016", length(filePaths))
Tag156LotNum <- rep("TS2877-1", length(filePaths))
ext24TagLotNum <- rep("TS4006", length(filePaths))
plateLotNum <- rep("41000274-0056_exp02.2016", length(filePaths))
DatList <- list()
for(h in 1 : length(filePaths)){
  files <- list.files(filePaths[h], full.names = TRUE)
  dat <- lapply(as.list(files), function(i){
    rccToDat(fileName = i)
  }) %>% do.call(rbind, .)
  dat$Set <- setName[h]
  dat$Date <- date[h]
  dat$cartridgeLotNum <- cartridgeLotNum[h]
  dat$Tag156LotNum <- Tag156LotNum[h]
  dat$ext24TagLotNum <- ext24TagLotNum[h]
  dat$plateLotNum <- plateLotNum[h]
  DatList[[h]] <- dat
}
Dat <- do.call(rbind, DatList)
write.csv(Dat, paste0(WhereAmI, "data/HBA2Attenutation_asthmaBiomarkersDiscovery.csv"))

#------------
# import data 
#------------
sampleMap0 <- read.csv("~/Dropbox/RNA-SeqAnalysis/validationStudy/reCalibration/attenuation/AttenuationnanoElementsDemo_labSheet.csv")
sampleMap <- c(paste0("E", sampleMap0$EorV_number))
names(sampleMap) <- paste(sampleMap0$fileName, unlist(lapply(strsplit(as.character(sampleMap0$Set_lane), "_"), function(i) i[2])), sep = "_")
data <- read.csv(paste0(WhereAmI, "data/HBA2Attenutation_asthmaBiomarkersDiscovery.csv"), row.names = 1)
data$Sample <- sampleMap[gsub(".RCC", "", unlist(lapply(strsplit(as.character(data$fileName), "/"), function(i) i[length(i)])))]
GeneNames <- as.character(data$Name)
GeneNames[data$Accession == "NM_001788.5"] <- "SEPT7"
data$Name <- GeneNames
codeClass <- as.character(data$CodeClass)
codeClass[data$Name == "comp56964_c0_seq1_hk"] <- "Housekeeping_Code_Summary"
codeClass[data$Name == "comp56975_c0_seq1_hk"] <- "Housekeeping_Code_Summary"
codeClass[data$Name == "comp56957_c0_seq1_hk"] <- "Housekeeping_Code_Summary"
data$CodeClass <- factor(codeClass)
asthmaDis <- data

#----------------------------------------------
# 2) HBA attenutation (confirmatory phase)
#----------------------------------------------
filePaths <- c("~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/20160210_asthma biomarkersvalidation set1 v1-v12 maxFOV Feb-10-2016_RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/20160210_asthma biomarkersvalidation set2 v13-v24 maxFOV Feb-10-2016_RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/20160210_asthma biomarkersvalidation set3 v25-v36 maxFOV Feb-10-2016_RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/20160211_asthma biomarkersvalidation set4 v37-v48 maxFOV Feb-11-2016_RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/20160211_asthma biomarkersvalidation set5 v49-v60 maxFOV Feb-11-2016_RCC",
               "~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/20160211_asthma biomarkersvalidation set6 v61-v72 maxFOV Feb-11-2016_RCC")
setName <- paste("Set", 1:length(filePaths))
date <- rep(c("2016-02-10","2016-02-11"), each = 3)
cartridgeLotNum <- rep("31020635467_exp12.2016", length(filePaths))
Tag156LotNum <- rep("TS2877-1", length(filePaths))
ext24TagLotNum <- rep("TS4006", length(filePaths))
plateLotNum <- rep("41000274-0056_exp02.2016", length(filePaths))
DatList <- list()
for(h in 1 : length(filePaths)){
  files <- list.files(filePaths[h], full.names = TRUE)
  dat <- lapply(as.list(files), function(i){
    rccToDat(fileName = i)
  }) %>% do.call(rbind, .)
  dat$Set <- setName[h]
  dat$Date <- date[h]
  dat$cartridgeLotNum <- cartridgeLotNum[h]
  dat$Tag156LotNum <- Tag156LotNum[h]
  dat$ext24TagLotNum <- ext24TagLotNum[h]
  dat$plateLotNum <- plateLotNum[h]
  DatList[[h]] <- dat
}
Dat <- do.call(rbind, DatList)
write.csv(Dat, "~/Dropbox/Asthma/biomarkerPanels/data/HBA2Attenutation_asthmaBiomarkersValidation.csv")

#------------
# import data 
#------------
mappingFile <- read.csv("~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/validationCohort_mappngFile.csv", row.names = 7)
data <- read.csv("~/Dropbox/Asthma/biomarkerPanels/data/validation/HBA2_attenuation_replication/HBA2Attenutation_asthmaBiomarkersValidation.csv", row.names = 1)
sampleMap <- gsub(".RCC", "", paste(unlist(lapply(strsplit(as.character(data$fileName), " "), function(i) i[3])),
                                    unlist(lapply(strsplit(as.character(data$fileName), "_"), function(i) i[8])), sep = "_"))
data$Sample <- mappingFile[sampleMap, "V.ID"]
GeneNames <- as.character(data$Name)
GeneNames[data$Accession == "NM_001788.5"] <- "SEPT7"
data$Name <- GeneNames
codeClass <- as.character(data$CodeClass)
codeClass[data$Name == "comp56964_c0_seq1_hk"] <- "Housekeeping_Code_Summary"
codeClass[data$Name == "comp56975_c0_seq1_hk"] <- "Housekeeping_Code_Summary"
codeClass[data$Name == "comp56957_c0_seq1_hk"] <- "Housekeeping_Code_Summary"
data$CodeClass <- factor(codeClass)
asthmaVal <- data

asthmaDat <- rbind(asthmaDis, asthmaVal)
asthmaDat$Cohort <- rep(c("Discovery", "Confirmatory"), c(nrow(asthmaDis), nrow(asthmaVal)))

##################################################
## quality control
##################################################
qc <- asthmaDat %>% mutate(FOVratio = FovCounted_Lane_Attributes/FovCount_Lane_Attributes) %>% 
  dplyr::select(FOVratio, BindingDensity = BindingDensity_Lane_Attributes, Sample, Set, ID_Lane_Attributes, Cohort) %>%
  gather(QCVar, QCmetrics, -c(Sample, Set, ID_Lane_Attributes, Cohort)) %>% dplyr::group_by(Sample, QCVar, Set, Cohort) %>% 
  #distinct(Sample, Set, ID_Lane_Attributes, Cohort) %>% 
  dplyr::mutate(Sample_Lane = paste(Sample, ID_Lane_Attributes, sep="_"))
qc$Sample_Lane <- factor(qc$Sample_Lane, levels = unique(qc$Sample_Lane))

hline.data <- data.frame(yint = c(0.75, 0.05, 2.25),
  QCVar = c("FOVratio", rep("BindingDensity", 2)))
ggplot(data = qc, aes(x = Sample_Lane, y = QCmetrics)) + geom_point() + 
  facet_grid(QCVar~Cohort, scale = "free") + 
  geom_hline(data = hline.data, aes(yintercept = yint), lty = 2, col = "red") + 
  customTheme(sizeStripFont=15, xAngle=90, hjust = 1, 
    vjust = 0.5, xSize=5, ySize=15, xAxisSize=15, yAxisSize=15)

## 3) Internal Controls
## Positive Controls
posConc <- c(128, 32, 8, 2, 0.5, 0.125)
names(posConc) <- paste0(paste("POS", LETTERS[1:6], sep = "_"), paste0("(", posConc, ")"))
posDat <- asthmaDat %>% filter(CodeClass == "Positive_Code_Summary") %>% 
  mutate(posConc = posConc[as.character(Name)],
    Sample_Lane = paste(Sample, ID_Lane_Attributes, sep="_"))

mods <- posDat %>% 
  group_by(ID_Lane_Attributes, Set, Cohort) %>% do(., r2 = summary(lm(Count ~ posConc, data = .))$r.squared)
r.squared <- data.frame(ID_Lane_Attributes = mods$ID_Lane_Attributes, R2 = round(unlist(mods$r2), 2),
  Set = mods$Set)
r.squared$x <- 20
r.squared$y <- 100
all(unlist(mods$r2) > 0.9)

#### Negative Controls
negSD <- asthmaDat %>% dplyr::filter(CodeClass == "Negative_Code_Summary") %>% group_by(fileName, Set, Sample, Cohort) %>% dplyr::summarise(negThres = mean(Count)+2*sd(Count))
pos0.5 = posDat %>% filter(posConc == 0.5)
all(levels(negSD$Sample) == levels(pos0.5$Sample))   # TRUE
negPos <- as.data.frame(rbind(negSD[, c("Sample", "Set", "negThres", "Cohort")], setNames(pos0.5[, c("Sample", "Set", "Count", "Cohort")], c("Sample", "Set", "negThres", "Cohort"))))
negPos$Type = rep(c("Neg(Mean+2SD)","Pos0.5fM"), each=nrow(negSD))
negPos$Sample <- factor(negPos$Sample, levels = unique(negPos$Sample)[order(as.numeric(gsub("E","",unique(negPos$Sample))))])
ggplot(negPos, aes(x = Sample, y = negThres, color = Type)) + geom_point() + ylab("Count") + 
  customTheme(sizeStripFont=15, xAngle=90, hjust = 1, vjust = 0.5, xSize=5, ySize=8, xAxisSize=5, yAxisSize=8) + 
  facet_grid(Cohort~.)

## 4) Counts for each type of tag
asthmaDat$CodeClass <- gsub("_Code_Summary", "", asthmaDat$CodeClass)
meanFeature <- asthmaDat %>% group_by(CodeClass, Name) %>% dplyr::summarise(mean = mean(Count)) %>% dplyr::arrange(mean)
asthmaDat$Name <- factor(as.character(asthmaDat$Name), levels = as.character(meanFeature$Name))

asthmaDat %>% group_by(fileName) %>% dplyr::summarise(TotalSum = sum(Count), Cohort=unique(Cohort), Sample=unique(Sample)) %>% 
  dplyr::select(TotalSum, Cohort) %>% 
  ggplot(aes(x = TotalSum, fill = Cohort)) + geom_density(alpha = 0.5) + 
  customTheme(sizeStripFont=15, xAngle=0, hjust = 0.5, vjust = 0.5, xSize=10, 
    ySize=10, xAxisSize=10, yAxisSize=10) +
  theme(legend.position = c(0.7,0.9))
## plot assay counts
Dat2 <- asthmaDat
Dat2$Name <- as.character(Dat2$Name)
Dat2$Name[Dat2$Name != "HBA2"] <- "other genes"
Dat3 <- Dat2 %>% group_by(fileName, Name) %>% dplyr::summarise(TotalSum = sum(Count), Cohort=unique(Cohort)) %>% 
  dplyr::select(TotalSum, Cohort, Name) 

Dat3 %>% 
  ggplot(aes(x = Name, y = TotalSum)) + geom_boxplot() +
  facet_grid(.~Cohort) +
  customTheme(sizeStripFont=9, xAngle=0, hjust = 0.5, vjust = 0.5, xSize=7, 
    ySize=10, xAxisSize=10, yAxisSize=10) + 
  scale_y_log10(breaks = c(1e5,2e5,seq(0,2e6, 5e5))) + xlab("")+
  ylab("Total counts for each assay") +
  theme(legend.position = c(0.5,0.2), legend.title=element_blank(),
    legend.text=element_text(size=7), legend.key=element_blank())

## plot assay counts per gene
ggplot(asthmaDat, aes(x = Name, y = Count, color = CodeClass)) + geom_point() + scale_y_log10() + 
  customTheme(sizeStripFont=15, xAngle=90, hjust = 1, vjust = 0.5, xSize=6, ySize=8, xAxisSize=8, yAxisSize=8) + 
  theme(legend.position = c(0.1, 0.85)) + xlab("CodeClass") +
  annotate("text", label = "HBA2", x = 100, y = 1500000, size = 4, colour = "black")


################################################################
#
# Normalization
#
################################################################
table(as.character(asthmaDat$Sample), as.character(asthmaDat$Date))
posConc <- c(128, 32, 8, 2, 0.5, 0.125)
names(posConc) <- paste0(paste("POS", LETTERS[1:6], sep = "_"), paste0("(", posConc, ")"))
posDat <- asthmaDat %>% filter(CodeClass == "Positive") %>% 
    mutate(posConc = posConc[as.character(Name)],
      Sample_Lane = paste(Sample, ID_Lane_Attributes, sep="_"))
  
## Positive Control Normalization
posNF <- posDat %>% group_by(fileName) %>% dplyr::summarise(posGeoMean = geometric.mean(Count)) %>% 
    dplyr::mutate(posNF = mean(posGeoMean)/posGeoMean)
all(posNF$posNF<3 & posNF$posNF > 0.3)
  
dataNorm <- full_join(asthmaDat, posNF, by = "fileName") %>% 
  mutate(posNorm = Count*posNF)
  
# check if the sum of the positive controls is the same for each assay
dataNorm %>% group_by(fileName) %>% filter(CodeClass == "Positive") %>% dplyr::summarise(sumPosSample = geometric.mean(posNorm), Sample = unique(Sample), Set = unique(Set)) %>% 
    ggplot(aes(x = Sample, y = sumPosSample)) + geom_point() + theme_bw()

## turn data into matrix of expression
dat <- dataNorm %>% dplyr::select(Sample, posNorm, Name, CodeClass) %>% 
       group_by(Sample, Name, CodeClass) %>% 
       summarise_all(funs(mean))
  
## house-keeper normalization
hk <- filter(dat, CodeClass == "Housekeeping")
hkMean <- hk %>% 
  dplyr::group_by(Name) %>% 
  dplyr::summarise(meanCount = mean(posNorm))
hk$Name <- factor(as.character(hk$Name), levels = as.character(hkMean$Name))
hk$Sample <- factor(as.character(hk$Sample), levels = unique(as.character(hk$Sample))[order(as.numeric(gsub("E|V","", unique(as.character(hk$Sample)))))])

ggplot(hk, aes(x = as.numeric(Sample), y = posNorm, color = Name)) + geom_point() + geom_line() +
  scale_y_log10() + customTheme(sizeStripFont=15, xAngle=0, hjust = 0.5, vjust = 0.5, xSize=10, ySize=15, xAxisSize=15, yAxisSize=15) + xlab("36 samples") + ylab("Counts")

hkDat0 <- hk %>% dplyr::select(Name, Sample, posNorm) %>% spread(Name, posNorm)
hkDat <- as.matrix(hkDat0[, -1])
hkDat <- log2(hkDat)
rownames(hkDat) <- as.character(hkDat0$Sample)

## Use geNorm to select house-keepers
res.BM <- selectHKs(hkDat, method = "geNorm", Symbols = colnames(hkDat), minNrHK = 2, log = FALSE)
hkRanking <- data.frame(Gene = as.character(res.BM$ranking), 
  Rank = 1:length(res.BM$ranking)) #,
#  Dataset = c(rep("UCSC genes", 3), rep("Trinity", 3)))

cvList <- list()
for(i in 3 : nrow(hkRanking)){
  genes <- as.character(hkRanking$Gene[1:i])
  
  ## House-keeping normalization factor
  hkNF <- dat %>% 
    filter(Name %in% genes) %>% 
    group_by(Sample) %>% dplyr::summarise(hkGeoMean = geometric.mean(posNorm)) %>% 
    mutate(hkNF = mean(hkGeoMean)/hkGeoMean)
  
  ## normalize data
  hkNorm <- full_join(dataNorm, hkNF, by = "Sample") %>% mutate(hkNorm = posNorm*hkNF)
  
  ## calculate CVs
  cvDat <- hkNorm %>% filter(CodeClass == "Endogenous") %>% 
    group_by(Name) %>% dplyr::summarise(cv = sd(hkNorm)/mean(hkNorm))
  cv <- cvDat$cv
  names(cv) <- cvDat$Name
  cvList[[i-2]] <- cv
}
names(cvList) <- paste0("hk", "1-", 3:nrow(hkRanking))
cvDatList <- as.data.frame(do.call(cbind, cvList))

## positive control normalization
cvDat <- dataNorm %>% filter(CodeClass == "Endogenous") %>% 
  group_by(Name) %>% dplyr::summarise(cv = sd(posNorm)/mean(posNorm))
cv <- cvDat$cv
names(cv) <- cvDat$Name
cvDatList$PosNorm <- cv
## raw data
cvDat <- data %>% filter(CodeClass == "Endogenous_Code_Summary") %>% 
  group_by(Name) %>% dplyr::summarise(cv = sd(Count)/mean(Count))
cv <- cvDat$cv
names(cv) <- cvDat$Name
cvDatList$Raw <- cv

boxplot(cvDatList, ylab = "coefficient of variation")
lapply(cvDatList, mean)

cvDatList %>% as.data.frame %>% gather(hkComb, cv) %>% 
  ggplot(aes(x = cv, color = hkComb)) + geom_density(alpha = 0.5)


##### plot normalize with top 3 house-keeping genes
genes <- as.character(hkRanking$Gene[1:3])

## House-keeping normalization factor
hkNF <- dat %>% 
  filter(Name %in% genes) %>% 
  group_by(Sample) %>% dplyr::summarise(hkGeoMean = geometric.mean(posNorm)) %>% 
  mutate(hkNF = mean(hkGeoMean)/hkGeoMean)

## normalize data
hkNorm <- full_join(dat, hkNF, by = "Sample") %>% mutate(hkNorm = posNorm*hkNF) %>% 
          dplyr::filter(CodeClass == "Endogenous")

eset0 <- hkNorm %>% dplyr::select(Name, Sample, hkNorm) %>% 
        spread(Name, hkNorm)
eset <- as.matrix(eset0[, -1])
eset <- log2(eset)
rownames(eset) <- as.character(eset0$Sample)

## PCA plot of expression data and the relationships of PCs with demographics
vars <- setdiff(colnames(demo), c("SPONSOR", "DRUG", "ID", "NAME", "AIC_YMD", "rnaseqID_quebec",
  "Allergen", "UniqueID", "SubjInitials_mislabelled", "MST_LST_numbers", 
  "Box", "concealedID"))
compVar(demo = demo[colnames(combat_eset), ], eset = t(combat_eset), variables = vars, ncomp = 10)

demo[setdiff(rownames(demo), rownames(demoDis)), c("Leukocyte_Counts.x10.9.", "Neu_percent", "lym_percent", "mono_percent", "eos_percent", "baso_percent")]

date <- as.character(dataNorm$Date)
names(date) <- dataNorm$Sample

## PCA
result <- mixOmics::pca(eset, scale = TRUE, center = TRUE)
plotIndiv(result, group = date[rownames(eset)], ellipse = TRUE,
  star = TRUE, ind.names = FALSE, legend = TRUE)

## Correct for batch
library(sva)
## correct for batch effect
batch <- factor(unlist(lapply(strsplit(date[rownames(eset)], "-"), function(i) i[2])))
modcombat = as.data.frame(matrix(1, nr = nrow(eset), nc = 1))
rownames(modcombat) <- names(batch)
all(names(batch) == rownames(eset))
combat_eset = ComBat(dat=t(eset), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

## before batch correction
result <- mixOmics::pca(eset, scale = TRUE, center = TRUE)
plotIndiv(result, group = date[rownames(eset)], ellipse = TRUE,
  star = TRUE, ind.names = FALSE, legend = TRUE)

## after batch correction
result <- mixOmics::pca(t(combat_eset), scale = TRUE, center = TRUE)
plotIndiv(result, group = batch, ellipse = TRUE,
  star = TRUE, ind.names = FALSE, legend = TRUE)

## compare samples with response
demo0 <- readRDS("~/Dropbox/Asthma/biomarkerPanels/data/demo/allsitesDemo.rds")
demo <- demo0[as.character(demo0$EorV_number) %in% rownames(eset), ]
rownames(demo) <- as.character(demo$EorV_number)
demo <- demo[rownames(eset), ]
all(rownames(eset) == rownames(demo))
dim(demo)

## compute area under the curve
require(pracma)
percentfev1 <- scale(t(demo[, c("BLFEV", "F10L", "F20L", "F30L","F45L", "F60L", "F90L", "F120L", "F180L", "F240L", "F300L", "F360L", "F420L")]),
  center = demo$BLFEV, scale = demo$BLFEV)
time <- c(0, 10, 20, 30, 45, 60, 90, 120, 180, 240, 300, 360, 420)
demo$AUC <- apply(percentfev1, 2, function(i){
  trapz(time, -i)
})


table(as.numeric(demo$rnaseq_coreSet_biomarkerAnalysis))
disDemo <- subset(demo, rnaseq_coreSet_biomarkerAnalysis == "Y")
table(disDemo$calculated_Response)
## number of unique individuals
all(table(droplevels(disDemo$NAME)) == 1)  ## all are unique

## validation cohort
valDemo <- subset(demo, rnaseq_coreSet_biomarkerAnalysis != "Y")
nrow(valDemo)
## number of unique subjects
length(table(droplevels(valDemo$NAME)))
## healthy controls
valDemo.hc <- subset(valDemo, calculated_Response == "Control")
nrow(valDemo.hc); length(table(droplevels(valDemo.hc$NAME)))
valDemo.er <- subset(valDemo, calculated_Response == "ER")
nrow(valDemo.er); length(table(droplevels(valDemo.er$NAME)))
valDemo.dr <- subset(valDemo, calculated_Response == "DR")
nrow(valDemo.dr); length(table(droplevels(valDemo.dr$NAME)))

intersect(as.character(disDemo$NAME), as.character(valDemo$NAME))

plot(disDemo$LAR, col = factor(disDemo$calculated_Response))
abline(h = -15, lty = 2)

plot(demo$LAR, col = factor(demo$calculated_Response))
library(lattice)
xyplot(demo$LAR ~ demo$EAR | rnaseq_coreSet_biomarkerAnalysis, 
  group = factor(demo$calculated_Response), data = demo)

plot(LAR ~ AUC, data = demo)
abline(lm(LAR ~ AUC, data = demo))
summary(lm(LAR ~ AUC+NAME, data = demo))

ggplot(demo, aes(x = NAME, y = LAR, fill = calculated_Response,
  color = calculated_Response,
  group = NAME, pch = rnaseq_coreSet_biomarkerAnalysis)) + 
  geom_point(size=3) + geom_line() +
  geom_text(data = demo, aes(x=NAME, y=LAR+1, 
    label=Allergen_cleanLabel, color = Allergen_cleanLabel)) +
  geom_hline(yintercept = -15, linetype = "dashed") +
  customTheme(sizeStripFont = 10, xAngle = 90, hjust = 1, vjust = 0.5, xSize =8, ySize = 10, xAxisSize = 10, yAxisSize = 10)

###










## flippers
multiples <- names(table(droplevels(demo$NAME)))[table(droplevels(demo$NAME)) >1]
multiples <- rowSums(table(droplevels(demo[demo$NAME %in% multiples, "NAME"]), 
  demo[demo$NAME %in% multiples, "calculated_Response"]) == 0)
multiples[multiples == 1]
demo.flipper <- subset(demo, NAME %in% names(multiples[multiples == 1]))

## how many flippers in each cohort
length(intersect(rownames(disDemo), rownames(demo.flipper))) # discovery
disDemo[intersect(rownames(disDemo), rownames(demo.flipper)), "NAME"]

length(intersect(rownames(valDemo.hc), rownames(demo.flipper))) # validation - controls

length(intersect(rownames(valDemo.er), rownames(demo.flipper))) # validation - er
valDemo.er[intersect(rownames(valDemo.er), rownames(demo.flipper)), "NAME"]

length(intersect(rownames(valDemo.dr), rownames(demo.flipper))) # validation - dr
valDemo.dr[intersect(rownames(valDemo.dr), rownames(demo.flipper)), "NAME"]

# analysis 1
load("~/Dropbox/Asthma/biomarkerPanels/data/validation/discovery&validation_expressionDatasets_Demo.RDATA")
biomarkers <- colnames(preReDis$ensembl)

X.train <- t(combat_eset)[rownames(disDemo), ]
Y.train = factor(disDemo$calculated_Response, c("ER", "DR"))
X.test <- t(combat_eset)[rownames(valDemo), ]
Y.test = valDemo$calculated_Response
Y.test[Y.test == "Control"] <- "DR"
Y.test = factor(Y.test, c("ER", "DR"))

fit <- enet(X = X.train, Y = Y.train, alpha = 0, lambda = NULL, family = "binomial", X.test = X.test,
  Y.test = Y.test, filter = "none", topranked = 50, keepVar = NULL)
cv <- perf.enet(object=fit, validation = "Mfold", M = 5, iter = 5, threads = 5, progressBar = FALSE)
cv$perf
fit$perfTest

plot(fit$probs ~ factor(valDemo$calculated_Response))
plot(fit$probs, col = factor(valDemo$calculated_Response), pch = 19)





table(demo.flipper$calculated_Response, demo.flipper$givenResponse)
table(demo.flipper$Allergen_cleanLabel, demo.flipper$calculated_Response, droplevels(demo.flipper$NAME))

demo.replicatesSampleResponse <- subset(demo, NAME %in% names(multiples[multiples > 1]))
table(demo.replicatesSampleResponse$calculated_Response, demo.replicatesSampleResponse$givenResponse)
demo.keep <- subset(demo, !(NAME %in% names(multiples)))
table(demo.keep$calculated_Response, demo.keep$givenResponse, demo.keep$Cohort)

## which samples have demographics data missing
vars <- c("SITE", "RACE", "SEX", "AGE",
  "Allergen_cleanLabel", "Leukocyte_Counts.x10.9.", "Neu_percent", "lym_percent", "mono_percent",
  "eos_percent", "baso_percent", "PC20", "NAME", "calculated_Response", "rnaseq_coreSet_biomarkerAnalysis")

colSums(is.na(demo[, vars]))
demoComplete <- na.omit(demo[intersect(rownames(demo), rownames(eset)), vars])


disDemo <- subset(demoComplete, rnaseq_coreSet_biomarkerAnalysis == "Y")
table(disDemo$calculated_Response)
## number of unique individuals
all(table(droplevels(disDemo$NAME)) == 1)  ## all are unique

## validation cohort
valDemo <- subset(demoComplete, rnaseq_coreSet_biomarkerAnalysis != "Y")
nrow(valDemo)
## number of unique subjects
length(table(droplevels(valDemo$NAME)))
## healthy controls
valDemo.hc <- subset(valDemo, calculated_Response == "Control")
nrow(valDemo.hc); length(table(droplevels(valDemo.hc$NAME)))
valDemo.er <- subset(valDemo, calculated_Response == "ER")
nrow(valDemo.er); length(table(droplevels(valDemo.er$NAME)))
valDemo.dr <- subset(valDemo, calculated_Response == "DR")
nrow(valDemo.dr); length(table(droplevels(valDemo.dr$NAME)))































table(demoComplete$rnaseq_coreSet_biomarkerAnalysis, demoComplete$calculated_Response)
disDemo[setdiff(rownames(disDemo), rownames(demoComplete)), vars]

demoComplete <- na.omit(demo[intersect(rownames(demo), rownames(eset)), c("SITE", "RACE", "SEX", "Wt..Kg.", "HT.cm.", "AGE", "BLFEV",
  "Allergen_cleanLabel", "Leukocyte_Counts.x10.9.", "Neu_percent", "lym_percent", "mono_percent",
  "eos_percent", "baso_percent", "PC20", "NAME", "calculated_Response", "rnaseq_coreSet_biomarkerAnalysis")])

table(demoComplete$rnaseq_coreSet_biomarkerAnalysis)

table(demoComplete$calculated_Response, droplevels(demoComplete$NAME))
table(demoComplete$calculated_Response, demoComplete$Cohort)

compVar(demo = demoComplete, eset = t(combat_eset[, rownames(demoComplete)]), variables = colnames(demoComplete), ncomp = 10)

## before batch correction
result <- mixOmics::pca(eset[rownames(demoComplete), ], scale = TRUE, center = TRUE)
plotIndiv(result, group = demoComplete$SITE, legend = TRUE)

## remove site effect
## correct for batch effect
batch <- droplevels(demoComplete$SITE)
modcombat = model.matrix(~1, data = demoComplete)
combat_eset2 = ComBat(dat=t(eset[rownames(demoComplete), ]), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

compVar(demo = demoComplete, eset = t(combat_eset2[, rownames(demoComplete)]), variables = colnames(demoComplete), ncomp = 10)

demoDis <- subset(demoComplete, Cohort == "Discovery")
demoVal <- subset(demoComplete, Cohort == "Validation")

setdiff(colnames(combat_eset2), rownames(demoDis))
setdiff(rownames(demoDis), colnames(combat_eset2))

load("~/Dropbox/Asthma/biomarkerPanels/data/validation/discovery&validation_expressionDatasets_Demo.RDATA")
biomarkers <- colnames(preReDis$ensembl)

X.train <- cbind(data.matrix(demoDis[, c("SEX", "Wt..Kg.", "HT.cm.", "AGE", "BLFEV", "Allergen_cleanLabel",
           "Leukocyte_Counts.x10.9.", "Neu_percent", "lym_percent", "mono_percent", "eos_percent", "baso_percent", "PC20")]),
  t(combat_eset2)[rownames(demoDis), biomarkers])
Y.train = factor(demoDis$calculated_Response, c("ER", "DR"))
X.test <- cbind(data.matrix(demoVal[, c("SEX", "Wt..Kg.", "HT.cm.", "AGE", "BLFEV", "Allergen_cleanLabel",
  "Leukocyte_Counts.x10.9.", "Neu_percent", "lym_percent", "mono_percent", "eos_percent", "baso_percent", "PC20")]),
  t(combat_eset2)[rownames(demoVal), biomarkers])
Y.test = factor(demoVal$calculated_Response, c("ER", "DR"))

fit <- enet(X = X.train, Y = Y.train, alpha = 0, lambda = NULL, family = "binomial", X.test = X.test,
  Y.test = Y.test, filter = "none", topranked = 50, keepVar = NULL)
cv <- perf.enet(object=fit, validation = "Mfold", M = 5, iter = 5, threads = 5, progressBar = FALSE)
cv$perf
fit$perfTest

plot(fit$probs ~ factor(demoVal$calculated_Response))
plot(fit$probs, col = factor(demoVal$calculated_Response), pch = 19)

performance(weights=as.numeric(fit$probs), trueLabels=factor(demoVal$calculated_Response, levels(Y.train)), pop.prev=0.6)

## Compute sensitivity and specificity
calcPerf = function(pred, truth, prev){
  mat <- table(pred, truth)
  if(rownames(mat)[1] != colnames(mat)[1])
    stop("check levels of inputs")
  
  spec = mat[1,1]/sum(mat[, 1])
  sens = mat[2,2]/sum(mat[, 2])
  npv = ((1-prev)*spec)/(prev*(1-sens)+(1-prev)*spec)
  ppv = (prev*sens)/((prev*sens)+(1-prev)*spec)
  perf = c(sens, spec, npv, ppv)
  names(perf) <- c("sens", "spec", "npv", "ppv")
  perf
}
## use a 0.5 decision boundary
testLabel <- rep(NA, nrow(X.test))
testLabel[fit$probs > 0.5] <- "DR"
testLabel[fit$probs < 0.5] <- "ER"
names(testLabel) <- rownames(X.test)

table(pred = factor(testLabel, levels(Y.test)), truth = Y.test)
calcPerf(pred = factor(testLabel, levels(Y.test)), truth = Y.test, prev = 0.6)






####################
##
# discovery and validation
#
####################
demoDis <- subset(demo.keep, Cohort == "Discovery")
table(demoDis$calculated_Response)
demoVal <- demo[setdiff(rownames(demo), rownames(demoDis)), ]
dim(demoDis); dim(demoVal);

fit <- enet(X = t(combat_eset)[rownames(demoDis), ], Y = factor(demoDis$calculated_Response, c("ER", "DR")),
           alpha = 0, lambda = NULL, family = "binomial", X.test = NULL,
           Y.test = NULL,
           filter = "none", topranked = 50, keepVar = NULL)
cv <- perf.enet(object=fit, validation = "Mfold", M = 5, iter = 5, threads = 5, progressBar = FALSE)
cv$perf

demoVal <- subset(demo, Cohort == "Validation", calculated_Response != "Control")
pred <- predict(fit$fit, newx = eset[rownames(demoVal), ], type = "response", s = fit$lambda)

plot(pred, col = as.numeric(factor(demoVal$calculated_Response)), pch = 19)

plot(pred ~ factor(demoVal$calculated_Response))
summary(lm(pred ~ factor(demoVal$calculated_Response)))

















