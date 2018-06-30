## snpEnrichment
devtools::install_github("mcanouil/snpEnrichment")
library(snpEnrichment)

### 1. prep data

snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment")
signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment")

initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)

writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL,
        ldThresh = 0.8, depth = 1000, mc.cores = 1)



### 2. read data

snpListDir <- system.file("extdata/List", package = "snpEnrichment")
data(transcript)
transcriptFile <- transcript

toyData <- readEnrichment(
  pattern = "Chrom", signalFile, transcriptFile, snpListDir, snpInfoDir,
  distThresh = 1000, sigThresh = 0.05, LD = TRUE, ldDir = NULL,
  mc.cores = 1
)
toyData



### 3. compute

reSample(object = toyData,
         nSample = 10,
         empiricPvalue = TRUE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = TRUE)




## exclude SNP from original list

excludeFile <- c(
  "rs4376885", "rs17690928", "rs6460708", "rs13061537", "rs11769827",
  "rs12717054", "rs2907627", "rs1380109", "rs7024214", "rs7711972",
  "rs9658282", "rs11750720", "rs1793268", "rs774568", "rs6921786",
  "rs1699031", "rs6994771", "rs16926670", "rs465612", "rs3012084",
  "rs354850", "rs12803455", "rs13384873", "rs4364668", "rs8181047",
  "rs2179993", "rs12049335", "rs6079926", "rs2175144", "rs11564427",
  "rs7786389", "rs7005565", "rs17423335", "rs12474102", "rs191314",
  "rs10513168", "rs1711437", "rs1992620", "rs283115", "rs10754563",
  "rs10851727", "rs2173191", "rs7661353", "rs1342113", "rs7042073",
  "rs1567445", "rs10120375", "rs550060", "rs3761218", "rs4512977"
)

# OR
excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment")
toyData_exclude <- excludeSNP(toyData, excludeFile, mc.cores = 1)
# in development
compareResults <- compareEnrichment(object.x = toyData,
                                    object.y = toyData_exclude,
                                    pattern = "Chrom",
                                    nSample = 10,
                                    empiricPvalue = TRUE,
                                    mc.cores = 1,
                                    onlyGenome = TRUE)




### visualize

show(toyData)
print(toyData)
head(getEnrichSNP(toyData, type = "xSNP"))

show(toyData_exclude)
print(toyData_exclude)
head(getEnrichSNP(toyData_exclude, type = "eSNP"))



