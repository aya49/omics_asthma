## incomplete

library("ggbio")
library("Homo.sapiens")
library("GenomicRanges")

## ideogram track : Plot single chromosome with cytoband

#chr 1 is automatically drawn by default (subchr="chr1")
p.ideo <- Ideogram(genome = "hg38")
p.ideo
#Highlights a region on "chr2"
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))



## gene model track
# Gene model is composed of genetic features CDS, UTR, introns, exons and non-genetic region

#load gene symbol : GRanges, one gene/row
data(genesymbol, package = "biovizBase")
#retrieve information of the gene of interest
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)
#Plot the different transcripts  for our gene of interest
p.txdb <- autoplot(Homo.sapiens, which = wh)
p.txdb
#Change inton geometry, use gap.geom
autoplot(Homo.sapiens, which = wh, gap.geom = "chevron")
#Change colours
autoplot(Homo.sapiens, which = wh, label.color = "black", color = "brown",
         fill = "brown")
columns(Homo.sapiens)

#Flexible label
autoplot(Homo.sapiens, which = wh, columns = c("GENENAME", "GO"), names.expr = "GENENAME::GO")


## Gene model from TranscriptDb object
# TranscriptDb doesn't contain any gene symbol information, so we use tx id as default for label.

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
autoplot(txdb, which = wh)


## Add a reference track
library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
p.bg <- autoplot(bg, which = wh)
## no geom
p.bg
## segment
p.bg + zoom(1/100)
## rectangle
p.bg + zoom(1/1000)
## text
p.bg + zoom(1/2500)


# To override a zemantic zoom threshold, you simply provide a geom explicitly.

library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
## force to use geom 'segment' at this level
autoplot(bg, which = resize(wh, width = width(wh)/2000), geom = "segment")



## To override a zemantic zoom threshold, you simply provide a geom explicitly.

library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
## force to use geom 'segment' at this level
autoplot(bg, which = resize(wh, width = width(wh)/2000), geom = "segment")

## Building your tracks
gr17 <- GRanges("chr17", IRanges(41234415, 41234569))
tks <- tracks(p.ideo, mismatch = p.mis, dbSNP = p.vr, ref = p.bg, gene = p.txdb,
              heights = c(2, 3, 3, 1, 4)) + xlim(gr17) + theme_tracks_sunset()
tks




## circular plots

#Load the data
data("CRC", package = "biovizBase")
head(hg19sub)

p <- ggbio() + circle(hg19sub, geom = "ideo", fill = "gray70") + #Ideogram
  circle(hg19sub, geom = "scale", size = 2) + #Scale
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3) # label
p #print plot

# show somatic mutation
head(mut.gr)
p <- ggbio() + circle(mut.gr, geom = "rect", color = "steelblue") + #somatic mutation
  circle(hg19sub, geom = "ideo", fill = "gray70") +#Ideogram
  circle(hg19sub, geom = "scale", size = 2) +#Scale
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)#label
p



