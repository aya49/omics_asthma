## 23andme data analysis

libr("gwascat") # interface to the [NHGRI's](http://www.genome.gov/) database of gwas
libr("ggplot2")
libr("biomaRt") #ensembl
libr("rtracklayer") #ucsc
libr("TxDb.Hsapiens.UCSC.hg19.knownGene")
libr("org.Hs.eg.db") #org.* annotation packages; can forge own and interact with using library("AnnotationDbi")
libr("ggbio")

d = get(load("temp"))

d$chrom = ordered(d$chrom, levels=c(seq(1, 22), "X", "Y", "MT"))
ggplot(d) + geom_bar(aes(chrom))


#genomicranges annotation packages
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene #transcriptDb; behind the scenes, everythign is SQLite
class(txdb) ## do some digging around!
transcripts(txdb)


## gene look up
# group individual ranges into groups: transcriptsBy, cdsBy, exonsBy
tx.by.gene <- transcriptsBy(txdb, "gene") #list names are Entrez gene ID's
# gene names; Entrez Gene ID
tx.by.gene # GRangesList object
columns(org.Hs.eg.db)
# keys: APOE gene
select(org.Hs.eg.db, keys="APOE", columns=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="SYMBOL") #keytypes()
# look up Gene ID
tx.by.gene["348"]


## SNP fall in region APOE? NOTE: check that chromosome format is the same
levels(d$chrom) <- paste("chr", c(1:22, "X", "Y", "M"), sep="") #convert to bioconductor format
my.snps <- with(d, GRanges(seqnames=chrom, 
                           IRanges(start=position, width=1), 
                           rsid=rsid, genotype=genotype)) # this goes into metadata
apoe.i <- findOverlaps(tx.by.gene["348"], my.snps) #RangesMatching class; if don't give chr name, warning sequence names don't match
hits <- matchMatrix(apoe.i)[, "subject"]
hits #see output

my.snps[hits]
# can also open a UCSC browser to this spot from R using `rtracklayer`
# ApoE4 allele is rs429358(C) + rs7412(C)
# most common allele (ApoE3, or e3/e3) is rs429358(T) + rs7412(C)


## check meta data on these SNP; elementMetadata(gwrngs) (scan); Strongest.SNP.Risk.Allele (at risk for what)

# join SNP with gwrngs meta data
dm = gwrngs.emd <- as.data.frame(get(data(gwrngs38)))
dm <- merge(d, gwrngs.emd, by.x="rsid", by.y="SNPs")
risk.alleles <- gsub("[^\\-]*-([ATCG?])", "\\1", dm$Strongest.SNP.Risk.Allele)

i.have.risk <- mapply(function(risk, mine) {
  risk %in% unlist(strsplit(mine, ""))
}, risk.alleles, dm$genotype)
dm$i.have.risk <- i.have.risk

# see data summary e.g. $Risk.Allele.Frequency
my.risk <- dm[dm$i.have.risk, ]
# rel.cols <- c(colnames(d), "Disease.Trait", "Risk.Allele.Frequency",
#               "p.Value", "i.have.risk", "X95..CI..text.")
rel.cols <- c("Disease.Trait", "Risk.Allele.Frequency",
              "p.Value", "X95..CI..text.")
# make sure $Initial.Sample.Size is the desired RACE!
head(my.risk[grep("European", my.risk$Initial.Sample.Size), rel.cols], 30)
head(dm[grepl("European", dm$Initial.Sample.Size), rel.cols], 30)


## plot
# location of all SNPs that `gwascat`tells me my allele is the "risk" allele (well, some "Disease.Traits" are height)
p <- plotStackedOverview(txdb, cytoband=FALSE)

(data(gwrngs38)$my.genotype <- 
    d$genotype[(match(data(gwrngs38)$SNPs, d$rsid))])

data(gwrngs38)$my.risk <- with(data(gwrngs38), 
                                        mapply(function(risk, mine) {
                                          risk %in% unlist(strsplit(mine, ""))
                                        }, gsub("[^\\-]*-([ATCG?])", "\\1", Strongest.SNP.Risk.Allele), my.genotype))

p + geom_hotregion(data(gwrngs38), aes(color=my.risk))






