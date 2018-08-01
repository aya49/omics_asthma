## helper functions
## aya43@sfu.ca
## created 20180517
## last modified 20180517



## input: Sys.time() value
## output: time elapsed since input
time_output = function(start, tz="GMT", message="") {
  start = as.POSIXct(start)
  end = Sys.time()
  time_elapsed = difftime(end, start, units="secs")
  cat(message, ifelse(message=="","",": "), ft(start,tz=tz), "-", ft(end,tz=tz), ">", ft(time_elapsed,tz=tz), "\n")
}

## input: Sys.time() value
## output: formatted time as string; used in time_output function
ft = function(time, tz="GMT") {
  return( format(.POSIXct(time,tz=tz), "%H:%M:%S") )
}

## input: x=vector of indices; n=cores or number of parts you want to split x into
## ouput: list of n indices vector
loop_ind_f <- function(x,n) {
  if (n==1) return(list(x))
  return(split(x, ceiling(seq_along(x)/ceiling(length(x)/n))))
}

## input: matrix
## output: returns u1 (col index with only 1 unique element), ua (col index where every row is a unique element), prints colnames and its unique elements if there is less than n unique elements
col_probe = function(m,n=15) {
  u1 = NULL
  ua = NULL
  nm = nrow(m)
  for (col in 1:ncol(m)) {
    coln = colnames(m)[col]
    if (is.data.table(m)) {
      a = unique(as.matrix(m[,..col]))
    } else {
      a = as.vector(unique(as.matrix(m[,col])))
    }
    la = length(a)
    if (la == 1) u1 = append(u1,col)
    if (la == nm) ua = append(ua,col)
    
    cat("\n", coln, ": ", la)
    if (la<n) cat("; ",a)
  }
  return(list(u1=u1,ua=ua))
}






## input: file path and a file extension 
## output: List of all file names in given path with the given extension
fileNames = function(pathList, ext="fcs") {
  temp.str = sapply(str_split(pathList, "/"), function(x) x[length(x)])
  pathList = sub(paste0(".", ext,"$"), "", temp.str, ignore.case=T)
  return(pathList)
}




## input: x=vector of indices; n=cores or number of parts you want to split x into
## ouput: list of n indices vector
loopInd <- function(x,n) {
  if (n==1) return(list(x))
  return(split(x, ceiling(seq_along(x)/ceiling(length(x)/n))))
}



## input: list of package names to load
## output: none; load/install package
libr <- function(pkgs) {
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) 
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}


## by amrit
function(genExp) {
  lib.size <- rowSums(genExp)
  genExpNorm <- t(log2(t(genExp + 0.5)/(lib.size + 1) * 1e+06))
  return(genExpNorm)
}

## by amrit
function (fileName) {
  library(dplyr)
  library(tidyr)
  lines <- data.frame(values = readLines(fileName))
  dat <- suppressWarnings(separate(data = lines, col = values, 
                                   sep = ",", into = c("CodeClass", "Name", "Accession", 
                                                       "Count")))
  ind <- grep("<[A-Z]", dat$CodeClass)
  attr <- rep(NA, nrow(dat))
  for (i in 1:length(ind)) attr[ind[i]:nrow(dat)] <- grep("<[A-Z]", 
                                                          dat$CodeClass, value = TRUE)[i]
  dat <- dat %>% mutate(CodeClass = paste(CodeClass, gsub(" ", 
                                                          "", chartr("<>", "  ", attr)), sep = "_"), fileName = fileName)
  dat <- dat[-grep("<", dat$CodeClass), ]
  dat <- dat[!is.na(dat$Name), ]
  techDat <- dat[1:(grep("CodeClass", dat$CodeClass) - 1), 
                 ] %>% dplyr::select(-c(Accession:Count)) %>% spread(CodeClass, 
                                                                     Name)
  bioDat <- dat[(grep("CodeClass", dat$CodeClass) + 1):nrow(dat), 
                ]
  Dat <- full_join(techDat, bioDat, by = "fileName")
  return(Dat)
}


## input: matrix
## output: matrix without all NA col/row
delna <- function(m)
  m[apply(m,1,function(x) !all(is.na(x))), apply(m,2,function(x) !all(is.na(x)))]






## input:
# chrom: chromosome
# pos: SNP position
# val: association p-value
# val_max1: where "genomewide significance" threshold should be drawn
# val_max2: a sub-genomewide-sig "grey zone" where SNPs are shown with a larger point size
# draw_line: draw thresholds? T/F
# val_min: p-vlaues less than val_min are forced equal to val_min
# label.x: threshold; if there's more than label.x elements in val, label xx axis

## output: none; manhattan plot
## adapted from http://bioinfo-mite.crb.wsu.edu/Rcode/wgplot.R

manhattan_plot = function(val, chrom=rep(1,length(val)), pos=c(1:length(val)), label.x=40,
                          val_thres=-log(.025), val_max2=-log(1e-4), val_max1=-log(5e-8), 
                          draw_line=T, val_min=0, guideline_interval=1, lines=NULL,
                          xlab="chromosome/position", ylab=paste0("-log10(p value) (thres=",round(val_thres,3),")"), 
                          main="gwas") {
  # val's name can be label
  
  ## prep input
  val = as.numeric(val)
  val[val<val_min] = val_min
  val[val>val_max1] = val_max1
  pos = as.numeric(pos)
  chrom = as.character(chrom)
  
  chrom = as.character(chrom)
  chroml = chrom
  chroml[chroml=="23"]="X"
  chroml[chroml=="24"]="Y"
  chroml[chroml=="25"]="XY"
  chroml[chroml=="26"]="MT" 
  
  chrom[chrom=="X"]="23"
  chrom[chrom=="Y"]="24"
  chrom[chrom=="XY"]="25"
  chrom[chrom=="MT"]="26"
  
  ord = order(as.numeric(chrom),pos)
  
  val = val[ord]
  pos = pos[ord]
  chrom = chrom[ord]
  
  chrom_unique = as.character(sort(as.numeric(unique(chrom))))
  
  chrom_un = length(chrom_unique)
  chrom_table = table(as.numeric(chrom))
  chrom_cumsum = cumsum(chrom_table)
  
  ## get colours
  require(colorspace)
  colour = rainbow_hcl(chrom_un)
  # colour = rainbow_hcl(25)
  # colour = rep(colour, ceiling(chrom_un/length(colour)))
  
  # p = -log(p,10)
  
  ## make positions cumulatve
  if ( any(diff(pos)<0) ) {
    pos_cum = cumsum(c(0,pos[which(!duplicated(chrom))-1]))
    pos = pos + rep(pos_cum, chrom_table)
    mids = pos_cum + diff(c(pos_cum,max(pos)))/2
  }
  
  guidelines = seq(val_min, val_max1, guideline_interval)
  
  
  
  # ## plot
  # require(ggplot2)
  # inds = val>val_thres
  # df = data.frame(val=val, pos=pos, chrom=chrom, sig=inds)
  # g = ggplot(df, aes(x=pos, y=val)) +
  #   geom_point(aes(colour=chrom)) + #, shape=factor(inds))) + #scatterplot
  #   # geom_point(data=df[inds,], aes(x=pos, y=val), inherit.aes=F)
  #   coord_cartesian(ylim=c(0,9)) + #zoom in # + ylim(c(0,9)) #delete points
  #   geom_hline(yintercept=val_max, linetype="dashed", color = "red") +
  #   labs(title=main, subtitle="manhattan plot", y=ylab, x=xlab)
  # # plot(g)
  # return(g)
  
  
  
  ## plot
  par(xaxt = "n", yaxt = "n")
  plot(c(max(pos),min(pos)), c(val_min,val_max1), type="n", xlab=xlab, ylab=ylab,
       axes=F, main=main, cex.lab=1.5)
  
  for (i in 1:chrom_un) {
    end = chrom_cumsum[i]
    start = chrom_cumsum[i] - chrom_table[i] + 1

    x = pos[start:end]
    y = val[start:end]

    inds = y>val_thres
    # inds1 = y>val_max1
    # inds2 = y>val_max2 & y<val_max1

    points(x, y, col=colour[i], pch=16, cex=1)
    points(x[inds], y[inds], col=colour[i], cex=2)
    # points(x[inds2], y[inds2], col=colour[i], pch="x", cex=0.5)
    # points(x[inds1], y[inds1], col=colour[i], pch=20)
  }

  par(xaxt="s", yaxt="s")
  if (length(val)<label.x) {
    axis(side=1, at=1:length(val), labels=names(val))
  } else {
    axis(side=1, at=c(0, pos[round(chrom_cumsum)], max(pos)), F)
  }
  text(mids, par("usr")[3] - .5, srt=0, pos=2, cex=1.1, offset=-0.2,
       labels=chrom_unique[1:chrom_un], xpd=T)
  axis(side=2, at=guidelines)

  for (i in guidelines) abline(h=i, col="grey", lty="dotted")

  
  if (draw_line) {
    abline(h=(val_max2), col="black", lty="dotted")
    abline(h=(val_thres), col="red", lty="dotted")
    
    if (!is.null(lines)) {
      for (lin in lines) {
        abline(h=lin, col="blue", lty="dotted")
      }
    }
  }
}



#From http://bioinfo-mite.crb.wsu.edu/Rcode/wgplot.R
#See also https://stat.ethz.ch/pipermail/r-help/2008-November/180812.html
###############################################################################
###
### Whole Genome Significance plot 
### Matt Settles
### Bioinformatics Core
### Washington State University, Pullman, WA
### 
### Created July 7, 2008
###
### July 8, 2008 - fixed color goof
###############################################################################
##############
### things to add
### marker name on plot for significant markers
##############

### THERE ARE ERRORS IN GAPS MHTPLOT, SO THIS IS A FIX
## data 	a data frame with three columns representing chromosome, position and p values logged or unlogged
## logscale a flag to indicate if p value are to be log-transformed, FALSE means already logtransformed
## base 	the base of the logarithm, when logscale =TRUE
## guidelines 	the cutt-offs where horizontal line(s) are drawn
## color 	the color for different chromosome(s), and random if unspecified
## chrom_unique 	chrom_unique for the x-axis, length = number of chromosomes
## xlab   label to be placed on the X axis
## ylab   lable to be placed on the Y axis
## ... 	other options in compatible with the R plot function

## USAGE
# source("http://bioinfo-mite.crb.wsu.edu/Rcode/wgplot.R")
## fake example with Affy500k data
# affy =c(40220, 41400, 33801, 32334, 32056, 31470, 25835, 27457, 22864, 28501, 26273, 
#          24954, 19188, 15721, 14356, 15309, 11281, 14881, 6399, 12400, 7125, 6207)
# chrom_cumsum = cumsum(affy)
# n.markers = sum(affy)
# chrom_un = length(affy)
# test = data.frame(chr=rep(1:chrom_un,affy),pos=1:n.markers,p=runif(n.markers))
# png("wgplot.png",units="in",width=8,height=5,res=300)
# par(las="2",cex=0.6,pch=21,bg="white")
# wgplot(test,guidelines = c(1,3, 5, 7, 9),color=palette()[2:5],chrom_unique=as.character(1:22))
# title("Whole Genome Associaton Plot of Significance for Chromosomes 1 to 22")
# dev.off()
##
# "wgplot" = function (
#   data, 
#   logscale = TRUE, 
#   base = 10, 
#   guidelines = c(3, 5, 7, 9),
#   siglines = NULL,
#   sigcolors = "red", 
#   color = sample(colors(), 26),
#   chrom = as.character(c(1:22,"X","Y","XY","MT")),
#   startbp = NULL,
#   endbp = NULL,
#   chrom_unique = as.character(c(1:22,"X","Y","XY","MT")),
#   xlab = "Chromosome",
#   ylab = "-Log10(p-value)", ...) {
#   if (any(is.na(data)))
#     data = data[-unique(which(is.na(data))%%nrow(data)),]
#   keep = which(data[,1] %in% chrom)
#   data = data[keep,]
#   if (!is.null(startbp) & !is.null(endbp) & length(chrom) == 1){
#     keep = which(data[,2] >= startbp & data[,2] <= endbp) 
#     data = data[keep,]       
#   }
#   
#   
#   chr  = data[, 1]
#   pos  = data[, 2]
#   p    = data[, 3]
#   
#   ### remove any NAs
#   which(is.na(data[,2]))
#   chr  = replace(chr,which(chr == "X"),"100")
#   chr  = replace(chr,which(chr == "Y"),"101")
#   chr  = replace(chr,which(chr == "XY"),"102")
#   chr  = replace(chr,which(chr == "MT"),"103")	
#   
#   ord  = order(as.numeric(chr),as.numeric(pos))
#   chr  = chr[ord]
#   pos  = pos[ord]
#   p    = p[ord]
#   
#   chrom_table = as.vector(table(as.numeric(chr)))
#   chrom_cumsum = cumsum(chrom_table)
#   n.markers = sum(chrom_table)
#   chrom_un = length(chrom_table)
#   id = 1:chrom_un
#   color = rep(color,ceiling(chrom_un/length(color)))
#   if (logscale)
#     p = -log(p,base)        
#   if ( any(diff(pos) < 0) ) {
#     pos_cum =  cumsum(c(0,pos[which(!duplicated(chr))-1]))
#     pos = pos + rep(pos_cum,chrom_table)
#     
#     mids = pos_cum + diff(c(pos_cum,max(pos)))/2
#   }
#   
#   par(xaxt = "n", yaxt = "n")
#   plot(c(pos,pos[1]), c(9,p), type = "n", xlab = xlab, ylab = ylab, axes = FALSE,  ...)
#   for (i in 1:chrom_un) {
#     u = chrom_cumsum[i]
#     l = chrom_cumsum[i] - chrom_table[i] + 1
#     cat("Plotting points ", l, "-", u, "\n")
#     points(pos[l:u], p[l:u], col = color[i], ...)
#   }
#   par(xaxt = "s", yaxt = "s")
#   axis(1, at = c(0, pos[round(chrom_cumsum)],max(pos)),FALSE)
#   text(mids, par("usr")[3] - 0.5, srt = 0, pos=2,cex=0.5,offset= -0.2,
#        chrom_unique = chrom_unique[1:chrom_un], xpd = TRUE)
#   #axis(side=1, at =  pos[round(chrom_cumsum-chrom_table/2)],tick=FALSE, chrom_unique= chrom_unique[1:chrom_un])
#   #abline(h = guidelines)
#   axis(side=2, at = guidelines )
#   if (!is.null(siglines))
#     abline(h = -log(siglines,base),col=sigcolors)
#   
#   #mtext(eval(expression(guidelines)), 2, at = guidelines)
#   
# }


## by amrit
# sizeStripFont	font of size of facet labels
# xAngle	angle of x-axis labels
# hjust	horizontal justification 0-left, 0.5-center, 1-right
# vjust	vertical justification 0-low, 0.5-middle, 1-high
# xSize	font size of x-axis label
# ySize	font size of y-axis label
# xAxisSize	font size of x-axis label title
# yAxisSize	fotn size of y-axis label title
# weights	are the predicted scores/probablities of test data
# trubeLabels	are the true labels associated with the test data
# direction	= "auto", ">", "<"
customTheme = function (sizeStripFont, xAngle, hjust, vjust, xSize, ySize, 
                        xAxisSize, yAxisSize) {
  theme(strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1), strip.text.x = element_text(size = sizeStripFont), 
        strip.text.y = element_text(size = sizeStripFont), axis.text.x = element_text(angle = xAngle, 
                                                                                      hjust = hjust, vjust = vjust, size = xSize, color = "black"), 
        axis.text.y = element_text(size = ySize, color = "black"), 
        axis.title.x = element_text(size = xAxisSize, color = "black"), 
        axis.title.y = element_text(size = yAxisSize, color = "black"), 
        panel.background = element_rect(fill = "white", color = "black"))
}




## by amrit
## input: demo (file x col matrix); groups=class_col; variables=cols;
## output:
descriptiveStat = function(demo, groups, variables, paired=F, pairing=F) {
  require(dplyr)
  require(tidyr)
  require(broom)
  
  if (all(paired)) {
    X <- demo[, c(variables, groups, pairing), drop = FALSE]
    colnames(X) <- c(variables, "Group", "Pairing")
    lvls <- levels(X$Group)
    meanSD <- X %>% gather(Variable, Value, -c(Group, Pairing)) %>% 
      dplyr::group_by(Variable, Group) %>% dplyr::summarise(MEAN = mean(Value, 
                                                                        na.rm = TRUE), SD = sd(Value, na.rm = TRUE))
    pval0 <- X %>% gather(Variable, Value, -c(Group, Pairing)) %>% 
      dplyr::group_by(Variable) %>% nest() %>% dplyr::mutate(model = purrr::map(data, 
                                                                                ~tryCatch(lme(Value ~ Group, random = ~1 | Pairing, 
                                                                                              data = .), error = function(e) NA)))
    pval <- do.call(rbind, lapply(pval0$model, function(i) {
      tryCatch(summary(i)$tTable[2, ], error = function(e) NA)
    })) %>% data.frame %>% mutate(Variable = variables, term = paste("Group", 
                                                                     lvls[2]), BH.FDR = p.adjust(p.value, "BH"))
  } else {
    X <- demo[, c(variables, groups), drop = FALSE]
    colnames(X) <- c(variables, "Group")
    lvls <- levels(X$Group)
    meanSD <- X %>% gather(Variable, Value, -Group) %>% dplyr::group_by(Variable, 
                                                                        Group) %>% dplyr::summarise(MEAN = mean(Value, na.rm = TRUE), 
                                                                                                    SD = sd(Value, na.rm = TRUE))
    pval <- X %>% gather(Variable, Value, -Group) %>% dplyr::group_by(Variable) %>% 
      nest() %>% dplyr::mutate(model = purrr::map(data, 
                                                  ~lm(Value ~ Group, data = .))) %>% unnest(model %>% 
                                                                                              purrr::map(broom::tidy)) %>% group_by(Variable) %>% 
      slice(2)
    pval$BH.FDR <- p.adjust(pval$p.value, "BH")
  }
  return(list(meanSD = meanSD, pval = pval))
}

## by amrit
normalizelibSum = function(genExp) {
  lib.size <- colSums(genExp)
  genExpNorm <- t(log2(t(genExp + 0.5)/(lib.size + 1) * 1e+06))
  return(genExpNorm)
}

## by amrit
annotate_transcripts = function (features, filter, mart, trinity_map_dir, attr = c("description", "ucsc", "chromosome_name", "strand", 
                                                                 "hgnc_symbol", "refseq_mrna")) {
  
  if (filter %in% c("ucsc", "trinity")) {
    features = features
  }
  if (filter == "ensembl_gene_id") {
    features <- unlist(lapply(strsplit(features, "\\."), 
                              function(i) i[1]))
  }
  gene <- rep(NA, length(features))
  if (filter %in% c("ucsc", "ensembl_gene_id")) {
    hk.known <- getBM(attributes = attr, filters = filter, 
                      values = features, mart = mart)$hgnc_symbol
    gene <- unique(hk.known)
  }
  else {
    trinityMapFile <- read.delim(trinity_map_dir)
    trinityMapFile$Contig <- unlist(lapply(strsplit(as.character(trinityMapFile$query_id), 
                                                    "_"), function(i) paste(i[1], i[2], sep = "_")))
    trinityMapFile$UniProt <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(trinityMapFile$subject_id), 
                                                                            "\\|"), function(i) i[[2]])), split = "_"), function(x) x[1]))
    trinityMapFile$GenSym <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(trinityMapFile$subject_id), 
                                                                           "\\|"), function(i) i[[3]])), split = "_"), function(x) x[1]))
    gene <- trinityMapFile$GenSym[trinityMapFile$query_id %in% 
                                    features]
  }
  gene
}

## by casey
## input: vector of gene symbols
## output: 
sear = function (input, type = c("mrna", "mirna")) {
  data("collections", envir = environment())
  type <- match.arg(type)
  tbl <- switch(type, 
                mrna = dplyr::select(collections, collection:geneset, members = members_mrna), 
                mirna = dplyr::select(collections, collection:geneset, members = members_mirna))
  uni <- tbl$members %>% unlist() %>% unique()
  recognized <- input[input %in% uni]
  if (length(recognized) < 10) {
    warning(sprintf("Submitted %s symbols, but only %s are recognized.", 
                    length(input), length(recognized)))
  }
  input <- recognized
  tbl %>% dplyr::rowwise(.) %>% dplyr::mutate(n_input = length(input), 
                                              n_geneset = length(members), intersect = length(intersect(input, 
                                                                                                        members)), p_value = phyper(intersect - 1, n_geneset, 
                                                                                                                                    length(uni) - n_geneset, n_input, lower.tail = F)) %>% 
    dplyr::ungroup(.) %>% dplyr::group_by(collection) %>% 
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
    dplyr::ungroup(.)
}

## by amrit
plotSampleHist = function (data = data, main = NULL, xlim = NULL, ylim = NULL) {
  for (i in 1:ncol(data)) {
    idx <- data[, i] > -1
    shist(data[idx, i], unit = 0.25, col = i, plotHist = FALSE, 
          add = i != 1, main = main, ylim = ylim, xlim = xlim, 
          xlab = expression("log"[2] ~ "cpm"))
  }
}