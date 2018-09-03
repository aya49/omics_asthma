################################################
#
# visualizationFunctions.R
# Author: Amrit Singh
# Date: April 01, 2016
#
# functions
#     1) plotIndiv_diablo; 2 sub-functions; splotMatPlot() and panel.ellipses
#     2) heatmap_diablo
#     3) enrichPathwayNetwork_diablo
#     4) circosPlot_diablo
#
################################################
library(RColorBrewer)
#col.mixo <- brewer.pal(n = 8, name = "PuOr")
col.mixo <- c("goldenrod2", "aquamarine3", "salmon", "black")
#col.mixo <- c("#7570B3", "#E7298A", "#C2C2C2", "black",   "#D55E00", "#0072B2", "#999999")
#plot(1:12, pch = 19, col=c(col.mixo1, col.mixo))

###---------------------------------------------------------------------------------------
################################################
#
## 1) plotIndiv_diablo
#
################################################
plotIndiv_diablo = function(object, ncomp = 1, groupOrder = NULL){
  
  VarX <- do.call(cbind, lapply(object$variates, function(i) i[, ncomp]))
  datNames <- colnames(VarX)
  Y <- object$Y
  if(is.null(groupOrder)){
    groupOrder = levels(Y)
  } else {
    groupOrder
  }
  
  if (!is.factor(Y))
    stop(gettextf("Y must be a factor!"))
  if (length(ncomp) != 1)
    stop(gettextf("You can only choose one component"))
  
  numberOfCols <- ncol(VarX)
  numberOfRows <- numberOfCols - 1
  
  mat <- matrix(0, nrow = numberOfRows, ncol = numberOfRows)
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      mat[i,j] <- paste(i,j, sep="_")
    }
  }
  plotType = list(cor=mat[lower.tri(mat)], scatter=mat[upper.tri(mat)],
                  lab=diag(mat), 
                  bar=paste(1:(numberOfRows-1), numberOfCols, sep="_"),
                  stackedbar=paste(numberOfRows, numberOfCols, sep="_"))
  
  par(mfrow = c(numberOfRows, numberOfCols), mar = c(0.5, 0.5, 0.5, 1.5), oma = c(2,2,2,2))
  for(i in 1:numberOfRows){
    for(j in 1:numberOfCols){
      ptype <- unlist(lapply(plotType, function(x){
        intersect(paste(i,j,sep="_"), x)
      }))
      splotMatPlot(x=VarX[, i], y=VarX[, j], datNames, Y, ptype, groupOrder)
      if(i == 1 & j %in% seq(2, numberOfRows, 1)){Axis(side = 3, x=VarX[, i])}
      if(j == numberOfRows & i %in% seq(1, numberOfRows-1, 1)){Axis(side = 4, x=VarX[, i])}
    }
  }
}

splotMatPlot = function(x, y, datNames, Y, ptype, groupOrder){
  if(names(ptype) == "cor"){
    plot(1, type = "n", axes = FALSE)
    r = round(cor(x, y), 2)
    text(1, 1, labels=r, cex = abs(0.6/strwidth(r)*r))
    box()
  }
  if(names(ptype) == "scatter"){
    panel.ellipses(x=x, y=y, Y = Y)
  }
  if(names(ptype) == "lab"){
    plot(1, type = "n", axes = FALSE)
    ind = as.numeric(unlist(lapply(strsplit(ptype, "_"), unique)))
    text(x=1, y=1, labels=datNames[ind], cex = 2)
    box()
  }
  if(names(ptype) == "bar"){
    Y2 <- factor(as.character(Y), levels = groupOrder)
    par(las=2)
    boxplot(x ~ Y2, horizontal=TRUE, axes = FALSE, ylim = c(min(x)-3, max(x)),
            col= col.mixo[match(levels(Y2), levels(Y))])
    axis(4, at=1:nlevels(Y2), labels=levels(Y2))
  }
  if(names(ptype) == "stackedbar"){
    Y2 <- factor(as.character(Y), levels = groupOrder)
    bars <- table(Y2)
    par(las=1)
    barplot(bars, col= col.mixo[match(levels(Y2), levels(Y))], 
            axes = FALSE)
    axis(4, at=seq(0,max(bars),length.out=5), labels=seq(0,max(bars),length.out=5))
  }
}
  

panel.ellipses = function(x, y, Y = Y, pch = par("pch"), col.lm = "red", axes = FALSE,
                           ...) {
  ind.gp = matrice = cdg = variance = list()
  for(i in 1:nlevels(Y)){
    ind.gp[[i]] = which(as.numeric(Y)==i)
  }
  
  matrice = lapply(ind.gp, function(z){matrix(c(x[z], y[z]), ncol = 2)})
  cdg = lapply(matrice, colMeans)
  variance = lapply(matrice, var) 
  
  library(ellipse)
  coord.ellipse = lapply(1:nlevels(Y), function(x){ellipse(variance[[x]], centre = cdg[[x]], ellipse.level = 0.95)})
  max.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, max)})
  min.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, min)})
  ind.names <- names(Y)
  cex = 0.5
  
  plot(x, y, xlab = "X.label", ylab = "Y.label", col=col.mixo[as.numeric(Y)], pch=20, axes=axes,
       xlim = c(min(x, min.ellipse[1, ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, ]), max(y, max.ellipse[2, ])))
  #text(x, y, ind.names, col = col, cex = cex)  
  box()
  for (z in 1:nlevels(Y)){
    points(coord.ellipse[[z]], type = "l", col = col.mixo[c(1:nlevels(Y))[z]])
  }
}
###---------------------------------------------------------------------------------------

################################################
#
## 2) heatmap_diablo
#
################################################
heatmap_diablo = function(object, Y, legend1.x.y = c(0.75, 0.8), 
                          legend2.x.y=c(0.75, 0.55), margins = c(2, 15)){
  X <- object$X
  Y <- object$Y
  
  keepA = lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)
  XDatList <- mapply(function(x, y){
    x[, y]
  }, x=X, y=keepA[-length(keepA)], SIMPLIFY = FALSE)
  XDat <- do.call(cbind, XDatList)
  XDat[which(XDat > 2)] <- 2
  XDat[which(XDat < -2)] <- -2

  dark <- brewer.pal(n = 12, name = 'Paired')[seq(2, 12, by = 2)]
  VarLabels <- factor(rep(names(X), lapply(keepA[-length(keepA)], sum)), levels = names(X)[order(names(X))])
  
  ## Plot heatmap
  cim(XDat, row.names = rep("", nrow(XDat)), col.names = rep("", ncol(XDat)),
      col.sideColors = dark[as.numeric(VarLabels)],  clust.method = c("average", "average"),
      row.sideColors = col.mixo[as.numeric(Y)], margins = margins)
  legend(x=legend1.x.y[1], y=legend1.x.y[2], levels(Y), 
         col=col.mixo[1:nlevels(Y)], pch = 19, bty="n", title = "Rows")
  legend(x=legend2.x.y[1], y=legend2.x.y[2], names(X), 
         col=dark[1:nlevels(VarLabels)][match(levels(VarLabels), names(X))], pch = 19, bty="n", title = "Columns")
  
  }

cimDiablo2 = function (object, margins = c(2, 15), pos.legend = "topright",
cex.legend = 1.5)
{
    if (!any(class(object) == "block.splsda"))
    stop("heatmapDiablo is only available for 'block.splsda' objects")
    if (length(object$X) <= 1)
    stop("This function is only available when there are more than 3 blocks")
    if (ncomp > min(object$ncomp))
    stop("'ncomp' needs to be higher than object$ncomp")
    X <- object$X
    Y <- object$Y
    keepA <- lapply(object$loadings, function(x)
      apply(x, 2, function(i) names(i)[which(i != 0)]))
    XDatList <- mapply(function(x, y) {
      x[, y]
    }, x = X, y = keepA[-length(keepA)], SIMPLIFY = FALSE)

    XDat <- do.call(cbind, XDatList)
    XDat[which(XDat > 2)] <- 2
    XDat[which(XDat < -2)] <- -2
    dark <- brewer.pal(n = 12, name = "Paired")[seq(2, 12, by = 2)]
    col = rep(1:length(X), unlist(lapply(keepA[-length(keepA)], length)))+3
    
    #VarLabels <- factor(rep(names(X), lapply(keepA[-length(keepA)],
    #sum)), levels = names(X)[order(names(X))])
    opar <- par()[!names(par()) %in% c("cin", "cra", "csi", "cxy",
    "din", "page")]
    par(mfrow = c(1, 1))
    cim(cor(XDat), col.names = rep("", ncol(XDat)), row.sideColors = color.mixo(col),
    col.sideColors = color.mixo(col), margins = margins)
    legend("topright", names(X), col = color.mixo(4:6), pch = 19, bty = "n")
    par(opar)
}

################################################
#
## 3) enrichPathwayNetwork_diablo
#
################################################
enrichPathwayNetwork_diablo = function(enrich, cutoff = 0.25, trim = TRUE){
  library(networkD3)
  jaccard=function (s1, s2){
    length(intersect(s1, s2))/length(union(s1, s2))
  }
  
  rebase_links = function (nodes){
    ref <- nodes$members
    names(ref) <- nodes$rowid
    t(combn(as.numeric(names(ref)), 2)) %>% as.data.frame(.) %>% 
      dplyr::tbl_df(.) %>% dplyr::rename(source = V1, target = V2) %>% 
      dplyr::mutate(jaccard = unlist(purrr::map2(ref[as.character(source)], 
                                                 ref[as.character(target)], jaccard))) %>% dplyr::mutate(source = match(source, 
                                                                                                                        nodes$rowid) - 1, target = match(target, nodes$rowid) - 
                                                                                                           1)
  }

  ## make nodes
  nodes <- enrich %>% dplyr::add_rownames("rowid") %>% 
    dplyr::mutate(rowid = as.numeric(rowid) - 1, group = -log10(fdr +  .Machine$double.xmin))
  links <- rebase_links(nodes) %>% filter(jaccard >= cutoff)
  if (trim) {
    selection <- c(links$source, links$target) %>% unique() %>% 
      match(., (1:nrow(nodes) - 1), nomatch = NA) %>% sort()
    nodes <- nodes %>% dplyr::slice(selection) %>% dplyr::mutate(rowid = (0:(n() - 
                                                                               1)))
    links <- rebase_links(nodes) %>% dplyr::filter(jaccard >= 
                                                     cutoff)
  }
  
  create_colorscale = function (nodes, palette) 
{
  cols <- RColorBrewer::brewer.pal(9, palette)[-c(1:2)]
  cols <- paste0("'", paste(cols, collapse = "', '"), "'")
  max <- -log10(nodes$fdr) %>% max() %>% ceiling()
  range <- quantile(0:max) %>% ceiling()
  range <- paste(range, collapse = ", ")
  networkD3::JS(paste0("d3.scale.linear().domain([", range, 
                       "]).range([", cols, "])"))
}
  
  # Create graph with node text faintly visible when no hovering
  #forceNetwork(Links = links, Nodes = nodes, Source = "source", Nodesize = "intersect",
  #             Target = "target", Value = "jaccard", NodeID = "rnaet",
  #             Group = "subcollection", bounded = FALSE,
  #             opacityNoHover = TRUE, 
  #             fontSize=12, zoom=T, legend=T,
  #             opacity = 0.8, charge=-300)
  
  networkD3::forceNetwork(Links = links, Nodes = nodes, NodeID = "geneset", 
                          Nodesize = "n_rnaset", Group = "group", Source = "source", 
                          Target = "target", Value = "jaccard", linkDistance = networkD3::JS("function(d) { return d.value * 100; }"), 
                          colourScale = create_colorscale(nodes, "BuPu"), fontSize = 16, 
                          fontFamily = "sans-serif", opacity = 0.85, zoom = F, 
                          legend = T, bounded = T, opacityNoHover = 1, charge = -300)
}

################################################
#
# 4) circosPlot
#
################################################
circosPlot_diablo = function(object, corThreshold,cex.label,
                             figSize = 800,
                             segmentWidth = 25,
                             linePlotWidth = 90,
                             showIntraLinks = FALSE)
{
  options(stringsAsFactors = FALSE);
  set.seed(1234);
  
  ##############################
  ###   networkDiagram_core.R
  ###
  ###   Authors: Michael Vacher (minor changes by Amrit :)
  ###
  ###   Parts of this src has been modified from the original OmicCircos package obtained from: 
  ###   Ying Hu Chunhua Yan <yanch@mail.nih.gov> (2015). OmicCircos: High-quality circular visualization of omics data. R package version 1.6.0.
  ############################## 
  
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  
  par(mfrow  = c(1,1), mar = rep(1, 4))
  X <- object$X
  for(i in 1 : length(X)){
    colnames(X[[i]]) <- paste(names(X)[i], colnames(X[[i]]), sep = "_")
  }
  Y <- object$Y

  keepA = lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)
  cord = mapply(function(x, y, keep){
    cor(x[, keep], y, use = "pairwise")
  }, x=X, y=object$variates[-length(object$variates)], keep = keepA[-length(keepA)], SIMPLIFY = FALSE)
  
  simMatList <- vector("list", length(X))
  for(i in 1:length(cord)){
    for(j in 1:length(cord)){
      simMatList[[i]][[j]] <- cord[[i]] %*% t(cord[[j]])
    }
  }
  corMat <- do.call(rbind, lapply(simMatList, function(i) do.call(cbind, i)))
  
  ## Expression levels
  Xdat <- as.data.frame(do.call(cbind, X)[, colnames(corMat)])
  
  AvgFeatExp0 <- Xdat %>% mutate(Y = Y) %>% gather(Features, Exp, -Y) %>%
    group_by(Y, Features) %>% dplyr::summarise(Mean = mean(Exp), SD = sd(Exp))
  AvgFeatExp0$Dataset <- unlist(lapply(strsplit(AvgFeatExp0$Features, "_"), function(i) i[1]))
  featExp <- AvgFeatExp0 %>% group_by(Dataset, Y) %>% arrange(Mean)
  #featExp$Features <- unlist(lapply(strsplit(as.character(featExp$Features), "_"), function(i) i[2]))
    # Generate a circular plot (circos like) from a correlation matrix (pairwise)
    # 
    # Args:
    #   corMat: the main correlation matrix.
    #         -> colnames == rownames (pairwise); values = correlations
    #   featExp: data.frame holding the expression data.
    #   corThreshold: minimum value for correlations (<threshold will be ignored)
    #   figSize: figure size
    #   segmentWidth: thickness of the segment (main circle)
    #   linePlotWidth: thickness of the line plot (showing expression data)
    #   showIntraLinks = display links intra segments
    
    # 1) Generate karyotype data
    chr <- genChr(featExp);
    chr.names <- unique(chr$chrom);# paste("chr", 1:seg.num, sep="");
    # Calculate angles and band positions
    db <- segAnglePo(chr, seg=chr.names);
    db <- data.frame(db);
    
    # 2) Generate Links
    links <- genLinks(chr, corMat, corThreshold);
    if(nrow(links) < 1)
      stop("Choose a lower correlation threshold")

      # 3) Plot
      # Calculate parameters
      circleR <- (figSize / 2.0) -  segmentWidth - linePlotWidth;
      linksR <- circleR - segmentWidth;
      linePlotR <- circleR + segmentWidth
      chrLabelsR <- (figSize / 2.0)  ;
      
      par(mar=c(2, 2, 2, 2));
      
      plot(c(1,figSize), c(1,figSize), type="n", axes=FALSE, xlab="", ylab="", main="");
      
      # Plot ideogram
      drawIdeogram(cex.label=cex.label, R=circleR, cir=db, W=segmentWidth,  show.band.labels=TRUE, show.chr.labels=TRUE, chr.labels.R= chrLabelsR, chrData=chr);
      
      # Plot links
      drawLinks(R=linksR, cir=db,   mapping=links,   col=linkColors, drawIntraChr=showIntraLinks);
      
      # Plot expression values
      #cTypes <- unique(featExp[,1]) #Get the different disease/cancer types (lines)
      #lineCols <- rainbow(nrow(cTypes), alpha=0.5);          ## removed by Amrit
      cTypes <- levels(Y)                    ## added by amrit
      lineCols <- col.mixo[1:nlevels(Y)]         ## added by amrit
      
      # Fixme: remove this loop and send the whole expr dframe to drawLinePlot
      for (i in 1:length(chr.names)){
        seg.name <- gsub("chr","",chr.names[i]);
        #Get data for each segment
        expr <- subset(featExp,featExp$Dataset==seg.name)
        
        expr <- dcast(expr, formula = Features ~ Y, value.var="Mean")   ## changed PAM50 to Y
        expr <- merge(expr, chr, by.x="Features", by.y="name");
        expr$po <- (as.numeric(expr$chromStart) + as.numeric(expr$chromEnd)) / 2.0;
        expr <- dplyr::rename(expr, seg.name = chrom, seg.po = po);
        
        # Reorder columns
        cOrder <- c(c(grep("seg.name", colnames(expr)),
                      grep("seg.po", colnames(expr))), c(1:length(cTypes)+1))     ## changed 1:nrow to 1:length  removed c(1:length(cTypes)+1)
        expr <- expr[, cOrder];
        
        # Plot data on each sub segment
        subChr <- subset(db, db$seg.name == chr.names[i] )
        drawLinePlot(R=linePlotR, cir=subChr,   W=linePlotWidth, lineWidth=1, mapping=expr, col=lineCols, scale=FALSE);
      }
      
      # Plot legend
      # First legeng bottom left corner
      legend(x=5, y = (circleR/4), title="Correlations", c("Positive Correlation", "Negative Correlation"), 
             col = c(colors()[134], colors()[128]), pch = 19, cex=0.8, bty = "n")
      # Second legend bottom righ corner
      legend(x=figSize-(circleR/4), y = (circleR/3), title="Expression", legend=cTypes,  ## changed PAM50 to Y, changed by Amrit
             col = lineCols, pch = 19, cex=0.8, bty = "n")
      # third legend top left corner
      legend(x=figSize-(circleR/2), y = figSize, title="Correlation cut-off", legend=paste("r", corThreshold, sep = "="),
             col = "black", cex=0.8, bty = "n") 
}

drawIdeogram <- function(R, xc=400, yc=400, cir, W, cex.label,
                         show.band.labels = FALSE,  
                         show.chr.labels = FALSE, chr.labels.R = 0,
                         chrData)
{
  # Draw the main circular plot: segments, bands and labels
  chr.po    <- cir;
  chr.po[,1]  <- gsub("chr","",chr.po[,1]);
  chr.num     <- nrow(chr.po);
  
  dat.c     <- chrData;
  dat.c[,1] <- gsub("chr", "", dat.c[,1]);
  
  for (chr.i in c(1:chr.num)){
    chr.s  <- chr.po[chr.i,1];
    
    v1 <- as.numeric(chr.po[chr.i,2]);
    v2 <- as.numeric(chr.po[chr.i,3]);
    v3 <- as.numeric(chr.po[chr.i,6]);
    v4 <- as.numeric(chr.po[chr.i,7]);
    
    dat.v <- subset(dat.c, dat.c[,1]==chr.s);
    dat.v <- dat.v[order(as.numeric(dat.v[,2])),];
    for (i in 1:nrow(dat.v)){  
      
      #col.v <- which(colors()==dat.v[i,5]); #get color index
      #col <- colors()[col.v]
      dark.clear <- brewer.pal(n = 12, name = 'Paired')
      col.v <- which(dark.clear==dat.v[i,5]); #get color index
      col <- dark.clear[col.v]
      
      w1 <- scale.v(as.numeric(dat.v[i,2]), v1, v2, v3, v4);
      w2 <- scale.v(as.numeric(dat.v[i,3]), v1, v2, v3, v4);
      
      draw.arc.s(xc, yc, R, w1, w2, col=col, lwd=W);
      
      if (show.band.labels){
        band.text <- unlist(lapply(strsplit(as.character(dat.v[i,4]), "_"), function(i) paste(i[-1], collapse = "_")));    ## changed by amrit
        
        band.po <- ((w1+w2)/2)# - ((w2-w1)/3) #position around the circle
        # print(c(band.po, w1, w2, (w2-w1)/3))
        band.po.in <- R-(W/3.0) #position on the band (middle)
        draw.text.rt(xc, yc,band.po.in  , band.po , band.text , cex=cex.label, segmentWidth = W, side="in" );
      } 
    } #End for row
    if (show.chr.labels){
      w.m <- (v1+v2)/2;   
      chr.t <- gsub("chr", "", chr.s);
      draw.text.rt(xc, yc, chr.labels.R, w.m, chr.t, cex=1, segmentWidth = W, parallel=TRUE);
    }
  } #End for
}

drawLinks <- function(R, xc=400, yc=400, cir, W,  
                      mapping=mapping,
                      lineWidth=1, col=rainbow(10, alpha=0.8)[7],  drawIntraChr=FALSE) 
{
  # Draw the links (computed correlation) between features
  chr.po    <- cir;
  chr.po[,1]  <- gsub("chr","",chr.po[,1]);
  chr.num     <- nrow(chr.po);
  
  chr.po[,4] <- gsub("chr", "", chr.po[,4]);
  dat.in <- mapping;
  dat.in[,1] <- gsub("chr", "", dat.in[,1]);
  dat.in[,4] <- gsub("chr", "", dat.in[,4]);
  
  dat    <- dat.in;
  
  for (i in 1:nrow(dat)){
    chr1.s   <- dat[i,1];
    chr2.s   <- dat[i,4];
    po1      <- dat[i,2];
    po2      <- dat[i,5];
    
    chr1     <- which(chr.po[,1]==chr1.s);
    chr2     <- which(chr.po[,1]==chr2.s);
    
    v1 <- as.numeric(chr.po[chr1,2]);
    v2 <- as.numeric(chr.po[chr1,3]);
    v3 <- as.numeric(chr.po[chr1,6]);
    v4 <- as.numeric(chr.po[chr1,7]);
    
    w1 <- scale.v(as.numeric(po1), v1, v2, v3, v4);
    
    v1 <- as.numeric(chr.po[chr2,2]);
    v2 <- as.numeric(chr.po[chr2,3]);
    v3 <- as.numeric(chr.po[chr2,6]);
    v4 <- as.numeric(chr.po[chr2,7]);
    
    w2 <- scale.v(as.numeric(po2), v1, v2, v3, v4);
    # Set the link width depending on the correlation coefficient
    lwd <- abs(as.numeric(dat[i,7]));
    
    # Set link color
    if (as.numeric(dat[i,7]) < 0.0){
      linkCol <- colors()[128]; #pale red
    } else {
      linkCol <- colors()[134]; # blue
    }
    linkCol <- add.alpha(linkCol, alpha=0.4);
    
    if (chr1 == chr2){
      if (drawIntraChr == TRUE){
        draw.link(xc, yc, R, w1, w2, col=linkCol, lwd=lineWidth);
      }
    } else {
      draw.link(xc, yc, R, w1, w2, col=linkCol, lwd=lineWidth);
    }
  } ### End for
}

drawLinePlot <- function(mapping=mapping, xc=400, yc=400, col.v=3,
                         R, cir,   W, col='black', scale=FALSE, lineWidth=1,
                         background.lines=FALSE,axis.width=1)
{
  # Generate a linear plot around the main ideogram.
  #
  # fixme: the function writes the same line multiple times
  # it needs to be called only once and process the data/segment
  # separately
  chr.po    <- cir;
  chr.po[,1]  <- gsub("chr","",chr.po[,1]);
  chr.num     <- nrow(chr.po);
  
  dat.in   <- mapping;
  dat.in[,1] <- gsub("chr", "", dat.in[,1]);
  
  # data set for the chromosome
  for (chr.i in 1:chr.num){
    chr.s <- chr.po[chr.i,1];
    chr.s <- gsub("chr","",chr.s);      
    dat   <- subset(dat.in, dat.in[,1]==chr.s);
    dat   <- dat[order(as.numeric(dat[,2])),];
    v1 <- as.numeric(chr.po[chr.i,2]);
    v2 <- as.numeric(chr.po[chr.i,3]);
    v3 <- as.numeric(chr.po[chr.i,6]);
    v4 <- as.numeric(chr.po[chr.i,7]);
    
    # background line
    if (background.lines){
      draw.arc.pg(xc, yc, v1, v2, R, R+W-5, col=colors()[245]);
    } else {
      draw.arc.s(xc, yc, R, v1, v2, col=colors()[245], lwd=axis.width);
    }
  }
  
  my.R1 <- R + W/5;
  my.R2 <- R + W - W/5;
  
  ## for the matrix colors
  num.col <- ncol(dat[,col.v:ncol(dat)]);
  num.row <- nrow(dat.in);
  
  if (length(col) == num.col){
    colors <- col;
  } else {
    colors <- rainbow(num.col, alpha=0.5);
  }
  
  for (chr.i in 1:chr.num){
    chr.s <- chr.po[chr.i,1];
    chr.s <- gsub("chr","",chr.s);
    
    dat   <- subset(dat.in, dat.in[,1]==chr.s);
    dat   <- dat[order(as.numeric(dat[,2])),];
    #print(head(dat))
    dat.i   <- c(col.v:ncol(dat));
    dat.m   <- dat.in[,dat.i];
    dat.m   <- as.matrix(dat.m);
    dat.min <- min(as.numeric(dat.m), na.rm=T);
    dat.max <- max(as.numeric(dat.m), na.rm=T);
    
    v1 <- as.numeric(chr.po[chr.i,2]);
    v2 <- as.numeric(chr.po[chr.i,3]);
    v3 <- as.numeric(chr.po[chr.i,6]);
    v4 <- as.numeric(chr.po[chr.i,7]);
    
    col.i <- 0;
    
    for (j in col.v:ncol(dat)){
      col.i <- col.i + 1;
      col   <- colors[col.i];
      
      my.v      <- as.numeric(dat[1,j]); 
      dat.i.old <- my.v;
      v.old   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
      
      po      <- as.numeric(dat[1,2]);
      w.from  <- scale.v(po, v1, v2, v3, v4);
      
      
      for (k in 1:nrow(dat)){
        
        dat.i <- as.numeric(dat[k,j]);
        
        if (is.na(dat.i)){
          next;
        }
        
        v    <- scale.v(dat.i, my.R1, my.R2, dat.min, dat.max);
        w.to <- scale.v(as.numeric(dat[k,2]), v1, v2, v3, v4);
        
        if (w.from > 0){      
          draw.line3(xc, yc, w.from, w.to, v.old, v, col=col, lwd=lineWidth)
        } 
        
        dat.i.old <- dat.i;
        w.from    <- w.to;
        v.old     <- v;
      } # end the row
    }   # end the col
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
  }     # end the chr/segment
}

genChr <-function (expr, bandWidth = 1.0)
{
  # Generate the segments and calculate the 
  # unique positions of the bands
  #
  # Args:
  #   expr : dataframe containing the features expression
  # example: colnames(concatFeatExp) "PAM50"    "Features" "Mean"     "SD"       "Dataset" 
  #   bandWidth: thickness of each band 
  #
  # Return:
  #   a data.frame that can be used with segAnglePo
  
  # expr can contains expression data for multiple diseases
  # here, we only use the Chrom and Dataset column and remove duplicates
  keeps <- c("Features","Dataset");
  expr <- expr[keeps];
  expr <- unique(expr);
  chrLengths <- data.frame(table(expr$Dataset)) ;
  rownames(chrLengths) <- chrLengths[,1];
  chrLengths[,1] <- NULL;
  colnames(chrLengths) <- c( "Freq");
  chrLengths[, "Count"] <- rep(0, nrow(chrLengths));
  
  # Last column contains the bands' color
  # Create a color scheme
  #dark <- c("brown3","darkgoldenrod","antiquewhite3","steelblue3")
  #clear <- c("brown1","darkgoldenrod1","antiquewhite1","steelblue1")
  #dark.clear <- brewer.pal(n = 12, name = 'Paired')
  dark.clear <- c("#FB9A99", "#E31A1C","#FDBF6F", "#FF7F00", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4")
  dark <- dark.clear[seq(2, 12, by = 2)]
  clear <- dark.clear[seq(1, 12, by = 2)]
  chrColScheme <- data.frame(dark, clear)
  n_datasets = length(unique(expr$Dataset))
  chrColScheme <- chrColScheme[c(1:n_datasets),]
  rownames(chrColScheme) <- unique(expr$Dataset)
  
  seg.out <- c();
  for (i in 1:nrow(expr)){
    chrName    <- paste("chr", expr[i,'Dataset'], sep="");
    dType <- toString(expr[i,'Dataset']);
    pStart <- chrLengths[dType,'Count'] * bandWidth;
    pStop <- chrLengths[dType,'Count'] * bandWidth + bandWidth ;
    chrLengths[dType,'Count'] <- chrLengths[dType,'Count'] + 1;
    fName <- as.character(as.matrix(expr[i,'Features']));    # added as.character() by amrit
    # Assign colors
    if (chrLengths[dType,'Count'] %% 2 == 0){
      chrCol <- chrColScheme[dType,]$clear;
    } else{
      chrCol <- chrColScheme[dType,]$dark;
    }
    seg.out <- rbind(seg.out, c(chrName, pStart, pStop,  fName, chrCol));
  }
  
  # Use the same names than in omicCircos
  colnames(seg.out) <- c("chrom", "chromStart", "chromEnd", "name", "color");
  seg.out <- as.data.frame(seg.out);
  
  return(seg.out);
}

genLinks <- function(chr, corMat, threshold)
{
  # Generates the links corresponding to pairwise correlations
  #
  # Args:
  #   chr: ideogram structure (generated from genChr)
  #   corMat: main correlation matrix
  #   threshold: minimum correlation value
  #
  # Return:
  #   the links data (see omicsCircos doc)
  #   
  linkList <- c();
  # Remove matrix diagonal and the the mat in a list
  linkList <- subset(melt(corMat), Var1!=Var2 );
  # Remove links below the threshold
  linkList <- subset(linkList, abs(linkList$value) >= threshold);
  
  #First merge
  linkList <- dplyr::rename(linkList, feat1=Var1, feat2=Var2);  # CHANGED BY AMRIT
  linkList <- merge(linkList, chr, by.x="feat1", by.y="name");
  # Set the position in the middle of the band
  linkList$po1 <- (as.numeric(linkList$chromStart) + as.numeric(linkList$chromEnd)) / 2.0;
  linkList <- dplyr::rename(linkList, chr1=chrom);  # CHANGED BY AMRIT
  keeps <- c("feat1","feat2","value","chr1","po1");
  linkList <- linkList[keeps];
  
  #Second merge
  linkList <- merge(linkList, chr, by.x="feat2", by.y="name");
  linkList$po2 <- (as.numeric(linkList$chromStart) + as.numeric(linkList$chromEnd)) / 2.0;
  linkList <- dplyr::rename(linkList, chr2=chrom);  # CHANGED BY AMRIT
  keeps <- c("chr1","po1","feat1","chr2","po2","feat2","value");
  linkList <- linkList[keeps];
  
  return(linkList);
}

bezierCurve <- function(x, y, n=10)  {  
  outx <- NULL  
  outy <- NULL   
  i <- 1	
  for (t in seq(0, 1, length.out=n))		{		
    b <- bez(x, y, t)		
    outx[i] <- b$x		
    outy[i] <- b$y 		
    i <- i+1		
  } 	
  return (list(x=outx, y=outy))	
} 

##
bez <- function(x, y, t)	{	
  outx <- 0	
  outy <- 0	
  n <- length(x)-1	
  for (i in 0:n)		{		
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]		
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]		
  } 	
  return (list(x=outx, y=outy))	
}

###########################################
# one value : from a to b 
scale.v <- function(v, a, b, min.v, max.v) {
  v <- v-min.v; 
  v <- v/(max.v-min.v); 
  v <- v*(b-a);  
  v+a
}

### draw.link
draw.link <- function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
  # for translocation
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}

### draw.link2
draw.link2 <- function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
  # for translocation
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x2  <- xc+r/2*cos(w3);
  y2  <- yc-r/2*sin(w3);
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0, x2, x2, x1);
  y <- c(y0, y2, y2, y1);
  points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}
### 

###
draw.link.pg <- function(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd) {
  w1 <- w1.1;
  w2 <- w2.2;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc1 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w1.1-w1.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1.1,w1.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.1.x <- xc + cos(ang.seq) * r;
  fan.1.y <- yc - sin(ang.seq) * r;
  
  ######################################################
  w1 <- w1.2;
  w2 <- w2.1;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc2 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w2.1-w2.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w2.1,w2.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.2.x <- xc + cos(ang.seq) * r;
  fan.2.y <- yc - sin(ang.seq) * r;
  
  polygon(c(bc1$x, fan.2.x, rev(bc2$x), rev(fan.1.x)), 
          c(bc1$y, fan.2.y, rev(bc2$y), rev(fan.1.y)), 
          fillOddEven=TRUE, border=col, col=col, lwd=lwd); 
}

###
draw.point.w <- function(xc, yc, r, w, col=col, cex=cex){
  w <- w/360*2*pi;
  x <- xc+r*cos(w);
  y <- yc-r*sin(w);
  points(x, y, pch=20, col=col, cex=cex);
}

###
draw.text.w <- function(xc, yc, r, w, n, col="black", cex=1){
  w <- w%%360;
  w <- w/360*2*pi;
  x <- xc+r*cos(w);
  y <- yc-r*sin(w);
  text(x,y,labels=n, col=col, cex=cex);
}

###
draw.text.rt <- function(xc, yc, r, w, n, col="black", cex=1, side="out", segmentWidth=20, parallel=FALSE){
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- (segmentWidth*2)/2.0;
  b <- the.w
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
      if (parallel == TRUE){
        the.w <- the.w -90;# 180;
      }
    } else if (the.w > 90 & the.w <= 180) {
      if (parallel == TRUE){
        the.w <- the.w -90;# 180;
      }
      else {
        the.w <- the.w + 180;
      }
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      if (parallel == TRUE){
        the.w <- the.w -90;# 180;
      }
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.pos <- 4;
      if (parallel == TRUE){
        the.w <- the.w + 90;
      }
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  
  if (side=="in"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w + 180;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+segmentWidth;
    }
    if (the.pos==4){
      x <- x-segmentWidth;
    }
  }
  
  
  text(x, y, adj=0, offset=1, labels=n, srt=the.w, 
       pos=the.pos, col=col, cex=cex);
}

###strokeLine2
draw.line <- function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
  w  <- (w/360)*2*pi;
  x1 <- xc+l1*cos(w);
  y1 <- yc-l1*sin(w);
  x2 <- xc+l2*cos(w);
  y2 <- yc-l2*sin(w);
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend);
}

###strokeLine3
draw.line2 <- function (xc, yc, w, r, l, col=col, lwd=lwd){
  line_w   <- l;
  theangle <- w;
  l1       <- r;
  theangle <- (theangle/360)*2*pi;
  x0       <- xc+l1*cos(theangle);
  y0       <- yc+l1*sin(theangle);
  w1       <- 45/360*2*pi;
  x1 = xc + sin(w1) * (x0);
  y1 = yc + cos(w1) * (y0);
  x2 = xc - sin(w1) * (x0);
  y2 = yc - cos(w1) * (y0);
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt");
}

###strokeLine by two angles
draw.line3 <- function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd){
  theangle1 <- w1;
  theangle2 <- w2;
  l1        <- r1;
  l2        <- r2;
  
  theangle1 <- (theangle1/360)*2*pi;
  x1        <- xc+l1*cos(theangle1);
  y1        <- yc-l1*sin(theangle1);
  
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt");
}

### plot fan or sector that likes a piece of doughnut (plotFan)
draw.arc.pg <- function (xc, yc, 
                         w1, w2, r1, r2, col="lightblue", border="lightblue", lwd=0.01
){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r1;
  fan.i.y <- yc - sin(ang.seq) * r1;
  
  
  fan.o.x <- xc + cos(ang.seq) * r2;
  fan.o.y <- yc - sin(ang.seq) * r2;
  
  polygon(c(rev(fan.i.x), fan.o.x ), c(rev(fan.i.y), fan.o.y), 
          fillOddEven=TRUE, border=border, col=col, lwd=lwd, lend=1)
  
}

draw.arc.s <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  # Draw circular arcs for the main ideogram)
  # s = simple
  # r = radius
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  ## lend=0(round); lend=1(butt); lend=2(square)
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
  #points(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
}

#########################################
## segment to angle and position
## segAnglePo
#########################################

# get angle if given seg and position
# seg should be ordered by user
segAnglePo <- function (seg.dat=seg.dat, seg=seg, angle.start=angle.start, 
                        angle.end=angle.end){
  
  if (missing(angle.start)){
    angle.start <- 0;
  }
  if (missing(angle.end)){
    angle.end <- 360;
  }
  ## check data.frame?
  colnames(seg.dat) <- c("seg.name","seg.Start","seg.End","name","gieStain");
  
  ## get length of the segomosomes
  seg.l   <- c();
  seg.min <- c();
  seg.sum <- c();
  seg.s   <- 0;
  seg.num   <- length(seg);
  seg.names <- seg;
  
  ########################################################
  ########################################################
  for (i in 1:seg.num){
    seg.n <- seg.names[[i]];
    
    dat.m <- subset(seg.dat, seg.dat[,1]==seg.n);
    seg.full.l <- max(as.numeric(dat.m[,"seg.End"]));
    seg.full.m <- min(as.numeric(dat.m[,"seg.Start"]));
    seg.l      <- cbind(seg.l, seg.full.l);
    seg.min    <- cbind(seg.min, seg.full.m);
    seg.s      <- seg.s + seg.full.l;
    seg.sum    <- cbind(seg.sum, seg.s);
  }
  
  ## initial parameters
  gap.angle.size <- 2;
  seg.angle.from <- angle.start + 270;
  
  seg.full.l  <- sum(as.numeric(seg.l));
  angle.range <- angle.end - angle.start;
  cir.angle.r <- (angle.range - seg.num * gap.angle.size)/seg.full.l;
  
  out.s     <- c();
  l.old     <- 0;
  gap.angle <- 0;
  for (i in 1:seg.num){
    seg.n <- seg.names[[i]];
    dat.m <- subset(seg.dat, seg.dat[,1]==seg.n);
    len   <- seg.sum[i];
    w1    <- cir.angle.r*l.old + gap.angle;
    w2    <- cir.angle.r*len   + gap.angle;
    out.s     <- rbind(out.s, c(seg.n, w1+seg.angle.from, w2+seg.angle.from, l.old, len, seg.min[i], seg.l[i]));
    gap.angle <- gap.angle + gap.angle.size;
    l.old     <- len;
  }
  
  colnames(out.s) <- c("seg.name","angle.start", "angle.end", "seg.sum.start", "seg.sum.end","seg.start", "seg.end");
  return(out.s);
}


### start do.scale
do.scale <- function(xc=xc, yc=yc, dat.min=dat.min, dat.max=dat.max, 
                     R=R, W=W, s.n=1, col="blue"){
  dat.m   <- round((dat.min+dat.max)/2, s.n);
  dat.min <- round(dat.min, s.n);
  dat.max <- round(dat.max, s.n);
  y1      <- yc + R ;
  y2      <- yc + R + W/2;
  y3      <- yc + R + W;
  x1      <- xc - W/20;
  x2      <- x1 - (W/20)*1.2;
  x3      <- x1 - (W/20)*3;
  segments(x1, y1, x1, y3, lwd=0.01, col=col);
  segments(x1, y1, x2, y1, lwd=0.01, col=col);
  segments(x1, y2, x2, y2, lwd=0.01, col=col);
  segments(x1, y3, x2, y3, lwd=0.01, col=col);
  text(x3, y1, dat.min, cex=0.2, col=col);
  text(x3, y2, dat.m,   cex=0.2, col=col);
  text(x3, y3, dat.max, cex=0.2, col=col);
}
### end do.scale

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
###---------------------------------------------------------------------------------------







################################################
#
# 4) circosPlot
#
################################################
circosPlot_diablo2 = function(X, Y, multiOmicPanel, corThreshold,
                             figSize = 800, cex.label,
                             segmentWidth = 25,
                             linePlotWidth = 90,
                                line=TRUE,
                             showIntraLinks = FALSE)
{

  options(stringsAsFactors = FALSE);
  set.seed(1234);
  
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  
  par(mfrow  = c(1,1), mar = rep(1, 4))
  for(i in 1 : length(X)){
    colnames(X[[i]]) <- paste(names(X)[i], colnames(X[[i]]), sep = "_")
    multiOmicPanel[[i]] <- paste(names(X)[i], multiOmicPanel[[i]], sep = "_")
  }
  
  corMat <- cor(do.call(cbind, mapply(function(x, y){
    x[, y]
  }, x = X, y = multiOmicPanel, SIMPLIFY = FALSE)))
  
  ## Expression levels
  Xdat <- as.data.frame(do.call(cbind, X)[, colnames(corMat)])
  
  AvgFeatExp0 <- Xdat %>% mutate(Y = Y) %>% gather(Features, Exp, -Y) %>%
    group_by(Y, Features) %>% dplyr::summarise(Mean = mean(Exp), SD = sd(Exp))
  AvgFeatExp0$Dataset <- unlist(lapply(strsplit(AvgFeatExp0$Features, "_"), function(i) i[1]))
  featExp <- AvgFeatExp0 %>% group_by(Dataset, Y) %>% arrange(Mean)
  #featExp$Features <- unlist(lapply(strsplit(as.character(featExp$Features), "_"), function(i) i[2]))
  # Generate a circular plot (circos like) from a correlation matrix (pairwise)
  # 
  # Args:
  #   corMat: the main correlation matrix.
  #         -> colnames == rownames (pairwise); values = correlations
  #   featExp: data.frame holding the expression data.
  #   corThreshold: minimum value for correlations (<threshold will be ignored)
  #   figSize: figure size
  #   segmentWidth: thickness of the segment (main circle)
  #   linePlotWidth: thickness of the line plot (showing expression data)
  #   showIntraLinks = display links intra segments
  
  # 1) Generate karyotype data
  chr <- genChr(featExp);
  chr.names <- unique(chr$chrom);# paste("chr", 1:seg.num, sep="");
  # Calculate angles and band positions
  db <- segAnglePo(chr, seg=chr.names);
  db <- data.frame(db);
  
  # 2) Generate Links
  links <- genLinks(chr, corMat, corThreshold);
  if(nrow(links) < 1)
    stop("Choose a lower correlation threshold")
  
  # 3) Plot
  # Calculate parameters
  circleR <- (figSize / 2.0) -  segmentWidth - linePlotWidth;
  linksR <- circleR - segmentWidth;
  linePlotR <- circleR + segmentWidth
  chrLabelsR <- (figSize / 2.0)  ;
  
  par(mar=c(2, 2, 2, 2));
  
  plot(c(1,figSize), c(1,figSize), type="n", axes=FALSE, xlab="", ylab="", main="");
  
  # Plot ideogram
  drawIdeogram(R=circleR, cir=db, W=segmentWidth,  show.band.labels=TRUE, show.chr.labels=TRUE, chr.labels.R= chrLabelsR, chrData=chr, cex.label = cex.label);
  
  # Plot links
  drawLinks(R=linksR, cir=db,   mapping=links,   col=linkColors, drawIntraChr=showIntraLinks);
  
  # Plot expression values
  #cTypes <- unique(featExp[,1]) #Get the different disease/cancer types (lines)
  #lineCols <- rainbow(nrow(cTypes), alpha=0.5);          ## removed by Amrit
  cTypes <- levels(Y)                    ## added by amrit
  lineCols <- col.mixo[1:nlevels(Y)]         ## added by amrit
  
  # Fixme: remove this loop and send the whole expr dframe to drawLinePlot
  if(line==TRUE)
  {
      for (i in 1:length(chr.names)){
        seg.name <- gsub("chr","",chr.names[i]);
        #Get data for each segment
        expr <- subset(featExp,featExp$Dataset==seg.name)
        
        expr <- dcast(expr, formula = Features ~ Y, value.var="Mean")   ## changed PAM50 to Y
        expr <- merge(expr, chr, by.x="Features", by.y="name");
        expr$po <- (as.numeric(expr$chromStart) + as.numeric(expr$chromEnd)) / 2.0;
        expr <- dplyr::rename(expr, seg.name = chrom, seg.po = po);
        
        # Reorder columns
        cOrder <- c(c(grep("seg.name", colnames(expr)),
                      grep("seg.po", colnames(expr))), c(1:length(cTypes)+1))     ## changed 1:nrow to 1:length  removed c(1:length(cTypes)+1)
        expr <- expr[, cOrder];
        
        # Plot data on each sub segment
        subChr <- subset(db, db$seg.name == chr.names[i] )
        drawLinePlot(R=linePlotR, cir=subChr,   W=linePlotWidth, lineWidth=1, mapping=expr, col=lineCols, scale=FALSE);
      }
  }
  # Plot legend
  # First legeng bottom left corner
  legend(x=5, y = (circleR/4), title="Correlations", c("Positive Correlation", "Negative Correlation"), 
         col = c(colors()[134], colors()[128]), pch = 19, cex=0.8, bty = "n")
  # Second legend bottom righ corner
  if(line==TRUE)
  legend(x=figSize-(circleR/4), y = (circleR/3), title="Expression", legend=cTypes,  ## changed PAM50 to Y, changed by Amrit
         col = lineCols, pch = 19, cex=0.8, bty = "n")
  # third legend top left corner
  legend(x=figSize-(circleR/2), y = figSize, title="Correlation cut-off", legend=paste("r", corThreshold, sep = "="),
         col = "black", cex=0.8, bty = "n") 
}


plotDiablo2 = function (x, ncomp = 1, groupOrder, ...) 
{
  object = x
  opar <- par()[!names(par()) %in% c("cin", "cra", "csi", "cxy", 
                                     "din", "page")]
  indY = object$indY
  object$variates = c(object$variates[-indY], object$variates[indY])
  object$loadings = c(object$loadings[-indY], object$loadings[indY])
  VarX <- do.call(cbind, lapply(object$variates, function(i) i[, 
                                                               ncomp]))
  datNames <- colnames(VarX)
  if (ncol(VarX) <= 2) 
    stop("This function is only available when there are more than 3 blocks")
  Y = object$Y
  if (missing(groupOrder)) 
    groupOrder = levels(Y)
  if (!is.factor(Y)) 
    stop(gettextf("Y must be a factor!"))
  if (length(ncomp) != 1) 
    stop(gettextf("You can only choose one component"))
  numberOfCols <- ncol(VarX)
  numberOfRows <- numberOfCols - 1
  mat <- matrix(0, nrow = numberOfRows, ncol = numberOfRows)
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      mat[i, j] <- paste(i, j, sep = "_")
    }
  }
  plotType = list(cor = mat[lower.tri(mat)], scatter = mat[upper.tri(mat)], 
                  lab = diag(mat), bar = paste(1:(numberOfRows - 1), numberOfCols, 
                                               sep = "_"), stackedbar = paste(numberOfRows, numberOfCols, 
                                                                              sep = "_"))
  par(mfrow = c(numberOfRows, numberOfCols), mar = rep.int(1/2, 
                                                           4), oma = c(2, 2, 2, 2))
  for (i in 1:numberOfRows) {
    for (j in 1:numberOfCols) {
      ptype <- unlist(lapply(plotType, function(x) {
        intersect(paste(i, j, sep = "_"), x)
      }))
      splotMatPlot2(x = VarX[, i], y = VarX[, j], datNames, 
                   Y, ptype, groupOrder)
      if (i == 1 & j %in% seq(2, numberOfRows, 1)) {
        Axis(side = 3, x = VarX[, i])
      }
      if (j == numberOfRows & i %in% seq(1, numberOfRows - 
                                         1, 1)) {
        Axis(side = 4, x = VarX[, i])
      }
    }
  }
  par(opar)
}

splotMatPlot2 = function(x, y, datNames, Y, ptype, groupOrder){
  if(names(ptype) == "cor"){
    plot(1, type = "n", axes = FALSE)
    r = round(cor(x, y), 2)
    text(1, 1, labels=r, cex = abs(0.6/strwidth(r)*r))
    box()
  }
  if(names(ptype) == "scatter"){
    panel.ellipses(x=x, y=y, Y = Y)
  }
  if(names(ptype) == "lab"){
    plot(1, type = "n", axes = FALSE)
    ind = as.numeric(unlist(lapply(strsplit(ptype, "_"), unique)))
    text(x=1, y=1, labels=datNames[ind], cex = 2)
    box()
  }
  if(names(ptype) == "bar"){
    Y2 <- factor(as.character(Y), levels = groupOrder)
    par(las=2)
    boxplot(x ~ Y2, horizontal=TRUE, axes = FALSE, ylim = c(min(x)-3, max(x)),
            col= col.mixo[match(levels(Y2), levels(Y))])
    axis(4, at=1:nlevels(Y2), labels=levels(Y2))
  }
  if(names(ptype) == "stackedbar"){
    Y2 <- factor(as.character(Y), levels = groupOrder)
    bars <- table(Y2)
    par(las=1)
    barplot(bars, col= col.mixo[match(levels(Y2), levels(Y))], 
            axes = FALSE)
    axis(4, at=seq(0,max(bars),length.out=5), labels=seq(0,max(bars),length.out=5))
  }
}




################################################
#
# 4) circosPlot
#
################################################
circosPlot_diabloModif = function(object, corThreshold, cex.label,
                             figSize = 800,
                             segmentWidth = 25,
                             linePlotWidth = 90,
                             showIntraLinks = FALSE)
{
  options(stringsAsFactors = FALSE);
  set.seed(1234);
  
  ##############################
  ###   networkDiagram_core.R
  ###
  ###   Authors: Michael Vacher (minor changes by Amrit :)
  ###
  ###   Parts of this src has been modified from the original OmicCircos package obtained from: 
  ###   Ying Hu Chunhua Yan <yanch@mail.nih.gov> (2015). OmicCircos: High-quality circular visualization of omics data. R package version 1.6.0.
  ############################## 
  
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  
  par(mfrow  = c(1,1), mar = rep(1, 4))
  X <- object$X
  for(i in 1 : length(X)){
    colnames(X[[i]]) <- paste(names(X)[i], colnames(X[[i]]), sep = "_")
  }
  Y <- object$Y
  
  keepA = lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)
  cord = mapply(function(x, y, keep){
    cor(x[, keep], y, use = "pairwise")
  }, x=X, y=object$variates[-length(object$variates)], keep = keepA[-length(keepA)])
  
  simMatList <- vector("list", length(X))
  for(i in 1:length(cord)){
    for(j in 1:length(cord)){
      simMatList[[i]][[j]] <- cord[[i]] %*% t(cord[[j]])
    }
  }
  corMat <- do.call(rbind, lapply(simMatList, function(i) do.call(cbind, i)))
  
  ## Expression levels
  Xdat <- as.data.frame(do.call(cbind, X)[, colnames(corMat)])
  #corMat <- cor(Xdat)
  
  geneNames <-   unlist(lapply(strsplit(unlist(lapply(strsplit(grep("rna", colnames(Xdat), value = TRUE), "_"), function(x) x[2])), "\\."), function(i) i[1]))
  library(biomaRt)
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", 
                 dataset="hsapiens_gene_ensembl")
  filterList <- 'ensembl_gene_id'
  attr = c('ensembl_gene_id', 'hgnc_symbol')
  
  colnames(Xdat)[grep("rna", colnames(Xdat))] <- paste(
  "rna", getBM(attributes = attr,
                    filters = filterList,
                    values = geneNames,
                    mart = mart)$hgnc_symbol, sep = "_")
  
  AvgFeatExp0 <- Xdat %>% mutate(Y = Y) %>% gather(Features, Exp, -Y) %>%
    group_by(Y, Features) %>% dplyr::summarise(Mean = mean(Exp), SD = sd(Exp))
  AvgFeatExp0$Dataset <- rep(names(X), unlist(lapply(cord, nrow)))
  featExp <- AvgFeatExp0 %>% group_by(Dataset, Y) %>% arrange(Mean)
  rownames(corMat) <- colnames(corMat) <- colnames(Xdat)
  #featExp$Features <- unlist(lapply(strsplit(as.character(featExp$Features), "_"), function(i) paste(i[-1], collapse="_")))
  # Generate a circular plot (circos like) from a correlation matrix (pairwise)
  # 
  # Args:
  #   corMat: the main correlation matrix.
  #         -> colnames == rownames (pairwise); values = correlations
  #   featExp: data.frame holding the expression data.
  #   corThreshold: minimum value for correlations (<threshold will be ignored)
  #   figSize: figure size
  #   segmentWidth: thickness of the segment (main circle)
  #   linePlotWidth: thickness of the line plot (showing expression data)
  #   showIntraLinks = display links intra segments
  
  # 1) Generate karyotype data
  chr <- genChr(featExp);
  chr.names <- unique(chr$chrom);# paste("chr", 1:seg.num, sep="");
  # Calculate angles and band positions
  db <- segAnglePo(chr, seg=chr.names);
  db <- data.frame(db);
  
  # 2) Generate Links
  links <- genLinks(chr, corMat, corThreshold);
  if(nrow(links) < 1)
    stop("Choose a lower correlation threshold")
  
  # 3) Plot
  # Calculate parameters
  circleR <- (figSize / 2.0) -  segmentWidth - linePlotWidth;
  linksR <- circleR - segmentWidth;
  linePlotR <- circleR + segmentWidth
  chrLabelsR <- (figSize / 2.0)  ;
  
  par(mar=c(2, 2, 2, 2));
  
  plot(c(1,figSize), c(1,figSize), type="n", axes=FALSE, xlab="", ylab="", main="");
  
  # Plot ideogram
  drawIdeogram(cex.label=cex.label, R=circleR, cir=db, W=segmentWidth,  show.band.labels=TRUE, show.chr.labels=TRUE, chr.labels.R= chrLabelsR, chrData=chr);
  
  # Plot links
  drawLinks(R=linksR, cir=db,   mapping=links,   col=linkColors, drawIntraChr=showIntraLinks);
  
  # Plot expression values
  #cTypes <- unique(featExp[,1]) #Get the different disease/cancer types (lines)
  #lineCols <- rainbow(nrow(cTypes), alpha=0.5);          ## removed by Amrit
  cTypes <- levels(Y)                    ## added by amrit
  lineCols <- col.mixo[1:nlevels(Y)]         ## added by amrit
  
  # Fixme: remove this loop and send the whole expr dframe to drawLinePlot
  for (i in 1:length(chr.names)){
    seg.name <- gsub("chr","",chr.names[i]);
    #Get data for each segment
    expr <- subset(featExp,featExp$Dataset==seg.name)
    
    expr <- dcast(expr, formula = Features ~ Y, value.var="Mean")   ## changed PAM50 to Y
    expr <- merge(expr, chr, by.x="Features", by.y="name");
    expr$po <- (as.numeric(expr$chromStart) + as.numeric(expr$chromEnd)) / 2.0;
    expr <- dplyr::rename(expr, seg.name = chrom, seg.po = po);
    
    # Reorder columns
    cOrder <- c(c(grep("seg.name", colnames(expr)),
                  grep("seg.po", colnames(expr))), c(1:length(cTypes)+1))     ## changed 1:nrow to 1:length  removed c(1:length(cTypes)+1)
    expr <- expr[, cOrder];
    
    # Plot data on each sub segment
    subChr <- subset(db, db$seg.name == chr.names[i] )
    drawLinePlot(R=linePlotR, cir=subChr,   W=linePlotWidth, lineWidth=1, mapping=expr, col=lineCols, scale=FALSE);
  }
  
  # Plot legend
  # First legeng bottom left corner
  legend(x=5, y = (circleR/4), title="Correlations", c("Positive Correlation", "Negative Correlation"), 
         col = c(colors()[134], colors()[128]), pch = 19, cex=0.8, bty = "n")
  # Second legend bottom righ corner
  legend(x=figSize-(circleR/4), y = (circleR/3), title="Expression", legend=cTypes,  ## changed PAM50 to Y, changed by Amrit
         col = lineCols, pch = 19, cex=0.8, bty = "n")
  # third legend top left corner
  legend(x=figSize-(circleR/2), y = figSize, title="Correlation cut-off", legend=paste("r", corThreshold, sep = "="),
         col = "black", cex=0.8, bty = "n") 
}

col.mixo = color.mixo(1:10)
plotDiablo3=function (x, ncomp = 1, groupOrder, ...) 
{
  object = x
  opar <- par()[!names(par()) %in% c("cin", "cra", "csi", "cxy", 
                                     "din", "page")]
  indY = object$indY
  object$variates = c(object$variates[-indY], object$variates[indY])
  object$loadings = c(object$loadings[-indY], object$loadings[indY])
  VarX <- do.call(cbind, lapply(object$variates, function(i) i[, 
                                                               ncomp]))
  datNames <- colnames(VarX)
  if (ncol(VarX) <= 2) 
    stop("This function is only available when there are more than 3 blocks")
  Y = object$Y
  if (missing(groupOrder)) 
    groupOrder = levels(Y)
  if (!is.factor(Y)) 
    stop(gettextf("Y must be a factor!"))
  if (length(ncomp) != 1) 
    stop(gettextf("You can only choose one component"))
  numberOfCols <- ncol(VarX)
  numberOfRows <- numberOfCols - 1
  mat <- matrix(0, nrow = numberOfRows, ncol = numberOfRows)
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      mat[i, j] <- paste(i, j, sep = "_")
    }
  }
  plotType = list(cor = mat[lower.tri(mat)], scatter = mat[upper.tri(mat)], 
                  lab = diag(mat), bar = paste(1:(numberOfRows - 1), numberOfCols, 
                                               sep = "_"), stackedbar = paste(numberOfRows, numberOfCols, 
                                                                              sep = "_"))
  par(mfrow = c(numberOfRows, numberOfCols), mar = rep.int(1/2, 
                                                           4), oma = c(2, 2, 2, 2))
  for (i in 1:numberOfRows) {
    for (j in 1:numberOfCols) {
      ptype <- unlist(lapply(plotType, function(x) {
        intersect(paste(i, j, sep = "_"), x)
      }))
      splotMatPlot3(x = VarX[, i], y = VarX[, j], datNames, 
                   Y, ptype, groupOrder)
      if (i == 1 & j %in% seq(2, numberOfRows, 1)) {
        Axis(side = 3, x = VarX[, i])
      }
      if (j == numberOfRows & i %in% seq(1, numberOfRows - 
                                         1, 1)) {
        Axis(side = 4, x = VarX[, i])
      }
    }
  }
  par(opar)
}
splotMatPlot3 = function(x, y, datNames, Y, ptype, groupOrder){
  if(names(ptype) == "cor"){
    plot(1, type = "n", axes = FALSE)
    r = round(cor(x, y), 2)
    text(1, 1, labels=r, cex = abs(0.6/strwidth(r)*r))
    box()
  }
  if(names(ptype) == "scatter"){
    panel.ellipses(x=x, y=y, Y = Y)
  }
  if(names(ptype) == "lab"){
    plot(1, type = "n", axes = FALSE)
    ind = as.numeric(unlist(lapply(strsplit(ptype, "_"), unique)))
    text(x=1, y=1, labels=datNames[ind], cex = 2)
    box()
  }
  if(names(ptype) == "bar"){
    Y2 <- factor(as.character(Y), levels = groupOrder)
    par(las=2)
    boxplot(x ~ Y2, horizontal=TRUE, axes = FALSE, ylim = c(min(x)-3, max(x)),
            col= col.mixo[match(levels(Y2), levels(Y))])
    axis(4, at=1:nlevels(Y2), labels=levels(Y2))
  }
  if(names(ptype) == "stackedbar"){
    Y2 <- factor(as.character(Y), levels = groupOrder)
    bars <- table(Y2)
    par(las=1)
    barplot(bars, col= col.mixo[match(levels(Y2), levels(Y))], 
            axes = FALSE)
    axis(4, at=seq(0,max(bars),length.out=5), labels=seq(0,max(bars),length.out=5))
  }
}