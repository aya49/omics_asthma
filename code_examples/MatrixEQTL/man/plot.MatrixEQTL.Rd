\name{plot.MatrixEQTL}
\alias{plot.MatrixEQTL}
\title{Plot histogram or QQ-plot of all p-values
}
\description{
This method plots a histogram or QQ-plot of p-values
for all tests performed by \code{\link{Matrix_eQTL_engine}}.
}
\usage{
\method{plot}{MatrixEQTL}(x, cex = 0.5, pch = 19, xlim = NULL, ylim = NULL,...)
}
\arguments{
	\item{x}{
		An object returned by \code{\link{Matrix_eQTL_engine}}.
	}
	\item{cex}{
		A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
	}
	\item{pch}{
		Plotting "character", i.e., symbol to use. See \code{\link[graphics]{points}}.
	}
	\item{xlim}{
		Set the range of the horisontal axis.
	}
	\item{ylim}{
		Set the range of the vertical axis.
	}
	\item{\dots}{
		further graphical parameters passed to \code{\link[graphics]{lines}} and \code{\link[graphics]{points}}.
	}
}
\details{
	The plot type (histogram vs. QQ-plot) is determined by the \code{pvalue.hist} parameter in the call of \code{\link{Matrix_eQTL_engine}} function.
}
\value{
	The method does not return any value.
}
\references{
	The package website: \url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\seealso{
	See \code{\link{Matrix_eQTL_engine}} for reference and sample code.
}
\note{
	The sample code below produces figures like these:
	
	Histogram: \figure{histogram.png}
		
	QQ-plot: \figure{QQplot.png}
}
\author{
	Andrey Shabalin \email{ashabalin@vcu.edu}
}
\examples{
library(MatrixEQTL)

# Number of samples
n = 100; 

# Number of variables
ngs = 2000;

# Common signal in all variables
pop = 0.2*rnorm(n);

# data matrices
snps.mat = matrix(rnorm(n*ngs), ncol = ngs) + pop;
gene.mat = matrix(rnorm(n*ngs), ncol = ngs) + pop + snps.mat*((1:ngs)/ngs)^9/2;

# data objects for Matrix eQTL engine
snps1 = SlicedData$new( t( snps.mat ) );
gene1 = SlicedData$new( t( gene.mat ) );
cvrt1 = SlicedData$new( );
rm(snps.mat, gene.mat)

# Slice data in blocks of 500 variables
snps1$ResliceCombined(500);
gene1$ResliceCombined(500);

# Produce no output files
filename = NULL; # tempfile()

# Perform analysis recording information for a histogram
meh = Matrix_eQTL_engine(
    snps = snps1, 
    gene = gene1, 
    cvrt = cvrt1, 
    output_file_name = filename, 
    pvOutputThreshold = 1e-100, 
    useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    pvalue.hist = 100);
    
plot(meh, col="grey")

# Perform analysis recording information for a QQ-plot
meq = Matrix_eQTL_engine(
    snps = snps1, 
    gene = gene1, 
    cvrt = cvrt1, 
    output_file_name = filename,
    pvOutputThreshold = 1e-6, 
    useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    pvalue.hist = "qqplot");
    
plot(meq)
}
\keyword{ QQ-plot }
\keyword{ histogram }
