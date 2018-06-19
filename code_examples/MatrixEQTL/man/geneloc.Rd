\name{geneloc}
\alias{geneloc}
\docType{data}
\title{
	Artificial data for Matrix eQTL sample code: Gene location file.
}
\description{
	Artificial Gene location file for 10 genes. 
	\tabular{llll}{
		geneid  \tab chr   		\tab left 	\tab right \cr
		Gene_01 \tab chr1     \tab 721289 \tab 731289 \cr
		Gene_02 \tab chr1     \tab 752565 \tab 762565 \cr
		\ldots  \tab \ldots   \tab \ldots \tab \ldots \cr
	}
}
\format{
  A data frame with 4 columns.
  \describe{
    \item{\code{geneid}}{A column with gene names. The order does not have to match the gene expression data set.}
    \item{\code{chr}}{Chromosome number, i.e. chr1.}
    \item{\code{left}}{Lower (smaller) coordinate of the gene begining/end.}
    \item{\code{right}}{Upper (larger) coordinate of the gene begining/end.}
  }
}
\references{
	The package website: \url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\seealso{
	See \code{\link{Matrix_eQTL_engine}} for reference and sample code.
}
\keyword{datasets}
