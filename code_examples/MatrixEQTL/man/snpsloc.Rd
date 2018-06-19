\name{snpsloc}
\alias{snpsloc}
\docType{data}
\title{
	Artificial data for Matrix eQTL sample code: SNP location file.
}
\description{
	Artificial SNP location file for 15 markers. 
	\tabular{llll}{
		snpid   \tab chr   		\tab pos 	   \cr
		Snp_01  \tab chr1     \tab 721289  \cr
		Snp_02  \tab chr1     \tab 752565  \cr
		\ldots  \tab \ldots   \tab \ldots  \cr
	}
}
\format{
  A data frame with 3 columns.
  \describe{
    \item{\code{snpid}}{A column with SNP names. The order does not have to match the genotype data set.}
    \item{\code{chr}}{Chromosome number, i.e. chr1.}
    \item{\code{pos}}{Coordinate of the SNP.}
  }
}
\references{
	The package website: \url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\seealso{
	See \code{\link{Matrix_eQTL_engine}} for reference and sample code.
}
\keyword{datasets}
