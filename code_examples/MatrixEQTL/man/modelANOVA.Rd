\name{modelANOVA}
\alias{modelANOVA}
\docType{data}
\title{
	Constant for \code{\link{Matrix_eQTL_engine}}.
}
\description{
	Set parameter \code{useModel = modelANOVA} in the call of \code{\link{Matrix_eQTL_main}} to indicate that the genotype should be treated as a categorical variable.
}
\note{
By default, the number of ANOVA categories is fixed to be 3.

To set it to a different number (say, 4) use the following command:

\code{options(MatrixEQTL.ANOVA.categories=4)}

To check the current settings run:

\code{getOption('MatrixEQTL.ANOVA.categories', 3);}
}
\examples{
library('MatrixEQTL')	
			
# Number of columns (samples)
n = 100;

# Number of covariates
nc = 10;

# Generate the standard deviation of the noise
noise.std = 0.1 + rnorm(n)^2;

# Generate the covariates
cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);

# Generate the vectors with single genotype and expression variables
snps.mat = floor(runif(n, min = 0, max = 3));
gene.mat = 1 + (snps.mat==1) + cvrt.mat \%*\% rnorm(nc) + rnorm(n) * noise.std;

# Create 3 SlicedData objects for the analysis
snps1 = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
gene1 = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrt1 = SlicedData$new( t(cvrt.mat) );

# name of temporary output file
filename = tempfile();

snps1
gene1

# Call the main analysis function
me = Matrix_eQTL_main(
  snps = snps1, 
  gene = gene1, 
  cvrt = cvrt1, 
  output_file_name = filename, 
  pvOutputThreshold = 1, 
  useModel = modelANOVA, 
  errorCovariance = diag(noise.std^2), 
  verbose = TRUE,
  pvalue.hist = FALSE );
# remove the output file
unlink( filename );

# Pull Matrix eQTL results - t-statistic and p-value

fstat = me$all$eqtls$statistic;
pvalue = me$all$eqtls$pvalue;
rez = c( Fstat = fstat, pvalue = pvalue)
# And compare to those from ANOVA in R
{
  cat('\n\n Matrix eQTL: \n'); 
  print(rez);
  cat('\n R anova(lm()) output: \n')
  lmodel = lm( gene.mat ~ cvrt.mat + factor(snps.mat), weights = 1/noise.std^2 );
  lmout = anova( lmodel )[2, 4:5];
  print( lmout )
}

# Results from Matrix eQTL and 'lm' must agree
stopifnot(all.equal(lmout, rez, check.attributes=FALSE))
}
\references{
	The package website: \url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\seealso{
	See \code{\link{Matrix_eQTL_engine}} for reference and sample code.
}
