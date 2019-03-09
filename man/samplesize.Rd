\name{samplesize}
\alias{samplesize}

\title{ Minimum Sample Size }

\description{
Function to determine the minimum sample size for calculating a
statistic based on its the confidence interval.
}

\usage{
samplesize(x, fun, sizes = NULL, lcl = NULL, ucl = NULL,
	nboot = 200, conf.level = 0.95, nrep = 500, graph = TRUE, ...)
}

\arguments{
  \item{x}{ a numeric vector. }
  \item{fun}{ an objective function at which to evaluate the sample size; see details.}
  \item{sizes}{ a numeric vector containing sample sizes; if \code{NULL} (default), \code{samplesize} 
	creates a vector ranging from 2 to \emph{n}-1. }
  \item{lcl}{ the lower confidence limit for the statistic defined in \code{fun}; 
	if \code{NULL} (default), \code{samplesize} estimates \code{lcl} based on bootstrap 
	percentile interval. }
  \item{ucl}{ the upper confidence limit for the statistic defined in \code{fun}; 
	if \code{NULL} (default), \code{samplesize} estimates \code{ucl} based on bootstrap 
	percentile interval. }
  \item{nboot}{ the number of bootstrap samples; it is used only if \code{lcl} or \code{ucl} is \code{NULL}. }
  \item{conf.level}{ the confidence level for calculating the \code{lcl} and \code{ucl}; 
	it is used only if \code{lcl} or \code{ucl} is \code{NULL}. }
  \item{nrep}{ the resampling (with replacement) number for each sample size in \code{sizes}; 
	default is 500. }
  \item{graph}{ logical; default is \code{TRUE}. }
  \item{\dots}{ further graphical arguments. }
}

\details{
If \code{ucl} or \code{lcl} is \code{NULL}, \code{fun} must be defined as in \code{\link{boot}}, i.e.,
the first argument passed will always be the original data and the second will be a vector of indices,
frequencies or weights which define the bootstrap sample. By now, \code{samplesize} considers the
second argument only as index.
}

\value{
A list of
  \item{CI }{a vector containing the lower and the upper confidence limit for the statistic evaluated.}
  \item{pointsOut }{a data frame containing the sample sizes (in \code{sizes}), the number of 
	points outside the CI (\code{n.out}) and the proportion of this number (\code{prop}). }
}

\section{Side Effects }{
If \code{graph = TRUE}, a graphic with the dispersion of the estimates for each sample size, 
as well as the graphic containing the number of points outside the confidence interval for 
the reference sample.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\examples{
cv <- function(x, i) sd(x[i]) / mean(x[i]) # coefficient of variation
x = rnorm(20, 15, 2)
cv(x)
samplesize(x, cv)

par(mfrow = c(1, 3), cex = 0.7, las = 1)
samplesize(x, cv, lcl = 0.05, ucl = 0.20)
abline(h = 0.05 * 500, col = "blue") # sample sizes with 5\% (or less) out CI

# End (not run)
}