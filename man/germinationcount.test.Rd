\name{germinationcount.test}
\alias{germinationcount.test}

\title{ Germination Count Range Test for Seed Sample Heterogeneity }

\description{
A test based on the exact probability distribution of the germination count range, i.e, the difference between germination count of seed samples.
}

\usage{
germinationcount.test(r, nsamples, n, N, K)
}

\arguments{
  \item{r}{ an integer representing the germination count difference between seed samples.} 
  \item{nsamples}{ an integer representing the number of seed samples.}
  \item{n}{ an integer representing the number of seeds per sample.}
  \item{N}{ an integer representing the size (number of seeds) of the seed lot.}
  \item{K}{ an integer representing the number of germinating seeds in the seed lot.}
}

\value{
A list of
  \item{R.value}{ integer; the input R-value (\code{r}).}
  \item{p.value}{ numeric; the exact p-value.}
  \item{germination.rate}{ numeric; the germination rate of the seed lot, calculated as the ration of \code{K} and \code{N}.}
}

\references{
Da Silva, A.R. (2020). On testing for seed sample heterogeneity with the exact probability distribution of the germination count range. 
	\emph{Seed Science Research}, \bold{30}(1): 59--63. doi:10.1017/S0960258520000112
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\examples{
germinationcount.test(r = 6, nsamples = 4, n = 50, N = 2000, K = 1700)

# End (Not run)
}
