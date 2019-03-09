\name{singh}
\alias{singh}
\alias{singh.default}
\alias{plot.singh}

\title{ Importance of Variables According to the Singh (1981) Criterion }

\description{
A function to calculate the Singh (1981) criterion for importance
of variables based on the squared generalized Mahalanobis distance.
\deqn{S_{.j} = \sum_{i=1}^{n-1} \sum_{i'>i}^{n} (x_{ij} - x_{i'j}) * (\bold{x}_i - \bold{x}_{i'})' * \bold{\Sigma}_{j}^{-1} }
}

\usage{
\method{singh}{default}(data, cov, inverted = FALSE)
\method{plot}{singh}(x, ...)
}

\arguments{
  \item{data}{ a data frame or matrix of data (\emph{n x p}). }
  \item{cov}{ a variance-covariance matrix (\emph{p x p}). }
  \item{inverted}{ logical. If \code{FALSE} (default), \code{cov} is supposed to be a variance-covariance matrix. }
  \item{x}{ an object of class \code{"singh"}. }
  \item{\dots}{ further graphical arguments. }
}

\value{
\code{singh} returns a matrix containing the Singh statistic, the
importance proportion and the cummulative proprtion of each
variable (column) in data.
}

\references{
Singh, D. (1981) The relative importance of characters affecting genetic
   divergence. \emph{Indian Journal Genetics & Plant Breeding}, 41:237-245.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{D2.dist}}
}

\examples{
# Manly (2004, p.65-66)
x1 <- c(131.37, 132.37, 134.47, 135.50, 136.17)
x2 <- c(133.60, 132.70, 133.80, 132.30, 130.33)
x3 <- c(99.17, 99.07, 96.03, 94.53, 93.50)
x4 <- c(50.53, 50.23, 50.57, 51.97, 51.37)
x <- cbind(x1, x2, x3, x4)
Cov <- matrix(c(21.112,0.038,0.078,2.01, 0.038,23.486,5.2,2.844,
	0.078,5.2,24.18,1.134, 2.01,2.844,1.134,10.154), 4, 4)
(s <- singh(x, Cov))
plot(s)

# End (not run)
}
