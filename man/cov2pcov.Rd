\name{cov2pcov}
\alias{cov2pcov}

\title{ Partial Covariance Matrix }

\description{
Compute a matrix of partial (co)variances for a group of variables with respect to another.

Take \eqn{\Sigma} as the covariance matrix of dimension \emph{p}. Now consider dividing \eqn{\Sigma} into two groups
of variables. The partial covariance matrices are calculate by:
\deqn{ \Sigma_{11.2} = \Sigma_{11} - \Sigma_{12} \Sigma_{22}^{-1} \Sigma_{21} }
\deqn{ \Sigma_{22.1} = \Sigma_{22} - \Sigma_{21} \Sigma_{11}^{-1} \Sigma_{12} }
}

\usage{
cov2pcov(m, vars1, vars2 = seq(1, ncol(m))[-vars1])
}

\arguments{
  \item{m}{ a square numeric matrix.}
  \item{vars1}{ a numeric vector indicating the position (rows or columns in \code{m}) of the set of variables
	at which to compute the partial covariance matrix.}
  \item{vars2}{ a numeric vector indicating the position (rows or columns in \code{m}) of the set of variables
	at which to adjust the partial covariance matrix.}
}

\value{ A square numeric matrix.}

\author{
Anderson Rodrigo da Silva <anderson.agro at hotmail.com>
}

\seealso{
\code{\link[stats]{cov}}
}

\examples{
(Cl <- cov(longley))
cov2pcov(Cl, 1:2)

# End (Not run)
}
