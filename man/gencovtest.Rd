\name{gencovtest}
\alias{gencovtest}
\alias{gencovtest.manova}
\alias{print.gencovtest}
\alias{plot.gencovtest}

\title{
Testing Genetic Covariance
}

\description{
\code{gencovtest()} tests genetic covariance components from a MANOVA model. Two different approaches can
be used: (I) a test statistic that takes into account the genetic and environmental effects and (II) a test
statistic that only considers the genetic information. The first type refers to tests based on the mean
cross-products ratio, whose distribution is obtained via Monte Carlo simulation of Wishart matrices. The
second way of testing genetic covariance refers to tests based upon an adaptation of Wilks' and Pillai's
statistics for evaluating independence of two sets of variables. All these tests are described by Silva (2015).
}

\usage{
\method{gencovtest}{manova}(obj, geneticFactor, gcov = NULL,
	residualFactor = NULL, adjNrep = 1,
	test = c("MCPR", "Wilks", "Pillai"),
	nsim = 9999,
	alternative = c("two.sided", "less", "greater"))
\method{print}{gencovtest}(x, digits = 4, ...)
\method{plot}{gencovtest}(x, var1, var2, ...)
}

\arguments{
  \item{obj}{ an object of class \code{"manova"}.}
  \item{geneticFactor}{ a character indicating the genetic factor from which to test covariance components.
	It must be declared as a factor in the manova object.}
  \item{gcov}{ optional; a matrix containing estimates of genetic covariances to be tested. If
	\code{NULL} (default), an estimate is obtained via method of moments.}
  \item{residualFactor}{ optional; a character indicating a source in the manova model to be used as
	error term. If \code{NULL} (default), the usual term "Residuals" will be used.}
  \item{adjNrep}{ a correction index for dealing with unbalanced data. See details.}
  \item{test}{ a character indicating the test. It must be on of the following:
	\code{"MCPR"} - the empirical type-I test based on Mean Cross-Products Ratios via Wishart simulation,
	\code{"Wilks"} - a type-II test based on the partial Wilks' Lambda,
	\code{"Pillai"} - a type-II test based on the partial Pillai's statistic.}
  \item{nsim}{ the number of Monte Carlo simulations. Used only if \code{test = "MCPR"}.}
  \item{alternative}{ the type of alternative hypothesis. Used only if \code{test = "MCPR"}. So far,
	only the option \code{"two.sided"} is implemented.}
  \item{x}{ an object of class \code{"gencovtest"}.}
  \item{digits}{ the number of digits to be displayed by the print method.}
  \item{var1}{ a character of integer indicating one of the two response variable or its position.}
  \item{var2}{ a character of integer indicating one of the two response variable or its position.}
  \item{\dots}{ further arguments.}
}

\details{
The genetic covariance matrix is currently estimated via method of moments, following the equation:
\deqn{G = (Mg - Me) / (nrep * adjNrep)}
where \eqn{Mg} and \eqn{Me} are the matrices of mean cross-products associated with the genetic factor and
the residuals, respectively; \eqn{nrep} is the number of replications, calculated as the ratio between the
total number of observations and the number of levels of the genetic factor; \eqn{adjNrep} is supposed to
adjust nrep, specially when estimating \eqn{G} from unbalanced data.
}

\value{
An object of class \code{gencovtest}, a list of
  \item{gcov}{ a p-dimensional square matrix containing estimates of the genetic covariances.}
  \item{gcor}{ a p-dimensional square matrix containing estimates of the genetic correlations.}
  \item{test}{ the test (as input).}
  \item{statistics}{ a p-dimensional square matrix containing the \code{test} statistics.
	If \code{test = "MCPR"} the mean cross-products ratios are computed; if \code{test = "Wilks"}
	the Wilks' Lambda is; and \code{test = "Pillai"} results on Pillai's \eqn{Tn}.}
  \item{p.values}{ a p-dimensional square matrix containing the associated p-values.}
  \item{alternative}{ the type of alternative hypothesis (as input).}
  \item{X2}{ a p-dimensional square matrix containing the Chi-square (D.f. = 1) approximation for Wilks's
	and Pillai's statistics. Stored only if one of these two tests is chosen.}
  \item{simRatio}{ an array consisting of \code{nsim} p-dimensional matrices containing the simulated mean
	cross-products ratios.}
  \item{dfg}{ the number of degrees of freedom associated with the genetic factor.}
  \item{dfe}{ the number of degrees of freedom associated with the residual term.}
}

\references{
Silva, A.R. (2015) \emph{On Testing Genetic Covariance}. LAP Lambert Academic Publishing. ISBN 3659716553
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\section{Warning}{
When using the MCPR test, be aware that \code{dfg} should be equal or greater than the number of variables (p).
Otherwise the simulation of Wishart matrices may not be done.

A collinearity diagnosis is carried out using the condition number (CN), for the inferences may be affected by the
quality of \eqn{G}. Thus, if CN > 100, a warning message is displayed.
}

\seealso{
\code{\link[stats]{manova}}
}

\examples{
# MANOVA
data(maize)
M <- manova(cbind(NKPR, ED, CD, PH) ~ family + env, data = maize)
summary(M)

# Example 1 - MCPR
t1 <- gencovtest(obj = M, geneticFactor = "family")
print(t1)
plot(t1, "ED", "PH")

# Example 2 - Pillai
t2 <- gencovtest(obj = M, geneticFactor = "family", test = "Pillai")
print(t2)
plot(t2, "ED", "PH")

# End (not run)
}
