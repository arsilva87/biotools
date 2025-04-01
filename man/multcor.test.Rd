\name{multcor.test}
\alias{multcor.test}

\title{ Pairwise Correlation t-Test }

\description{
It performs multiple correlation t-tests from a correlation matrix based
on the statistic:
	\deqn{t = r * \sqrt(df / (1 - r^2))}
where, in general, \eqn{df = n - 2}.
}

\usage{
multcor.test(x, n = NULL, Df = NULL,
	alternative = c("two.sided", "less", "greater"), adjust = "none")
}

\arguments{
  \item{x}{ a correlation matrix.}
  \item{n}{ the number of observations; if \code{NULL} (default), the argument
	\code{Df} must be passed.}
  \item{Df}{ the number of degrees of freedom of the t statistic; if \code{NULL} (default),
	the argument \code{n} must be passed and, in this case, \code{multcor.test}
	considers \eqn{Df = n - 2}.}
  \item{alternative}{ the alternative hypothesis. It must be one of "two.sided",
	"greater" or "less". You can specify just the initial letter.
	"greater" corresponds to positive association,"less" to negative association.
	The default is "two.sided".}
  \item{adjust}{ The adjustment method for multiple tests. It must be one of
	"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
	"none" (default). For more information, see \code{\link[stats]{p.adjust}}.}
}


\value{
A list with class "multcor.test" containing the following components:
  \item{t.values }{the t-value calculated for each correlation.}
  \item{p.values }{the p.value for each t-test, adjusted for multiple tests.}
  \item{p.check }{a matrix containing the \code{p.values} for each t-test (lower triangular)
	and a symbol indicating the significance level at which one can to reject the null hypothesis
	(upper triangular).}
  \item{adjustemnt }{a character indicating the p-value adjustment method.}
  \item{df }{the degrees of freedom of the tests.}
  \item{alternative }{a character indicating the type of alternative hypothesis.}
  \item{data.name }{a character string giving the name of the data.}
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link[stats]{cor}}, \code{\link[stats]{cor.test}}, \code{\link[stats]{p.adjust}}
}

\examples{
data(peppercorr)
multcor.test(peppercorr, n = 20)

# End (not run)
}
