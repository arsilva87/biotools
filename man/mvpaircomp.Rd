\name{mvpaircomp}
\alias{mvpaircomp}
\alias{print.mvpaircomp}

\title{
Multivariate Pairwise Comparisons
}

\description{
Performs pairwise comparisons of multivariate mean vectors of factor levels, overall or nested.
The tests are run in the same spirt of \code{summary.manova()}, based on multivariate statistics such as Pillai's trace
and Wilks' lambda, which can be applied to test multivariate contrasts.
}

\usage{
mvpaircomp(model, factor1, nesting.factor = NULL,
   test = "Pillai", adjust = "none", SSPerror = NULL, DFerror = NULL)

\method{print}{mvpaircomp}(x, ...)
}

\arguments{
  \item{model}{ a multivariate analysis of variance (MANOVA) model, fitted using \code{lm()} or \code{manova()}; an object of class \code{"mlm"}.}
  \item{factor1}{ a character string indicating a factor declared in the \code{model}, whose levels will be compared in pairs.}
  \item{nesting.factor}{ optional; a character string indicating a factor also declared in \code{model} whose levels will nest the contrasts performed
	with \code{factor1}. \code{factor1} and \code{nesting.factor} can have pretty much any form or relationship, that is,
	not necessarily nested one another.}
  \item{test}{ a character string indicating the type of multivariate statistics to be calculated to perform the
	F-test approximation. Default is \code{"Pillai"}. Other options are: \code{"Wilks"}, \code{"Hotelling-Lawley"} and \code{"Roy"}.
	But they use to give very close results.}
  \item{adjust}{ a character string indicating the p-value adjustment method for multiple comparisons. Default is \code{"none"}.
	See \code{\link[stats]{p.adjust}}}
  \item{SSPerror}{ optional; a numeric matrix representing the residual sum of squares and cross-products, to be used to compute the
	multivariate statistics.}
  \item{DFerror}{ optional; a numeric value representing the residual degrees of freedom, to be used to compute the
	multivariate statistics.}
  \item{x}{ an object of class \code{mvpaircomp}.}
  \item{\dots}{ further arguments.}
}

\value{
An object of class \code{mvpaircomp}, a list of
  \item{st}{ an array containing the summary of the multivariate tests.}
  \item{SSPcontrast}{ an array containing p-dimensional square matrices of sum of squares and cross-products of the contrasts.}
  \item{adjust}{ a character string indicating the p-value adjustment method used.}
  \item{fac1}{ a character string indicating the factor being tested.}
  \item{fac2}{ a character string indicating the nesting factor.}
}

\references{
Krzanowski, W. J. (1988) \emph{Principles of Multivariate Analysis}. A User's Perspective. Oxford.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link[stats]{manova}}
}

\examples{
# Example 1
data(maize)
M <- lm(cbind(NKPR, ED, CD, PH) ~ family + env, data = maize)
anova(M)  # MANOVA table
mvpaircomp(M, factor1 = "family", adjust = "bonferroni")

# Example 2 (with nesting factor)
# Data on producing plastic film from Krzanowski (1998, p. 381)
tear <- c(6.5, 6.2, 5.8, 6.5, 6.5, 6.9, 7.2, 6.9, 6.1, 6.3,
          6.7, 6.6, 7.2, 7.1, 6.8, 7.1, 7.0, 7.2, 7.5, 7.6)
gloss <- c(9.5, 9.9, 9.6, 9.6, 9.2, 9.1, 10.0, 9.9, 9.5, 9.4,
           9.1, 9.3, 8.3, 8.4, 8.5, 9.2, 8.8, 9.7, 10.1, 9.2)
opacity <- c(4.4, 6.4, 3.0, 4.1, 0.8, 5.7, 2.0, 3.9, 1.9, 5.7,
             2.8, 4.1, 3.8, 1.6, 3.4, 8.4, 5.2, 6.9, 2.7, 1.9)
Y <- cbind(tear, gloss, opacity)
rate     <- gl(2, 10, labels = c("Low", "High"))
additive <- gl(2, 5, length = 20, labels = c("Low", "High"))

fit <- manova(Y ~ rate * additive)
summary(fit, test = "Wilks")  # MANOVA table
mvpaircomp(fit, factor1 = "rate", nesting.factor = "additive", test = "Wilks")
mvpaircomp(fit, factor1 = "additive", nesting.factor = "rate", test = "Wilks")

# End (not run)
}
