\name{mantelPower}
\alias{mantelPower}

\title{Power of Mantel's Test}

\description{Power calculation of Mantel's permutation test.}

\usage{
mantelPower(obj, effect.size = seq(0, 1, length.out = 50), alpha = 0.05)
}

\arguments{
  \item{obj}{ an object of class "mantelTest". See \code{\link[biotools]{mantelTest}}. }
  \item{effect.size}{ numeric; the effect size specifying the alternative hypothesis.}
  \item{alpha}{ numeric; the significance level at which to compute the power level.}
}

\value{
A data frame containing the effect size and its respective power level.
}

\references{
Silva, A.R.; Dias, C.T.S.; Cecon, P.R.; Rego, E.R. (2015). An alternative procedure for performing a
	power analysis of Mantel's test. \emph{Journal of Applied Statistics}, \bold{42}(9): 1984-1992.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link[biotools]{mantelTest}}
}

\examples{
# Mantel test
data(garlicdist)
garlic <- tocher(garlicdist)
coph <- cophenetic(garlic)
mt1 <- mantelTest(garlicdist, coph, xlim = c(-1, 1))

# Power calculation, H1: rho = 0.3
mantelPower(mt1, effect.size = 0.3)

# Power calculation, multiple H1s and different alphas
p01 <- mantelPower(mt1, alpha = 0.01)
p05 <- mantelPower(mt1, alpha = 0.05)
p10 <- mantelPower(mt1, alpha = 0.10)
plot(p01, type = "l", col = 4)
lines(p05, lty = 2, col = 4)
lines(p10, lty = 3, col = 4)
legend("bottomright", c("0.10", "0.05", "0.01"),
	title = expression(alpha), col = 4, lty = 3:1, cex = 0.8)

# End (Not run)
}
