\name{fitplotsize}
\alias{fitplotsize}

\title{ Parameter Estimation of the Plot Size Model }

\description{
Function to estimate the parameters of the nonlinear Lessman & Atkins (1963) model
for determining the optimum plot size as a function of the experimental
coefficient of variation (CV) or as a function of the residual standard error.
\deqn{CV = a * plotsize ^ {-b}.}
It creates initial estimates of the parameters \emph{a} and \emph{b} by log-linearization
and uses them to provide its least-squares estimates via \code{\link[stats]{nls}}.
}

\usage{
fitplotsize(plotsize, CV)
}

\arguments{
  \item{plotsize}{ a numeric vector containing estimates of plot size.}
  \item{CV}{ a numeric vector of experimental coefficient of variation or residual standard error. }
}

\value{
A \code{\link[stats]{nls}} output.
}

\references{
Lessman, K. J. & Atkins, R. E. (1963) Optimum plot size and relative
   efficiency of lattice designs for grain sorghum yield tests.
   \emph{Crop Sci}., 3:477-481.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link[biotools]{optimumplotsize}}
}

\examples{
ps <- c(1, 2, 3, 4, 6, 8, 12)
cv <- c(35.6, 29, 27.1, 25.6, 24.4, 23.3, 21.6)
out <- fitplotsize(plotsize = ps, CV = cv)
predict(out) # fitted.values
plot(cv ~ ps)
curve(coef(out)[1] * x^(-coef(out)[2]), add = TRUE)

# End (not run)
}
