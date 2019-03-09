\name{optimumplotsize}
\alias{optimumplotsize}

\title{ Maximum Curvature Point for Optimum Plot Size }

\description{
The Meier & Lessman (1971) method to determine the maximum curvature point
for optimum plot size as a function of the experimental coefficient of variation.
}

\usage{
optimumplotsize(a, b)
}

\arguments{
  \item{a}{ a parameter estimate of the plot size model; see \code{\link{fitplotsize}}. }
  \item{b}{ a parameter estimate of the plot size model; see \code{\link{fitplotsize}}. }
}

\value{
The (approximated) optimum plot size value.
}

\references{
Meier, V. D. & Lessman, K. J. (1971) Estimation of optimum field plot
   shape and size for testing yield in Crambe abyssinia Hochst.
   \emph{Crop Sci}., 11:648-650.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{fitplotsize}}
}

\examples{
ps <- c(1, 2, 3, 4, 6, 8, 12)
cv <- c(35.6, 29, 27.1, 25.6, 24.4, 23.3, 21.6)
out <- fitplotsize(plotsize = ps, CV = cv)
plot(cv ~ ps)
curve(coef(out)[1] * x^(-coef(out)[2]), add = TRUE)
optimumplotsize(a = coef(out)[1], b = coef(out)[2])

# End (not run)
}
