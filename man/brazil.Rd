\name{brazil}
\alias{brazil}
\docType{data}

\title{
Brazil Grid
}

\description{
Lat/Long coordinates within Brazil's limits.
}

\usage{data("brazil")}
\format{
  A data frame with 17141 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector (longitude)}
    \item{\code{y}}{a numeric vector (latitude)}
  }
}

\examples{
data(brazil)
plot(brazil, cex = 0.1, col = "gray")
}
\keyword{datasets}
