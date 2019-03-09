\name{maize}
\alias{maize}
\docType{data}

\title{
Maize Data
}

\description{
Data from and experiment with five maize families carried out in randomized block design, 
with four replications (environments).
}

\usage{data("maize")}

\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{\code{NKPR}}{a numeric vector containing values of Number of Kernels Per cob Row.}
    \item{\code{ED}}{a numeric vector containing values of Ear Diameter (in cm).}
    \item{\code{CD}}{a numeric vector containing values of Cob Diameter (in cm).}
    \item{\code{PH}}{a numeric vector containing values of Plant Heigth (in m).}
    \item{\code{family}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
    \item{\code{env}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
  }
}

\examples{
data(maize)
str(maize)
summary(maize)
}

\keyword{datasets}
