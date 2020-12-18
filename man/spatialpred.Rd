\name{spatialpred}
\alias{spatialpred}

\title{
Spatial Predictions Based on the Circular Variable-Radius Moving Window Method
}

\description{
A heuristic method to perform spatial predictions. The method consists of a local interpolator with stochastic features. 
It allows to build effective detailed maps and to estimate the spatial dependence without any assumptions on the spatial process.
}

\usage{
spatialpred(coords, data, grid)
}

\arguments{
  \item{coords}{ a data frame or numeric matrix containing columns with geographic coordinates}
  \item{data}{ a numeric vector of compatible dimension with \code{coords} containing the response variable data to be predicted at each \code{grid} point }
  \item{grid}{ a data frame or numeric matrix containing columns with geographic coordinates where \code{data} is to be predicted }
}

\value{
A data.frame containing spatial predictions, standard errors, the radius and the number of observations used in each prediction over the grid.
}

\details{
If \code{grid} receives the same input as \code{coords}, \code{spatialpred} will calculate the Percenntual Absolute Mean Error (PAME) of predictions.
}

\references{
Da Silva, A.R., Silva, A.P.A., Tiago-Neto, L.J. (2020) A new local stochastic method for predicting data with spatial heterogeneity. 
	\emph{ACTA SCIENTIARUM-AGRONOMY}, \bold{43}:e49947.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\section{Warning }{
Depending on the dimension of \code{coords} and/or \code{grid}, \code{spatialpred()} can be time demanding.
}

\seealso{
\code{\link{sHe}}
}

\examples{
# data(moco)
# p <- spatialpred(coords = moco[, 1:2], data = rnorm(206), grid = moco[, 1:2]) 
# note: using coords as grid to calculate PAME
# head(p)
# lattice::levelplot(pred ~ Lat*Lon, data = p)

# End (not run)
}