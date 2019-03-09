\name{sHe}
\alias{sHe}

\title{
Spatial Analysis of Gene Diversity
}

\description{
Estimate spatial gene diversity (expected heterozygozity - \emph{He}) through the individual-centred approach by Manel et al. (2007).
\code{sHe()} calculates the unbiased estimate of \emph{He} based on the information of allele frequency obtained from codominant or
dominant markers in individuals within a circular moving windows of known radius over the sampling area.
}

\usage{
sHe(x, coord.cols = 1:2, marker.cols = 3:4, 
	marker.type = c("codominant", "dominant"), 
	grid = NULL, latlong2km = TRUE, radius, nmin = NULL)
}

\arguments{
  \item{x}{ a data frame or numeric matrix containing columns with coordinates of individuals and marker genotyping}
  \item{coord.cols}{ a vector of integer giving the columns of coordinates in \code{x}}
  \item{marker.cols}{ a vector of integer giving the columns of markers in \code{x}}
  \item{marker.type}{ a character; the type of molecular marker}
  \item{grid}{ optional; a two-column matrix containing coordinates over which to predict \emph{He}}
  \item{latlong2km}{ logical; should coordinates be converted from lat/long format into kilometer-grid based?}
  \item{radius}{ the radius of the moving window. It must be in the same format as sampling coordinates}
  \item{nmin}{ optional; a numeric value indicating the minimum number of individuals used to calculate \emph{He}. If is
	the number of individuals in a certain location is less then \code{nmin}, \code{sHe} will consider \emph{He} as zero.}
}

\details{
The unbiased estimate of expected heterogygozity (Nei, 1978) is given by:

\deqn{ He = (1 - \sum_{i=1}^{n} p_{i}^{2}) \frac{2n}{2n - 1} }

where \eqn{p_{i}} is the frequency of the i-th allele per locus considering the \eqn{n} individuals in a certain location.
}

\value{
A list of
  \item{diversity}{ a data frame with the following columns: \emph{coord.x} - the x-axis coordinates of the predicion grid,
	\emph{coord.y} - the y-axis coordinates of the predicion grid, \emph{n} - the number of individuals in a certain points in the grid,
	\emph{MaxDist} - the maximum observed distance among these individuals, \emph{uHe} - the unbiased estimate
	of gene diversity (as expressed above), and \emph{SE} - the standard error of \emph{uHe}. }
  \item{mHe}{ a matrix containing the estimates of \emph{He} for every marker, on each point of the \code{grid}.}
  \item{locations}{ a numeric matrix containing the sampling coordinates, as provides as input.}
}

\references{
da Silva, A.R.; Malafaia, G.; Menezes, I.P.P. (2017) biotools: an R function to predict 
   spatial gene diversity via an individual-based approach. \emph{Genetics and Molecular Research}, 
   \bold{16}: gmr16029655.

Manel, S., Berthoud, F., Bellemain, E., Gaudeul, M., Luikart, G., Swenson, J.E., Waits, L.P., 
	Taberlet, P.; Intrabiodiv Consortium. (2007) A new individual-based spatial approach for 
	identifying genetic discontinuities in natural populations. \emph{Molecular Ecology}, \bold{16}:2031-2043.

Nei, M. (1978) Estimation of average heterozygozity and genetic distance from a small number of individuals. 
	\emph{Genetics}, \bold{89}: 583-590.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>

Ivandilson Pessoa Pinto de Menezes <ivan.menezes@ifgoiano.edu.br>
}

\section{Warning }{
Depending on the dimension of \code{x} and/or \code{grid}, \code{sHe()} can be time demanding.
}

\seealso{
\code{\link{levelplot}}
}

\examples{
data(moco)
data(brazil)

# check points
plot(brazil, cex = 0.1, col = "gray")
points(Lat ~ Lon, data = moco, col = "blue", pch = 20)

# using a retangular grid (not passed as input!)
# ex <- sHe(x = moco, coord.cols = 1:2,
#	marker.cols = 3:20, marker.type = "codominant",
#	grid = NULL, radius = 150)
#ex
# plot(ex, xlab = "Lon", ylab = "Lat")

# A FANCIER PLOT...
# using Brazil's coordinates as prediction grid
# ex2 <- sHe(x = moco, coord.cols = 1:2,
#	marker.cols = 3:20, marker.type = "codominant",
#	grid = brazil, radius = 150)
# ex2
#
# library(maps)
# borders <- data.frame(x = map("world", "brazil")$x, 
#	y = map("world", "brazil")$y)
#
# library(latticeExtra)
# plot(ex2, xlab = "Lon", ylab = "Lat", 
#	xlim = c(-75, -30), ylim = c(-35, 10), aspect = "iso") +
#   latticeExtra::as.layer(xyplot(y ~ x, data = borders, type = "l")) +
#   latticeExtra::as.layer(xyplot(Lat ~ Lon, data = moco))

# End (not run)
}