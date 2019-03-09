\name{D2.disc}
\alias{D2.disc}
\alias{D2.disc.default}
\alias{print.D2.disc}
\alias{predict.D2.disc}

\title{Discriminant Analysis Based on Mahalanobis Distance}

\description{
A function to perform discriminant analysis based on the squared 
generalized Mahalanobis distance (D2) of the observations to the center of the groups.
}

\usage{
\method{D2.disc}{default}(data, grouping, pooled.cov = NULL)
\method{print}{D2.disc}(x, ...)
\method{predict}{D2.disc}(object, newdata = NULL, ...)
}

\arguments{
  \item{data}{ a numeric \code{data.frame} or \code{matrix} (\emph{n x p}).}
  \item{grouping}{ a vector of length \emph{n} containing the class of each observation (row) in \code{data}.}
  \item{pooled.cov}{a \code{grouping}-pooled covariance matrix (\emph{p x p}). If \code{NULL}
	(default), \code{D2.disc} will automatically compute a pooled covariance matrix.}
  \item{x, object}{ an object of class \code{"D2.disc"}.}
  \item{newdata}{ numeric \code{data.frame} or \code{matrix} of observations to be classified. 
	If \code{NULL} (default), the input data used as argument in \code{D2.disc} will be used.}
  \item{\dots}{ further arguments.}
}

\value{
A list of
  \item{call }{the call which produced the result.}
  \item{data }{ numeric matrix; the input data.}
  \item{D2 }{a matrix containing the Mahalanobis distances between each row of \code{data}
 	and the center of each class of \code{grouping}. In addition, the original
	and the predicted (lowest distance) class are displayed, as well as a
	chacater vector indicating where the misclassification has occured.}
  \item{means }{a matrix containing the vector of means of each class in \code{grouping}.}
  \item{pooled }{the pooled covariance matrix.}
  \item{confusion.matrix }{an object of class \code{\link{confusionmatrix}}.}
}

\references{
Manly, B.F.J. (2004) \emph{Multivariate statistical methods}: a primer. CRC Press. (p. 105-106).

Mahalanobis, P.C. (1936) On the generalized distance in statistics. 
	\emph{Proceedings of The National Institute of Sciences of India}, 12:49-55.
}

\author{Anderson Rodrigo da Silva <anderson.agro@hotmail.com>}

\seealso{
\code{\link{D2.dist}}, \code{\link{confusionmatrix}}, \code{\link{lda}}
}

\examples{
data(iris)
(disc <- D2.disc(iris[, -5], iris[, 5]))
first10 <- iris[1:10, -5]
predict(disc, first10)
predict(disc, iris[, -5])$class

# End (not run)
}
