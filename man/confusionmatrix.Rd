\name{confusionmatrix}
\alias{confusionmatrix}

\title{ Confusion Matrix }

\description{
A function to compute the confusion matrix of two classification
vectors. It can be useful for evaluating discriminant analysis or
other classification systems.
}

\usage{
confusionmatrix(obs, predict)
}

\arguments{
  \item{obs}{ a vector containing the observed classes.}
  \item{predict}{ a vector with the same length of \code{obs} containing the predicted classes.}
}

\value{
A square matrix containing the number of objects in each class,
observed (rows) and predicted (columns). Diagonal elements
refers to agreement of \code{obs} and \code{predict}.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link[biotools]{aer}}, \code{\link[MASS]{lda}}
}

\examples{
data(iris)
da <- lda(Species ~ ., data = iris)
pred <- predict(da, dimen = 1)
confusionmatrix(iris$Species, pred$class)

# End (not run)
}
