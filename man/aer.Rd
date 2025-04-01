\name{aer}
\alias{aer}

\title{ Apparent Error Rate }

\description{
A function to calculate the apparent error rate of two classification
vectors, i.e., the proportion of observed cases incorrectly predicted.
It can be useful for evaluating discriminant analysis or other classification systems.
\deqn{aer = \frac{1}{n} \sum_{i=1}^{n} I(y_i \neq \hat{y}_i)}
}

\usage{
aer(obs, predict)
}

\arguments{
  \item{obs}{ a vector containing the observed classes. }
  \item{predict}{ a vector with the same length of \code{obs} containing the predicted classes. }
}

\value{
The apparent error rate, a number between 0 (no agreement) and 1 (thorough agreement).
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link[biotools]{confusionmatrix}}, \code{\link[MASS]{lda}}
}

\examples{
data(iris)
da <- lda(Species ~ ., data = iris)
pred <- predict(da, dimen = 1)
aer(iris$Species, pred$class)

# End (not run)
}
