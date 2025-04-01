\name{raise.matrix}
\alias{raise.matrix}

\title{ Raising a Square Matrix to a Power }

\description{
\code{raise.matrix} raises a square matrix to a power by using
spectral decomposition.
}

\usage{
raise.matrix(x, power = 1)
}

\arguments{
  \item{x}{ a square matrix. }
  \item{power}{ numeric; default is 1.}
}

\value{
An object of class "matrix".
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}


\seealso{
\code{\link[base]{eigen}}, \code{\link[base]{svd}}
}

\examples{
m <- matrix(c(1, -2, -2, 4), 2, 2)
raise.matrix(m)
raise.matrix(m, 2)

# End (not run)
}
