\name{distClust}
\alias{distClust}

\title{ Cluster Distance Matrix }

\description{
Function to compute a matrix of average distances within and between clusters.
}

\usage{
distClust(d, nobj.cluster, id.cluster)
}

\arguments{
  \item{d}{ an object of class "dist" containing the distances between objects.}
  \item{nobj.cluster}{ a numeric vector containing the numbers of objects per cluster. }
  \item{id.cluster}{ a numeric vector for identification of the objects per cluster. }
}

\value{
A squared matrix containing distances within (diagonal) and between
(off-diagonal) clusters.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}


\seealso{
\code{\link[biotools]{tocher}}, \code{\link[stats]{dist}}
}
