\name{garlicdist}
\alias{garlicdist}
\docType{data}

\title{ Distances Between Garlic Cultivars }

\description{
The data give the squared generalized Mahalanobis distances between 17 garlic cultivars.
The data are taken from the article published by Silva & Dias (2013).
}

\usage{data(garlicdist)}

\format{
An object of class "dist" based on 17 objects.
}

\source{
Silva, A.R. & Dias, C.T.S. (2013) A cophenetic correlation coefficient for
   Tocher's method. \emph{Pesquisa Agropecuaria Brasileira}, 48:589-596. 
}

\examples{
data(garlicdist)
tocher(garlicdist)

# End (not run)
}

\keyword{datasets}
