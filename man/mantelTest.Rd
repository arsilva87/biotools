\name{mantelTest}
\alias{mantelTest}

\title{ Mantel's Permutation Test }

\description{
Mantel's permutation test based on Pearson's correlation coefficient
to evaluate the association between two distance square matrices.
}

\usage{
mantelTest(m1, m2, nperm = 999, alternative = "greater",
	graph = TRUE, main = "Mantel's test", xlab = "Correlation", ...)
}

\arguments{
  \item{m1}{ an object of class "matrix" or "dist", containing distances among \emph{n} individuals.}
  \item{m2}{ an object of class "matrix" or "dist", containing distances among \emph{n} individuals.}
  \item{nperm}{ the number of matrix permutations.}
  \item{alternative}{ a character specifying the alternative hypothesis. It must be one of "greater" (default),
	"two.sided" or "less".}
  \item{graph}{ logical; if TRUE (default), the empirical distribution is plotted.}
  \item{main}{ opitional; a character describing the title of the graphic.}
  \item{xlab}{opitional; a character describing the \emph{x}-axis label.}
  \item{\dots}{ further graphical arguments. See \code{\link[graphics]{par}}.}
}

\value{
A list of
  \item{correlation }{ numeric; the observed Pearson's correlation between \code{m1} and \code{m2}.}
  \item{p.value }{ numeric; the empirical p-value of the permutation test.}
  \item{alternative }{ character; the alternative hypothesis used to compute \code{p.value}.}
  \item{nullcor }{ numeric vector containing randomized values of correlation, i.e., under the null hypothesis
	that the true correlation is equal to zero.}
}

\references{
Mantel, N. (1967). The detection of disease clustering and a generalized regression approach.
	\emph{Cancer Research}, 27:209--220.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link[biotools]{mantelPower}}
}

\examples{
# Distances between garlic cultivars
data(garlicdist)
garlicdist

# Tocher's clustering
garlic <- tocher(garlicdist)
garlic

# Cophenetic distances
coph <- cophenetic(garlic)
coph

# Mantel's test
mantelTest(garlicdist, coph,
	xlim = c(-1, 1))

# End (Not run)
}
