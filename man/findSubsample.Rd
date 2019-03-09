\name{findSubsample}
\alias{findSubsample}

\title{ Finding an Optimized Subsample}

\description{ It allows one to find an optimized (minimized or maximized) numeric subsample according to
	a statistic of interest. For example, it might be of interest to determine a subsample whose
	standard deviation is the lowest among all of those obtained from all possible subsamples of the same size.
}

\usage{
findSubsample(x, size, fun = sd, minimize = TRUE, niter = 10000)
}

\arguments{
  \item{x}{ a numeric vector.}
  \item{size}{ an integer; the size of the subsample.}
  \item{fun}{ an object of class \code{function}; the statistic at which to evaluate the subsample.}
  \item{minimize}{ logical; if TRUE (default) \code{findSubsample} will find a subsample that minimizes \code{stat}.}
  \item{niter}{ an integer indicating the number of iterations, i.e., the number of subsamples to be selected
	(without replacement) from the original sample, \code{x}. The larger is this number, the more optimized is
	the subsample to be found, but this also implies in time-consuming.}
}

%% \details{ }

\value{
A list of
  \item{dataname }{a \code{character}.}
  \item{niter }{the number of iterations.}
  \item{fun }{the objective function.}
  \item{stat }{the achieved statistic for the optimized subsample.}
  \item{criterion }{a \code{character} indicating the type of optimization.}
  \item{subsample }{a numeric vector; the optimized subsample.}
  \item{labels}{ a string containg the labels of the subsample values.}
}

\author{ Anderson Rodrigo da Silva <anderson.agro@hotmail.com> }

\seealso{
\code{\link{sample}}, \code{\link{creategroups}}
}

\examples{
# Example 1
y <- rnorm(40, 5, 2)
findSubsample(x = y, size = 6)

# Example 2
f <- function(x) diff(range(x)) # max(x) - min(x)
findSubsample(x = y, size = 6, fun = f, minimize = FALSE, niter = 20000)

# End (not run)
}