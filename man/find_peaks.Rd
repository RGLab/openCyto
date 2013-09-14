\name{find_peaks}
\alias{find_peaks}
\title{Finds the local maxima (peaks) in the given vector after smoothing the data
with a kernel density estimator.}
\usage{
  find_peaks(x, y = NULL, num_peaks = NULL, adjust = 2,
    ...)
}
\arguments{
  \item{x}{numeric vector}

  \item{y}{numeric vector. If given, it is treated as the
  density values for \code{x}. The length of \code{y} must
  equal that of \code{x}.}

  \item{num_peaks}{the number of peaks to find. By default,
  all peaks are returned. See details.}

  \item{adjust}{the bandwidth to use in the kernel density
  estimation. See \code{\link{density}} for more
  information.}

  \item{...}{additional arguments passed to the
  \code{\link{density}} function}
}
\value{
  the values where the peaks are attained. The peaks are
  sorted in descending order based on the density heights.
}
\description{
  First, we smooth the data using kernel density estimation
  (KDE) with the \code{\link{density}} function. Then, we
  find every local maxima such that the density is concave
  (downward).
}
\details{
  The \code{num_peaks} argument returns only the largest
  peaks. If \code{num_peaks} is greater than the number of
  peaks found, then all the peaks are returned, and a
  warning is issued.

  Effectively, we find the local maxima with a discrete
  analogue to a second derivative applied to the KDE. For
  details, see this StackOverflow post:
  \url{http://bit.ly/Zbl7LV}.
}
\examples{
library(flowClust)
set.seed(42)
# 2 peaks with a minor hump
y <- SimulateMixture(10000, c(.5, .3, .2), c(2, 5, 7), c(1, 1, 1), nu = 10)
plot(density(y))
peaks <- find_peaks(y)
abline(v = peaks, col = "red")
}

