\name{deriv_density}
\alias{deriv_density}
\title{Constructs the derivative specified of the kernel density estimate of a
numeric vector}
\usage{
  deriv_density(x, deriv = 1, bandwidth = NULL, adjust = 1,
    num_points = 10000, ...)
}
\arguments{
  \item{x}{numeric vector}

  \item{deriv}{a numeric value specifying which derivative
  should be calculated. By default, the first derivative is
  computed.}

  \item{bandwidth}{the bandwidth to use in the kernel
  density estimate. If \code{NULL} (default), the bandwidth
  is estimated using the plug-in estimate from
  \code{\link[ks]{hpi}}.}

  \item{adjust}{a numeric weight on the automatic
  bandwidth, analogous to the \code{adjust} parameter in
  \code{\link{density}}}

  \item{num_points}{the length of the derivative of the
  kernel density estimate}

  \item{...}{additional arguments passed to
  \code{\link[ks:drvkde]{drvkde}}}
}
\value{
  list containing the derivative of the kernel density
  estimate
}
\description{
  The derivative is computed with
  \code{\link[ks:drvkde]{drvkde}}.
}
\details{
  For guidance on selecting the bandwidth, see this
  CrossValidated post: \url{http://bit.ly/12LkJWz}
}

