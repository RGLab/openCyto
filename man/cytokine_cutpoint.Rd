\name{cytokine_cutpoint}
\alias{cytokine_cutpoint}
\title{Constructs a cutpoint for a flowFrame by using a derivative of the kernel
density estimate}
\usage{
  cytokine_cutpoint(flow_frame, channel, num_peaks = 1,
    ref_peak = 1,
    method = c("first_deriv", "second_deriv"), tol = 0.01,
    adjust = 1, ...)
}
\arguments{
  \item{flow_frame}{a \code{flowFrame} object}

  \item{channel}{the channel name}

  \item{num_peaks}{the number of peaks expected to see.
  This effectively removes any peaks that are artifacts of
  smoothing}

  \item{ref_peak}{After \code{num_peaks} are found, this
  argument provides the index of the reference population
  from which a gate will be obtained. By default, the peak
  farthest to the left is used.}

  \item{method}{the method used to select the cutpoint. See
  details.}

  \item{tol}{the tolerance value}

  \item{adjust}{the scaling adjustment applied to the
  bandwidth used in the first derivative of the kernel
  density estimate}

  \item{...}{additional arguments passed to
  \code{\link{deriv_density}}}
}
\value{
  the cutpoint along the x-axis
}
\description{
  We determine a gating cutpoint using either the first or
  second derivative of the kernel density estimate (KDE) of
  the \code{channel} specified within the
  \code{flow_frame}.
}
\details{
  By default, we compute the first derivative of the kernel
  density estimate from the channel specified within the
  given \code{flow_frame}. Next, we determine the lowest
  valley from the derivative, which corresponds to the
  density's mode for cytokines. We then contruct a gating
  cutpoint as the value less than the tolerance value
  \code{tol} in magnitude and is also greater than the
  lowest valley.

  Alternatively, if the \code{method} is selected as
  \code{second_deriv}, we select a cutpoint from the second
  derivative of the KDE. Specifically, we choose the
  cutpoint as the largest peak of the second derivative of
  the KDE density which is greater than the reference peak.
}

