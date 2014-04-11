\name{tailgate}
\alias{cytokine}
\alias{tailgate}
\title{Gates the tail of a density using the derivative of a kernel density estimate
after standardizing and collapsing flowFrames}
\usage{
  tailgate(fr, channel, filter_id = "", num_peaks = 1,
    ref_peak = 1, tol = 0.01, positive = TRUE,
    side = "right", ...)

  cytokine(fr, channel, filter_id = "", num_peaks = 1,
    ref_peak = 1, tol = 0.01, positive = TRUE,
    side = "right", ...)
}
\arguments{
  \item{fr}{a \code{flowFrame} object}

  \item{channel}{the channel from which the cytokine gate
  is constructed}

  \item{filter_id}{the name of the filter}

  \item{num_peaks}{the number of peaks expected to see.
  This effectively removes any peaks that are artifacts of
  smoothing}

  \item{ref_peak}{After \code{num_peaks} are found, this
  argument provides the index of the reference population
  from which a gate will be obtained.}

  \item{tol}{the tolerance value used to construct the
  cytokine gate from the derivative of the kernel density
  estimate}

  \item{positive}{If \code{TRUE}, then the gate consists of
  the entire real line to the right of the cutpoint.
  Otherwise, the gate is the entire real line to the left
  of the cutpoint. (Default: \code{TRUE})}

  \item{side}{On which side of the density do we want to
  gate the tail, the \code{'right'} (default) or
  \code{'left'}?}

  \item{...}{additional arguments.}
}
\value{
  a \code{filterList} containing the gates (cutpoints) for
  each sample
}
\description{
  Gates the tail of a density using the derivative of a
  kernel density estimate after standardizing and
  collapsing flowFrames
}
\examples{
\dontrun{
 gate <- tailgate(fr, Channel = "APC-A") # fr is a flowFrame
}
}

