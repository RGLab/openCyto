\name{standardize_flowset}
\alias{standardize_flowset}
\title{Standardizes a channel within a \code{flowSet} object using the mode of the
kernel density estimate and the Huber estimator of the standard deviation}
\usage{
  standardize_flowset(flow_set, channel = "FSC-A")
}
\arguments{
  \item{flow_set}{a \code{flowSet} object}

  \item{channel}{the channel to standardize}
}
\value{
  list containing the transformed \code{flowSet} object
  along with the \code{transformation} list, where each
  element contains the transformation parameters for each
  \code{flowFrame}
}
\description{
  Standardizes a channel within a \code{flowSet} object
  using the mode of the kernel density estimate and the
  Huber estimator of the standard deviation
}

