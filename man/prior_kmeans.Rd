\name{.prior_kmeans}
\alias{.prior_kmeans}
\title{Elicits data-driven priors from a flowSet object for specified channels using
the K-Means clustering algorithm}
\usage{
  .prior_kmeans(flow_set, channels, K, nu0 = 4, w0 = 10,
    nstart = 10, pct = 0.1, min = NULL, max = NULL, ...)
}
\arguments{
  \item{flow_set}{a \code{flowSet} object}

  \item{channels}{a character vector containing the
  channels in the \code{flowSet} from which we elicit the
  prior parameters for the Student's t mixture}

  \item{K}{the number of mixture components to identify}

  \item{nu0}{prior degrees of freedom of the Student's t
  mixture components.}

  \item{w0}{the number of prior pseudocounts of the
  Student's t mixture components.}

  \item{nstart}{number of random starts used by
  \code{kmeans} algorithm}

  \item{pct}{percentage of randomly selected cells in each
  \code{flowFrame} that is used to elicit the prior
  parameters. The value should must be greater than 0 and
  less than or equal to 1.}

  \item{min}{a numeric vector that sets the lower bounds
  for data filtering. If \code{NULL} (default), no
  truncation is applied.}

  \item{max}{a numeric vector that sets the upper bounds
  for data filtering. If \code{NULL} (default), no
  truncation is applied.}

  \item{...}{Additional arguments passed to \code{kmeans}}
}
\value{
  list of \code{flowClust} prior parameters
}
\description{
  We elicit data-driven prior parameters from a
  \code{flowSet} object for specified channels. For each
  sample in the \code{flowSet} object, we apply
  \code{kmeans} to obtain \code{K} clusters. From each
  cluster, we determine its centroid and the sample
  covariance matrix. We then aggregate these two sample
  moments across all samples for each cluster.
}
\details{
  Because the cluster labels returned from \code{kmeans}
  are arbitrary, we align the clusters based on the
  centroids that are closest to a randomly selected
  reference sample. We apply the Hungarian algorithm
  implemented using the \code{solve_LSAP} function from the
  \code{clue} package to assist with the alignment.

  If each frame within \code{flow_set} has a large number
  of cells, the computational costs of \code{kmeans} can be
  a burden. We provide the option to randomly select
  \code{pct}, a percentage of the cells from each flow
  frame to which \code{kmeans} is applied.
}

