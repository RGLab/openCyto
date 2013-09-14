\name{truncate_flowset}
\alias{truncate_flowset}
\title{Removes any observation from each flowFrame object within the flowSet that has
values outside the given range for the specified channels}
\usage{
  truncate_flowset(flow_set, channels, min = NULL,
    max = NULL)
}
\arguments{
  \item{flow_set}{a \code{flowSet} object}

  \item{min}{a numeric vector that sets the lower bounds
  for data filtering}

  \item{max}{a numeric vector that sets the upper bounds
  for data filtering}
}
\value{
  a \code{flowSet} object
}
\description{
  The minimum/maximum values are ignored if \code{NULL}.
}

