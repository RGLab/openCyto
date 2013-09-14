\name{quadGate2rectangleGates}
\alias{quadGate2rectangleGates}
\title{Extracts the quadrants of a quadGate as a list of rectangleGates}
\usage{
  quadGate2rectangleGates(quad_gate, markers, channels,
    quadrants = 1:4)
}
\arguments{
  \item{quad_gate}{a \code{quadGate} object}

  \item{markers}{character vector of the marker names for
  the x- and y-axes}

  \item{channels}{character vector of the channel names for
  the x- and y-axes}

  \item{quadrants}{a vector indicating the quadrants to
  extract}
}
\value{
  a \code{filters} object containing a list of the
  rectangle gates
}
\description{
  The quadrants are numbered in a clockwise manner with the
  top-left quadrant numbered 1, the top-right quadrant
  numbered 2, and so on.
}

