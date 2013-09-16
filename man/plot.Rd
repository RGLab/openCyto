\name{plot}
\alias{plot}
\title{plot the gating scheme}
\arguments{
  \item{x}{\code{gatingTemplate} object}

  \item{y}{either \code{character} specifying the root node
  which can be used to visualize only the subgraph or
  \code{missing} which display the entire gating scheme}

  \item{...}{other arguments

  graphAttr, nodeAttr: graph rendering attributes passed to
  \link[Rgraphviz:renderGraph]{renderGraph} showRef
  \code{logical}: whether to display the reference gates.
  Sometime it maybe helpful to hide all those reference
  gates which are not the cell population of interest and
  used primarily for generating other population nodes.}

  \item{x}{\code{fcFilterList}}

  \item{y}{\code{character} channel name}

  \item{samples}{\code{character} a vector of sample names
  to be plotted}

  \item{posteriors}{\code{logical} indicating whether
  posteriors should be plotted}

  \item{xlim,ylim}{scale settings for x,y axises}

  \item{lwd}{line width}

  \item{breaks}{passed to \link{hist}}

  \item{data}{\code{GatingSet} object}

  \item{node}{\code{character} population name associated
  with the \code{fcFilterList}}

  \item{...}{other arguments passed to base \code{plot}}

  \item{x}{\code{fcTree}}

  \item{y}{\code{character} node name in the \code{fcTree}}

  \item{channel}{\code{character} specifying the channel.}

  \item{data}{\code{GatingSet} that the \code{fcTree} is
  associated with}

  \item{...}{other arguments}
}
\description{
  plot the gating scheme using Rgraphviz

  It is usually called by \code{plot} method for
  \code{fcTree} instead of directly by users.

  This provides the priors and posteriors as well as the
  gates for the purpose of debugging flowClust gating
  algorithm
}

