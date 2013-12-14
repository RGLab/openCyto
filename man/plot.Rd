\name{plot}
\alias{plot}
\title{plot a \code{fcFilterList}}
\arguments{
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
  It is usually called by \code{plot} method for
  \code{fcTree} instead of directly by users.

  This provides the priors and posteriors as well as the
  gates for the purpose of debugging flowClust gating
  algorithm
}
\examples{
\dontrun{
 env1<-new.env(parent=emptyenv())
 #gt is a gatingTemplate, gs is a GatingSet
 gating(gt,gs,env1) #the flowClust gating results are stored in env1
 plot(env1$fct,"nonDebris",post=T) #plot the priors as well as posteriors for the "nonDebris" gate
}
}

