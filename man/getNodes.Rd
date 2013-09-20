\name{getNodes}
\alias{getNodes}
\alias{getNodes,fcTree-method}
\alias{getNodes,gatingTemplate-method}
\title{get nodes from \link{gatingTemplate} object}
\arguments{
  \item{x}{\code{gatingTemplate}}

  \item{y}{\code{character} node index. When
  \code{missing}, return all the nodes}

  \item{order}{\code{character} specifying the order of
  nodes. options are "default", "bfs", "dfs", "tsort"}

  \item{only.names}{\code{logical} specifiying whether user
  wants to get the entire \code{gtPopulation} object or
  just the name of the population node}

  \item{x}{\code{fcTree}}

  \item{y}{\code{character} node name}
}
\description{
  get nodes from \link{gatingTemplate} object

  get nodes from \code{fcTree}
}
\examples{
gt <- gatingTemplate(system.file("extdata/template_tcell.csv",package = "openCyto"))
getNodes(gt)[1:2]
getNodes(gt, only.names = TRUE)
getNodes(gt, '2')
}

