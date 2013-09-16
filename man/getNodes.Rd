\name{getNodes}
\alias{getNodes}
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
gt <- gatingTemplate(system.file("gatingTemplate/Cytotrol_Tcell.csv",package = "openCyto"))
getNodes(gt)
getNodes(gt, only.names = TURE)
getNodes(gt, '2')
}

