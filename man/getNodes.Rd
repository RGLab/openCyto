\name{getNodes}
\alias{getChildren}
\alias{getChildren,gatingTemplate,character-method}
\alias{getGate}
\alias{getGate,fcTree,character-method}
\alias{getGate,gatingTemplate,character-method}
\alias{getNodes}
\alias{getNodes,fcTree-method}
\alias{getNodes,gatingTemplate-method}
\alias{getParent}
\alias{getParent,gatingTemplate,character-method}
\alias{plot,fcFilterList,ANY-method}
\alias{plot,fcTree,character-method}
\alias{plot,filterList,ANY-method}
\alias{plot,gatingTemplate-method}
\alias{plot,gatingTemplate,ANY-method}
\alias{plot,gatingTemplate,character-method}
\alias{plot,gatingTemplate,missing-method}
\alias{ppMethod}
\alias{ppMethod,gatingTemplate,character-method}
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

  \item{x}{\code{fcTree}}

  \item{y}{\code{character} node name}

  \item{obj}{\code{fcTree}}

  \item{y}{\code{character} node name}
}
\description{
  get nodes from \link{gatingTemplate} object

  get children nodes

  get parent nodes

  get gating method from the node

  get preprocessing method from the node

  plot the gating scheme using Rgraphviz

  get nodes from \code{fcTree}

  get gates saved in \code{fcTree}
}
\examples{
gt <- gatingTemplate(system.file("extdata/template_tcell.csv",package = "openCyto"))
getNodes(gt)[1:2]
getNodes(gt, only.names = TRUE)
getNodes(gt, "/nonDebris")
getChildren(gt, "/nonDebris")
getParent(gt, "/nonDebris/singlets")
getGate(gt, "/nonDebris", "/nonDebris/singlets")
ppMethod(gt,  "/nonDebris/singlets",  "/nonDebris/singlets/lymph")
plot(gt) #plot entire tree
plot(gt, "lymph") #only plot the subtree rooted from "lymph"
}

