\name{fcTree}
\alias{fcTree}
\title{A class to represent a flowClust tree.}
\usage{
  fcTree(gt)
}
\arguments{
  \item{gt}{a \code{gatingTemplate} object}
}
\description{
  It is a graphNEL used as a container to store priors and
  posteriors for each flowClust gate that can be visualized
  for the purpose of fine-tunning parameters for flowClust
  algorithm

  It adds an extra node data slot "fList"(which is a
  \code{filterList} object) to the \code{gatingTemplate}
}

