\name{templateGen}
\alias{templateGen}
\title{generate a partially complete csv template from the existing gating hierarchy}
\usage{
  templateGen(gh)
}
\arguments{
  \item{gh}{a \code{GatingHierarchy} likely parsed from a
  xml workspace}
}
\value{
  a gating template in \code{data.frame} format that
  requires further edition after output to csv
}
\description{
  To ease the process of replicating the existing (usually
  a manual one) gating schemes, this function populate an
  empty gating template with the 'alias', 'pop', 'parent'
  and 'dims' columns that exacted from an
  \code{GatingHierarchy}, and leave the other columns (e.g.
  `gating_method`) blank. So users can make changes to that
  template instead of writing from scratch.
}

