\name{polyfunction_nodes}
\alias{polyfunction_nodes}
\title{Constructs a vector of all the combinations of A & B & C}
\usage{
  polyfunction_nodes(markers)
}
\arguments{
  \item{markers}{character vector of marker names}
}
\value{
  vector containing all combinations of the markers
}
\description{
  The \code{permutations} function is from the
  \code{gregmisc} package on CRAN.
}
\examples{
polyfunction_nodes(c('IFNg', 'IL2', 'TNFa', 'GzB', 'CD57'))
}

