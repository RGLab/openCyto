\name{registerGatingFunction}
\alias{registerGatingFunction}
\title{Register a gating function with OpenCyto}
\usage{
  registerGatingFunction(fun = NA,
    methodName = "myGatingMethod", dep = NA)
}
\arguments{
  \item{fun}{\code{function} to be registered}

  \item{methodName}{\code{character} name of the gating
  method}
}
\value{
  \code{logical} TRUE if successful and prints a message.
  FALSE otherwise.
}
\description{
  Function registers a new gating method with openCyto so
  that it may be used in the csv template.
}

